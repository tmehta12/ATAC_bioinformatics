#!/usr/bin/env python

# TSS enrich wrapper
# Adapted from 'ENCODE DCC TSS enrich wrapper' by Daniel Kim, Jin Lee (leepc12@gmail.com)
# Author: Tarang K. Mehta
# USAGE: python3 ATAC_Bioinf_pipeline_v2b_part5bD.py 'FINAL_BAM' 'TSS' 'OUTPUT_PREFIX' 'read_len' 'CHROMSIZES'

import sys
import os
from matplotlib import mlab
import logging
import pybedtools
# import metaseq
from matplotlib import pyplot as plt
import numpy as np
import plotutils
import subprocess

# These require the relevant python definition files in the same folder
from array_helpers import _array, _array_parallel, _local_coverage, \
    _local_count, _count_array, _count_array_parallel
import filetype_adapters
import helpers
from helpers import rebin
from bx.bbi.bigwig_file import BigWigFile


# os.environ['QT_QPA_PLATFORM']='offscreen'

####################################################################################################
## Define the input files

FINAL_BAM = sys.argv[1]
TSS = sys.argv[2]
OUTPUT_PREFIX = sys.argv[3]
read_len = sys.argv[4] # this will take the read lengths as a single integer as calculated in main script
CHROMSIZES = sys.argv[5]

# store the read length as an integer
for line in open(read_len):
   if line.strip():           # line contains eol character(s)
       read_length = int(line)          # assuming single integer on each line

####################################################################################################
## Define classes

class BaseSignal(object):
    """
    Base class to represent objects from which genomic signal can be
    calculated/extracted.

    `__getitem__` uses the underlying adapter the instance was created with
    (e.g., :class:`metaseq.filetype_adapters.BamAdapter` for
    a :class:`BamSignal` object).
    """
    def __init__(self, fn):
        self.fn = fn

    def array(self, features, processes=None, chunksize=1, ragged=False,
              **kwargs):
        """
        Creates an MxN NumPy array of genomic signal for the region defined by
        each feature in `features`, where M=len(features) and N=(bins or
        feature length)

        Parameters
        ----------
        features : iterable of interval-like objects
            An iterable of interval-like objects; see docstring for
            `local_coverage` method for more details.

        processes : int or None
            If not None, then create the array in parallel, giving each process
            chunks of length `chunksize` to work on.

        chunksize : int
            `features` will be split into `chunksize` pieces, and each piece
            will be given to a different process. The optimum value is
            dependent on the size of the features and the underlying data set,
            but `chunksize=100` is a good place to start.

        ragged : bool
            If False (default), then return a 2-D NumPy array.  This requires
            all rows to have the same number of columns, which you get when
            supplying `bins` or if all features are of uniform length.  If
            True, then return a list of 1-D NumPy arrays

        Notes
        -----
        Additional keyword args are passed to local_coverage() which performs
        the work for each feature; see that method for more details.
        """
        if processes is not None:
            arrays = _array_parallel(
                self.fn, self.__class__, features, processes=processes,
                chunksize=chunksize, **kwargs)
        else:
            arrays = _array(self.fn, self.__class__, features, **kwargs)
        if not ragged:
            stacked_arrays = np.row_stack(arrays)
            del arrays
            return stacked_arrays
        else:
            return arrays

    def local_coverage(self, features, *args, **kwargs):
        processes = kwargs.pop('processes', None)
        if not processes:
            return _local_coverage(self.adapter, features, *args, **kwargs)

        if isinstance(features, (list, tuple)):
            raise ValueError(
                "only single features are supported for parallel "
                "local_coverage")

        # we don't want to have self.array do the binning
        bins = kwargs.pop('bins', None)

        # since if we got here processes is not None, then this will trigger
        # a parallel array creation
        features = helpers.tointerval(features)
        x = np.arange(features.start, features.stop)
        features = list(helpers.split_feature(features, processes))
        ys = self.array(
            features, *args, bins=None, processes=processes, ragged=True,
            **kwargs)
        # now we ravel() and re-bin
        y = np.column_stack(ys).ravel()
        if bins:
            xi, yi = rebin(x, y, bins)
            del x, y
            return xi, yi
        return x, y

    local_coverage.__doc__ = _local_coverage.__doc__

class IntervalSignal(BaseSignal):
    def __init__(self, fn):
        """
        Abstract class for bed, BAM and bigBed files.
        """
        BaseSignal.__init__(self, fn)

    def local_count(self, *args, **kwargs):
        return _local_count(self.adapter, *args, **kwargs)

    local_count.__doc__ = _local_count.__doc__

    def count_array(self, features, processes=None, chunksize=1,  **kwargs):
        if processes is not None:
            arrays = _count_array_parallel(
                self.fn, self.__class__, features, processes=processes,
                chunksize=chunksize, **kwargs)
        else:
            arrays = _count_array(self.fn, self.__class__, features, **kwargs)
        return np.concatenate(arrays)

class BigWigSignal(BaseSignal):
    def __init__(self, fn):
        """
        Class for operating on bigWig files
        """
        super(BigWigSignal, self).__init__(fn)
        self.adapter = filetype_adapters.BigWigAdapter(fn)

class BamSignal(IntervalSignal):
    def __init__(self, fn):
        """
        Class for operating on BAM files.
        """
        BaseSignal.__init__(self, fn)
        self._readcount = None
        self.adapter = filetype_adapters.BamAdapter(self.fn)

    def genome(self):
        """
        "genome" dictionary ready for pybedtools, based on the BAM header.
        """
        # This gets the underlying pysam Samfile object
        f = self.adapter.fileobj
        d = {}
        for ref, length in zip(f.references, f.lengths):
            d[ref] = (0, length)
        return d

    def mapped_read_count(self, force=False):
        """
        Counts total reads in a BAM file.

        If a file self.bam + '.scale' exists, then just read the first line of
        that file that doesn't start with a "#".  If such a file doesn't exist,
        then it will be created with the number of reads as the first and only
        line in the file.

        The result is also stored in self._readcount so that the time-consuming
        part only runs once; use force=True to force re-count.

        Parameters
        ----------
        force : bool
            If True, then force a re-count; otherwise use cached data if
            available.
        """
        # Already run?
        if self._readcount and not force:
            return self._readcount

        if os.path.exists(self.fn + '.mmr') and not force:
            for line in open(self.fn + '.mmr'):
                if line.startswith('#'):
                    continue
                self._readcount = float(line.strip())
                return self._readcount

        cmds = ['samtools',
                'view',
                '-c',
                '-F', '0x4',
                self.fn]
        p = subprocess.Popen(
            cmds, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = p.communicate()
        if stderr:
            sys.stderr.write('samtools says: %s' % stderr)
            return None
        mapped_reads = int(stdout)

        # write to file so the next time you need the lib size you can access
        # it quickly
        if not os.path.exists(self.fn + '.mmr'):
            fout = open(self.fn + '.mmr', 'w')
            fout.write(str(mapped_reads) + '\n')
            fout.close()

        self._readcount = mapped_reads
        return self._readcount

class BigBedSignal(IntervalSignal):
    def __init__(self, fn):
        """
        Class for operating on bigBed files.
        """
        IntervalSignal.__init__(self, fn)
        self.adapter = filetype_adapters.BigBedAdapter(fn)


class BedSignal(IntervalSignal):
    def __init__(self, fn):
        """
        Class for operating on BED files.
        """
        IntervalSignal.__init__(self, fn)
        self.adapter = filetype_adapters.BedAdapter(fn)

####################################################################################################
## Define registry for extension identification

_registry = {
    'bam': BamSignal,
    'bed': BedSignal,
    'gff': BedSignal,
    'gtf': BedSignal,
    'vcf': BedSignal,
    'bigwig': BigWigSignal,
    'bigbed': BigBedSignal,
}

####################################################################################################
## Define internal functions for file format identifiers

def genomic_signal(fn, kind):
    """
    Factory function that makes the right class for the file format.

    Typically you'll only need this function to create a new genomic signal
    object.

    :param fn: Filename
    :param kind:
        String.  Format of the file; see
        metaseq.genomic_signal._registry.keys()
    """
    try:
        klass = _registry[kind.lower()]
    except KeyError:
        raise ValueError(
            'No support for %s format, choices are %s'
            % (kind, _registry.keys()))
    m = klass(fn)
    m.kind = kind
    return m

def data_dir():
    """
    Returns the data directory that contains example files for tests and
    documentation.
    """
    return os.path.join(os.path.dirname(__file__), 'test', 'data')


def example_filename(fn):
    """
    Return a bed file from the pybedtools examples directory.  Use
    :func:`list_example_files` to see a list of files that are included.
    """
    fn = os.path.join(data_dir(), fn)
    if not os.path.exists(fn):
        raise ValueError("%s does not exist" % fn)
    return fn


####################################################################################################
## Define functions for doing the TSS enrichment plot

def make_tss_plot(bam_file, tss, prefix, chromsizes, read_len, bins=400, bp_edge=2000,
                  processes=8, greenleaf_norm=True):
    '''
    Take bootstraps, generate tss plots, and get a mean and
    standard deviation on the plot. Produces 2 plots. One is the
    aggregation plot alone, while the other also shows the signal
    at each TSS ordered by strength.
    '''
    logging.info('Generating tss plot...')
    tss_plot_file = '{0}_tss-enrich.png'.format(prefix)
    tss_plot_large_file = '{0}_large_tss-enrich.png'.format(prefix)

    # Load the TSS file
    tss = pybedtools.BedTool(tss)
    tss_ext = tss.slop(b=bp_edge, g=chromsizes)

    # Load the bam file
    bam = genomic_signal(bam_file, 'bam') # Need to shift reads and just get ends, just load bed file?
    # bam = metaseq.genomic_signal(bam_file, 'bam') # Need to shift reads and just get ends, just load bed file?
    bam_array = bam.array(tss_ext, bins=bins, shift_width = -read_len/2, # Shift to center the read on the cut site
                          processes=processes, stranded=True)

    # Actually first build an "ends" file
    #get_ends = '''zcat {0} | awk -F '\t' 'BEGIN {{OFS="\t"}} {{if ($6 == "-") {{$2=$3-1; print}} else {{$3=$2+1; print}} }}' | gzip -c > {1}_ends.bed.gz'''.format(bed_file, prefix)
    #print(get_ends)
    #os.system(get_ends)

    #bed_reads = metaseq.genomic_signal('{0}_ends.bed.gz'.format(prefix), 'bed')
    #bam_array = bed_reads.array(tss_ext, bins=bins,
    #                      processes=processes, stranded=True)

    # Normalization (Greenleaf style): Find the avg height
    # at the end bins and take fold change over that
    if greenleaf_norm:
        # Use enough bins to cover 100 bp on either end
        num_edge_bins = int(100/(2*bp_edge/bins))
        bin_means = bam_array.mean(axis=0)
        avg_noise = (sum(bin_means[:num_edge_bins]) +
                     sum(bin_means[-num_edge_bins:]))/(2*num_edge_bins)
        bam_array /= avg_noise
    else:
        bam_array /= bam.mapped_read_count() / 1e6

    # Generate a line plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    x = np.linspace(-bp_edge, bp_edge, bins)

    ax.plot(x, bam_array.mean(axis=0), color='r', label='Mean')
    ax.axvline(0, linestyle=':', color='k')

    # Note the middle high point (TSS)
    tss_point_val = max(bam_array.mean(axis=0))

    ax.set_xlabel('Distance from TSS (bp)')
    if greenleaf_norm:
        ax.set_ylabel('TSS Enrichment')
    else:
        ax.set_ylabel('Average read coverage (per million mapped reads)')
    ax.legend(loc='best')

    fig.savefig(tss_plot_file)

    # Print a more complicated plot with lots of info
    # Find a safe upper percentile - we can't use X if the Xth percentile is 0
    upper_prct = 99
    if mlab.prctile(bam_array.ravel(), upper_prct) == 0.0:
        upper_prct = 100.0

    plt.rcParams['font.size'] = 8
    fig = plotutils.imshow(bam_array,
                                   x=x,
                                   figsize=(5, 10),
                                   vmin=5, vmax=upper_prct, percentile=True,
                                   line_kwargs=dict(color='k', label='All'),
                                   fill_kwargs=dict(color='k', alpha=0.3),
                                   sort_by=bam_array.mean(axis=1))
    # fig = metaseq.plotutils.imshow(bam_array,
    #                                x=x,
    #                                figsize=(5, 10),
    #                                vmin=5, vmax=upper_prct, percentile=True,
    #                                line_kwargs=dict(color='k', label='All'),
    #                                fill_kwargs=dict(color='k', alpha=0.3),
    #                                sort_by=bam_array.mean(axis=1))

    # And save the file
    fig.savefig(tss_plot_large_file)
    return tss_plot_file, tss_plot_large_file, tss_point_val

####################################################################################################
## Run and generate the plots

tss_plot_file, tss_plot_large_file, tss_point_val = make_tss_plot(FINAL_BAM, # Use final to avoid duplicates
                                                                  TSS,
                                                                  OUTPUT_PREFIX,
                                                                  CHROMSIZES,
                                                                  read_length)

# if __name__=='__main__':
#     tss_plot_file, tss_plot_large_file, tss_point_val = make_tss_plot(FINAL_BAM, # Use final to avoid duplicates
#                                                                       TSS,
#                                                                       OUTPUT_PREFIX,
#                                                                       CHROMSIZES,
#                                                                       read_length)
#     print '{0}  {1}'.format(FINAL_BAM, str(tss_point_val))
