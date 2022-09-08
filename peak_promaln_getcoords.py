#!/usr/bin/python
#=====================================================================================================================================================

# Using a bed file, this script will calculate the coordinates of those regions from an alignment fasta
# use case: calculate the start and end of an ATAC-seq peak, plus the peak summit according to the alignment coordinates (accounting for the alignment gaps of '-')
# you need to split the bed file of peaks into unique orthogroups and name the files by orthogroup.
# this allows to run as an array - will work on each orthogroup individually by only loading it's associated FASTA file
# best to only run for 1-to-1 orthogroups

#=====================================================================================================================================================

print('\n===================================================================================\n')
print('peak_promaln_getcoords.py')
import sys

if len(sys.argv) < 3:
    print('\nNot enough arguments entered.\nUsage: python3 peak_promaln_getcoords.py multilineFASTA.fa peakBEDfile outBEDfile\n')
    print('\n===================================================================================\n')
    sys.exit()


print("\nImporting...\n\n")
import os
import csv
import Bio
import time
import numpy as np
import itertools
import subprocess
import collections
import re

from Bio import SeqIO
from Bio import AlignIO
from collections import defaultdict
from collections import Counter


# 1. read in the multiline FASTA and store as dict
record_dict = SeqIO.to_dict(SeqIO.parse(sys.argv[1], "fasta"))

# record_dict = SeqIO.to_dict(SeqIO.parse('OMA00703_aln.fa', "fasta"))

if record_dict:
    print('Multiline FASTA file' + str(sys.argv[1]) + 'dict is not empty! Continuing ...')
else:
    print('Multiline FASTA file' + str(sys.argv[1]) + 'dict IS empty! Aborting ...')
    sys.exit()

print("\nmultiline FASTA recorded to dict record_dict...\n\n")

# record_dict["ENSMZEG00005017372"] # access a record using key
# str(record_dict["ENSMZEG00005017372"].seq) # access sequence value

# example processing
# for m in record_dict:
#     sequence_key1 = (record_dict[m].id) # access the key
#     print(sequence_key1)

# 2. open the bed file and for each line, strip by newline, split by tab and assign to bed_out list
with open(sys.argv[2], 'r') as bed_file:
    bed_out = [line.strip('\n').split('\t') for line in bed_file]

# with open('OMA00703.promap.bed', 'r') as bed_file:
#     bed_out = [line.strip('\n').split('\t') for line in bed_file]

# for line in bed_out:
#     print(line)

# print(list(bed_out)[:1])

# example processing:
# for x in bed_out:
#     if str(x[0]) == str(x[0]):
#         print(x[0],';',x[1])

print("\nbed file assigned to bed_out list...\n\n")

# 3. calculate the start and end of an ATAC-seq peak, plus the peak summit according to the alignment coordinates (accounting for the alignment gaps of '-')

# i. match the ensembl IDs (FASTA key and col0 in BED)
# ii. determine how many alignment gaps exist in the FASTA seq - need to ignore gaps and count along nucleotides (using the relevant coords below) and get new position with gaps
    # before the peak_start
    # before the peak_end
    # before the peak_summit

seq_range = [] # list of seq range
gap_range = [] # list of gap range
seqgapntpos_dict = {} # dict of gaps and nt's in seq
outbed = [] # out bed file as list

for m in record_dict:
    ensID_fasta = (record_dict[m].id) # access the key
    # print(sequence_key1)
    for row in bed_out:
        ensID_bed = row[0]
        peak_start = row[1]
        peak_end = row[2]
        peakID = row[3]
        score = row[4]
        strand = row[5]
        MACS2_signalValue = row[6]
        MACS2_pval = row[7]
        MACS2_qval = row[8]
        peak_summit = row[9]
        IDR = row[10]
        orthogroup = row[11]
        ens_gs1 = row[12]
        ens_gs2 = row[13]
        prom_chr = row[14]
        prom_start = row[15]
        prom_end = row[16]
        peaklength = row[17]
        tss_coord = row[18]
        peak_genchr = row[19]
        peak_genstart = row[20]
        peak_genend = row[21]
        if ensID_fasta == ensID_bed: # i. match the ensembl IDs (FASTA dict key and col0 in BED)
            # print(ensID_fasta,' : ',record_dict[m].seq) # prints the ID and sequence
            # record_dict[m].seq[0] # returns the first base
            # record_dict[m].seq.count('-') # counts the total number of gaps '-'
            # record_dict[m].seq[0:50].count('-') # counts the total number of gaps '-' between 0 and 50
            # re.findall("[a-zA-Z_]+", str(record_dict[m].seq[0:50])) # prints all a-z runs between 0 and 50
            # print(len(str(record_dict[m].seq))) # print's the length - note that 788 means it is 0-787
            gaps = list(re.finditer('-+', str(record_dict[m].seq))) # this lists all the consecutive gap regions
            nts = list(re.finditer(r'[A-Za-z]+', str(record_dict[m].seq))) # this lists all the consecutive nt regions
            # print('Number of gaps =', len(gaps)) # prints total number of gaps
            # print(gaps)
            # print(nts)
            for i in range(len(str(record_dict[m].seq))):
                seq_range.append(i) # add the seq range (all nt coords) to a list
            # for nts_region_number, nts_match in enumerate(nts, 1):
            #     print('nt region:',nts_region_number,';',nts_match.start(),':',nts_match.end()-1)
            for gap_region_number, gap_match in enumerate(gaps, 1):
                # print('gap region #',gap_region_number,';',gap_match.start(),':',gap_match.end()-1)
                for j in range(gap_match.start(),gap_match.end()):
                    gap_range.append(j) # add the gap coords to a list
                # print(gap_match.span()) # note that .span, .start, and .end are all part of the re package
                # print('Index Position of Gap region {} = {} to {}'.format(
                #         region_number,
                #         match.start(),
                #         match.end() - 1)) # prints the positions of each gap
                # print('Length of Gap region {} = {}'.format(
                #         region_number,
                #         match.end() - match.start())) # prints the length of each gap
                # for j in range(gap_match.start(),gap_match.end()):
                #     print(j,"gap") # this will print each position in the sequence that is a gap
            seq_range2 = set(seq_range)
            gap_range2 = set(gap_range)
            setgap = set(seq_range2 | gap_range2)
            for elem in sorted(setgap):
                if elem in seq_range2:
                    if elem in gap_range2:
                        seqgapntpos_dict[elem]='gap'
                    else:
                        seqgapntpos_dict[elem]='nt'
                else:
                    seqgapntpos_dict['gap']=elem # mark each position in the seq as either a nt or gap to find the new peak positions
            seqntpos_dict = {k:v for k, v in seqgapntpos_dict.items() if v == 'nt'} # just add the nt coords to a dict to get the nth position of peak in alignment
            if int(peak_end) > int(len(seqntpos_dict)): # the peak can extend beyond promoter - check if the peak end is more than the length of the number of nt's in the alignment; if so, store peak end as the end of the alignment (even if a gap)
                print('peak extends beyond promoter - peak end stored as the end of the alignment which will include gaps')
                peak_alnend = str(len(record_dict[m].seq) - 1) # this does NOT need a new position calculation as it includes gaps
                peak_alnstart = str(list(seqntpos_dict)[int(peak_start)]) # get the new peak_start position in alignment (taking gaps into consideration)
            else:
                print('peak does not extend beyond promoter - new peak end will be calculated')
                peak_alnend = str(list(seqntpos_dict)[int(peak_end) - 1]) # this is a new peak end position calculation that includes gaps as the peak ends in the gene promoter alignment, and doesn't extend beyond
                peak_alnstart = str(list(seqntpos_dict)[int(peak_start)]) # get the new peak_start position in alignment (taking gaps into consideration)
            if int(peak_summit) > int(len(seqntpos_dict)): # the peak summit can extend beyond promoter - check if the summit is more than the length of the number of nt's in the alignment; if so, store the summit as +n
                print('peak summit extends beyond promoter - peak summit will be stored as a + figure')
                peak_alnsummit = (str('+') + str(int(peak_summit) - int(len(seqntpos_dict)))) # get the new peak_start position in alignment (taking gaps into consideration)
            else:
                print('peak summit does not extend beyond promoter - new peak summit in alignment will be calculated')
                peak_alnsummit = str(list(seqntpos_dict)[int(peak_summit) - 1]) # get the new peak_start position in alignment (taking gaps into consideration)
            outbed.extend([[ensID_bed,peak_alnstart,peak_alnend,peakID,score,strand,MACS2_signalValue,MACS2_pval,MACS2_qval,peak_alnsummit,IDR,orthogroup,ens_gs1,ens_gs2,prom_chr,prom_start,prom_end,peaklength,tss_coord,peak_genchr,peak_genstart,peak_genend,peak_start,peak_end,peak_summit]])

with open(sys.argv[3], 'w') as file: # open a file for writing
    file.writelines('\t'.join(i) + '\n' for i in outbed) # join each element as tab-separated and trailing newlines for each line in the list
file.close()

print('\nOut bed file written with file name:\t' + str(sys.argv[3]))

## outfile format
# 1. ensID_bed
# 2. peak_alnstart
# 3. peak_alnend
# 4. peakID
# 5. score
# 6. strand
# 7. MACS2_signalValue
# 8. MACS2_pval
# 9. MACS2_qval
# 10. peak_alnsummit
# 11. IDR
# 12. orthogroup
# 13. ens_gs1
# 14. ens_gs2
# 15. prom_chr
# 16. prom_start
# 17. prom_end
# 18. peaklength
# 19. tss_coord
# 20. peak_genchr
# 21. peak_genstart
# 22. peak_genend
# 23. peak_start # peak start in native promoter seq
# 24. peak_end # peak end in native promoter seq
# 25. peak_summit # peak summit in native promoter seq
