#!/bin/sh

##################################################################################################

#### A. calliptera 3 dpf, 7dpf and 12 dpf ATAC and RNA-seq analysis
## 2025

##################################################################################################

## Steps carried out prior to this:

# 1. Raw ATAC-seq reads (50 bp paired end sequenced on NovaSeq 6000) generating an average of 46 million reads per sample and 8 million reads per sample naked DNA control. 
# 2 Reads were trimmed using Trim Galore! (v0.5.0) and FastQC (v0.11.9) using default parameters.
# 3. Read alignment, post alignment filtering and ATAC peak calling were performed according to the ENCODE projects ‘ATAC-seq Data Standards and Processing Pipeline’ for replicated data (https://www.encodeproject.org/atac-seq/). Briefly:
#      o	Trimmed reads were mapped to the Astatotilapia calliptera ‘fAstCal1.2’ genome using bowtie2 (v2.2.6) with parameters ‘–k 4 –X2000 –mm’ and outputted in BAM format using SAMtools (v1.9).
#      o	Since ATAC-seq can generate a high proportion (15-50% in a typical experiment) of mitochondrial mapped reads, any reads mapping to their respective mitochondrial genomes were identified using BLAST (v2.3.0) 13 and removed from the BAM file using SAMtools (v1.9).
#      o	The resulting BAM files were sorted, and duplicated reads were marked using Sambamba v0.6.5.
#      o	Duplicated, unmapped, non-primary alignment, and failing platform QC reads were filtered using SAMtools (v1.9), retaining reads mapped as proper pairs, and fragment length distributions were plotted using Picard (v1.140).
#      o	BAM files were converted to tagalign files using Bedtools (v2.30.0) and Tn5 shifting of ATAC mappings carried out prior to peak calling. 
#      o	Peaks were identified using macs2 (v2.1.1) with the shifted tag as test and corresponding control DNA as input with parameters ‘-f BED -p 0.05 --nomodel --shift -75 --extsize 150’.
#      o	Narrow peaks were used to create coverage tracks using bedClip and bedToBigBed in the UCSC-tools package (v333).
#      o	Following the ENCODE pipeline, Irreproducible Discovery Rate (IDR) peaks of true replicates were flagged as either true (<0.1) or false (≥0.1) using idr v2.0.4.
# 4. Narrow peaks either overlapping or most proximal to all annotated features (5 kb gene promoters annotated too) in the Astatotilapia calliptera ‘fAstCal1.2’ genome was mapped using the intersect function of Bedtools (v2.30.0).
# 5. TF footprints were characterised using HINT-ATAC in the Regulatory Genomic Toolbox (v0.13.0) 14 using a stringent false positive rate (FPR) of 0.0001, with both species-specific and cichlid-wide position weight matrices (PWMs) as defined in our previous study 15, as well as vertebrate PWMs from JASPAR (v9.0), HOCOMOCO, GTRD, and UniPROBE.

##################################################################################################

## Copying the required files to Midhukrishna home directory (hc-storage)

cd /home/tmehta12/hc-storage

nano cpfiles.sh

#!/bin/bash
#SBATCH -p low # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --job-name=cpfiles # job name
#SBATCH --output=cp_%A_%a.out # standard out file name
#SBATCH --error=cp_%A_%a.err # standard error file name
#SBATCH --mem=8G # memory for the job

src=(/home/tmehta12/hc-storage/ATAC_Acalliptera/)
tgt=(/home/hsmidhuk/hc-storage/ATAC_Acalliptera)
mkdir -p ${tgt}
rsync -av --checksum ${src} ${tgt}

# run the above script
sbatch cpfiles.sh # DONE

## Prepare dir structure and files for the analysis

# Unzip files in directories # DONE
module load unzip/6.0

cd /home/hsmidhuk/hc-storage/ATAC_Acalliptera/2.Annotation
for f in *.zip; do unzip "$f" -d "${f%.zip}"; done # unzip to a directory with the same name as the zip file (excluding .zip)

cd /home/hsmidhuk/hc-storage/ATAC_Acalliptera/3.TFfprint_SignalTrack
for f in *.zip; do unzip "$f" -d "${f%.zip}"; done # unzip to a directory with the same name as the zip file (excluding .zip)

cd /home/hsmidhuk/hc-storage/ATAC_Acalliptera/4.Peak_Orth
unzip 4.Peak_Orth.zip

cd /home/hsmidhuk/hc-storage/ATAC_Acalliptera/4b.peak_phylop
unzip 4b.peak_phylop.zip


##################################################################################################

## Steps to follow:

# 1.	Literature Review and Data Collection
# •	Familiarise with the scientific area and types of data by reading and making notes on relevant literature.
# •	Familiarise and understand the methods already applied above and prepare a flowchart to represent this (could use https://app.diagrams.net/) – this will form a good intro to the methods you will apply. 
# •	Prepare a flowchart to diagrammatically represent the methods you will apply below – you will add to this as you progress (could use https://app.diagrams.net/)

# 2.	Peak annotation analysis - use ggplot in R studio for plotting
# 2.1. (Supp Fig. 1) Pie chart showing the distribution of narrow open-chromatin peaks in annotated features across the three stages
# 2.2. (Fig. 1b) Barplot (including error bars) of mean number of peaks based on distance from TSS (so 2 bars for each stage - +/- 100 bp from TSS or >=100bp from TSS)
# 2.3. (Fig. 1c) Carry out peak enrichment in annotations using Genomic Association Tester, GAT: https://gat.readthedocs.io/en/latest/) with multiple test correction of FDR<0.05, and prepare a dotplot with facets (for each stage) - Circles show enriched annotated region (y-axis) of significance (FDR <0.05, to right) and fold enrichment (x-axis) values of all annotated peaks. Number of peaks overlapping each annotation shown by size of each circle. 

# Files to use:

# a. final annotation files (if you need)
dir=/home/hsmidhuk/hc-storage/ATAC_Acalliptera/4b.peak_phylop/final_annot
acannotbed=${dir}/Astatotilapia_calliptera.fAstCal1.2.annot.bed
acannotgff=${dir}/Astatotilapia_calliptera.fAstCal1.2.annot.gff

# b. distribution of peaks in annotated features - for 2.1 and 2.2 above

mkdir -p /home/hsmidhuk/hc-storage/ATAC_Acalliptera/Analysis/1.peakannot
cd /home/hsmidhuk/hc-storage/ATAC_Acalliptera/Analysis/1.peakannot

peakannotdir=/home/hsmidhuk/hc-storage/ATAC_Acalliptera/4b.peak_phylop/0.phylop_run/Astatotilapia_calliptera
peakannotdir=/home/hsmidhuk/hc-storage/ATAC_Acalliptera/4b.peak_phylop/0.phylop_run

# Cols in the *.narrowPeak.features.gff are:

# 1. chr/scaff
# 2. peak ID - as generated and assigned by MACS
# 3. annotation - intergenic, CNEs, 5kb gene promoter, 5' UTR, exon, intron, 3’ UTR
# 4. start - The starting position of the feature in the sequence. The first base is numbered 1.
# 5. end - The ending position of the feature (inclusive).
# 6. score - Contains the scaled IDR value, min(int(log2(-125IDR), 1000). e.g. peaks with an IDR of 0 have a score of 1000, idr 0.05 have a score of int(-125log2(0.05)) = 540, and idr 1.0 has a score of 0.
# 7. strand - Valid entries include "+", "-", or "." (for don't know/don't care).
# 8. frame - all entries will be '.'
# 9. grouped (separated by ; -  7 values in total)
  # MACS2 peak IDR output- signalValue float - fold-change at peak summit: Measurement of enrichment for the region for merged peaks. When a peak list is provided this is the value from the peak list.
  # MACS2 peak IDR output- p-value float: Merged peak p-value. When a peak list is provided this is the value from the peak list.
  # MACS2 peak IDR output- q-value float: Merged peak q-value. When a peak list is provided this is the value from the peak list.
  # MACS2 peak IDR output- summit int: Merged peak summit
  # MACS2 peak IDR output- IDR True (T) or False (F) of peaks passing IDR threshold of 10%
  # Ensembl Gene ID
  # Ensembl Gene Symbol

# b.1 - collate the peak count data
for i in ${peakannotdir}/Astatotilapia_calliptera/Astatotilapia_calliptera_ATAC_peaks.final.narrowPeak.features.gff; do
  sp=$(echo ${i} | sed 's|.*/||g' | sed 's|_ATAC_peaks.final.narrowPeak.features.gff||g')
  # echo $sp
  awk -v sp=$sp '{print sp,$3,$9}' OFS='\t' ${i} | awk -F';' '{print $1,$5}' OFS='\t' | awk '{print $1,$2":"$4}' OFS='\t' | awk '{a[$0]++}END{for(i in a){print i, a[i]}}' >> 1b-i.peakannot.counts.tmp1
done

sort -V -k1,1 -k2,2 1b-i.peakannot.counts.tmp1 | sed 's| |\t|g' > 1b-i.peakannot.counts
rm 1b-i.peakannot.counts.tmp1
# final file to use for IDR True and False stats: 1b-i.peakannot.counts

for i in ${peakannotdir}/Astatotilapia_calliptera/Astatotilapia_calliptera_ATAC_peaks.final.narrowPeak.features.gff; do
  sp=$(echo ${i} | sed 's|.*/||g' | sed 's|_ATAC_peaks.final.narrowPeak.features.gff||g')
  # echo $sp
  awk -v sp=$sp '{print sp,$3}' OFS='\t' ${i} | awk '{a[$0]++}END{for(i in a){print i, a[i]}}' >> 1b-i.peakannot.counts.tmp1
done
sort -V -k1,1 -k2,2 1b-i.peakannot.counts.tmp1 > 1b-i.peakannot.counts2 # NOTE, THIS ADDS FALSE AND TRUE IDR PEAKS FOR EACH ANNOT!!!
rm 1b-i.peakannot.counts.tmp1

# >> create a supplementary table with the above data > Supp Table S2)
# >> use '/home/hsmidhuk/hc-storage/ATAC_Acalliptera/Analysis/1.peakannot/1b-i.peakannot.counts2' to create a pie chart in R

# b.2 - two plots to prepare
# b.2.a - Fig 1b - barplot, with error bars for peaks that fall +-100bp OR >=100bp from TSS (main plot to include in figure) in each stage - Summit_to_TSS is col22 and IDR is col15 (note some are NULL so remove those as it means that OG has no peak!) here: ../4.Peak_Orth/*_ATAC_peaks.final.narrowPeak.promoverlap.orth5
# b.2.b - pie chart for the distribution of peaks in annotated features across the three stages

# create a single nine column file - species, tissue, total_peaks, total_+-100TSS, IDR_T_+-100TSS, IDR_F_+-100TSS, total_>100TSS, IDR_T_>100TSS, IDR_F_>100TSS
tssdir=/home/hsmidhuk/hc-storage/ATAC_Acalliptera/4.Peak_Orth

printf 'spID\tspecies\ttissue\ttotalpeaks\ttssnear\ttssnearidrt\ttssnearidrf\ttssdist\ttssdistidrt\ttssdistidrf\n' > 1b-ii.tss_stats.out1

for i in ${tssdir}/*_ATAC_peaks.final.narrowPeak.promoverlap.orth5; do
  sp=$(echo $(basename "${i}") | sed 's|_ATAC_peaks.final.narrowPeak.promoverlap.orth5||g' | awk -F'_' '{print $1}')
  spID=$(echo $(basename "${i}") | sed 's|_ATAC_peaks.final.narrowPeak.promoverlap.orth5||g' | awk -F'_' '{print $1}' | sed 's/1a//g' | sed 's/1b//g' | sed 's/2a//g' | sed 's/2b//g' | sed 's/3a//g' | sed 's/3b//g')
  tissue=$(echo $(basename "${i}") | sed 's|_ATAC_peaks.final.narrowPeak.promoverlap.orth5||g' | awk -F'_' '{print $2}')
  totalpeaks=$(awk '$15!="NULL" && $15!="IDR"' ${i} | wc -l | awk '{print $1}')
  tssnear=$(awk '$15!="NULL" && $15!="IDR"' ${i} | awk '$22<=100 && $22>=-100' | wc -l | awk '{print $1}')
  tssnearidrt=$(awk '$15!="NULL" && $15!="IDR"' ${i} | awk '$22<=100 && $22>=-100 && $15=="T"' | wc -l | awk '{print $1}')
  tssnearidrf=$(awk '$15!="NULL" && $15!="IDR"' ${i} | awk '$22<=100 && $22>=-100 && $15=="F"' | wc -l | awk '{print $1}')
  tssdist=$(awk '$15!="NULL" && $15!="IDR"' ${i} | awk '$22>100 || $22<-100' | wc -l | awk '{print $1}')
  tssdistidrt=$(awk '$15!="NULL" && $15!="IDR"' ${i} | awk '$22>100 || $22<-100' | awk '$15=="T"' | wc -l | awk '{print $1}')
  tssdistidrf=$(awk '$15!="NULL" && $15!="IDR"' ${i} | awk '$22>100 || $22<-100' | awk '$15=="F"' | wc -l | awk '{print $1}')
  echo -e ${spID}'_'${tissue}'\t'${sp}'\t'${tissue}'\t'${totalpeaks}'\t'${tssnear}'\t'${tssnearidrt}'\t'${tssnearidrf}'\t'${tssdist}'\t'${tssdistidrt}'\t'${tssdistidrf} >> 1b-ii.tss_stats.out2
done

# then calculate the average (based on number of samples per stage) of each column to use for plotting
for j in Ac; do
  grep $j 1b-ii.tss_stats.out2 | awk '$3=="3dpf"' > 1b-ii.tss_stats.out2.$j.3dpf
  grep $j 1b-ii.tss_stats.out2 | awk '$3=="7dpf"' > 1b-ii.tss_stats.out2.$j.7dpf
  grep $j 1b-ii.tss_stats.out2 | awk '$3=="12dpf"' > 1b-ii.tss_stats.out2.$j.12dpf
done

printf 'species\ttissue\ttotalpeaksavg\ttssnearavg\ttssnearidrtavg\ttssnearidrfavg\ttssdistavg\ttssdistidrtavg\ttssdistidrfavg\n' > 1b-ii.tss_stats.out3

for i in 1b-ii.tss_stats.out2.*.*; do
  sp=$(echo $(basename "${i}") | sed 's|1b-ii.tss_stats.out2.||g' | awk -F'.' '{print $1}')
  tissue=$(echo $(basename "${i}") | sed 's|1b-ii.tss_stats.out2.||g' | awk -F'.' '{print $2}')
  totalpeaksavg=$(awk '{ sum += $4 } END { if (NR > 0) print sum / NR }' ${i})
  tssnearavg=$(awk '{ sum += $5 } END { if (NR > 0) print sum / NR }' ${i})
  tssnearidrtavg=$(awk '{ sum += $6 } END { if (NR > 0) print sum / NR }' ${i})
  tssnearidrfavg=$(awk '{ sum += $7 } END { if (NR > 0) print sum / NR }' ${i})
  tssdistavg=$(awk '{ sum += $8 } END { if (NR > 0) print sum / NR }' ${i})
  tssdistidrtavg=$(awk '{ sum += $9 } END { if (NR > 0) print sum / NR }' ${i})
  tssdistidrfavg=$(awk '{ sum += $10 } END { if (NR > 0) print sum / NR }' ${i})
  echo -e ${sp}'\t'${tissue}'\t'${totalpeaksavg}'\t'${tssnearavg}'\t'${tssnearidrtavg}'\t'${tssnearidrfavg}'\t'${tssdistavg}'\t'${tssdistidrtavg}'\t'${tssdistidrfavg} >> 1b-ii.tss_stats.out3 # file not needed but creating for reference of IDR
done

# final file to use to create BOTH plots - standard error/standard deviation barplot (Fig. 1b) AND pie chart: '1b-ii.tss_stats.out2'



# c. calculate the average of peaks that are in up to 5 kb gene promoter regions - this is to report in main text write-up

tssdir=/home/hsmidhuk/hc-storage/ATAC_Acalliptera/4.Peak_Orth

printf 'spID\tspecies\ttissue\ttotalpeaks\ttss5kb\ttss5kbidrt\ttss5kbidrf\n' > 1b-ii.5kb_stats.out1

for i in ${tssdir}/*_ATAC_peaks.final.narrowPeak.promoverlap.orth5; do
  sp=$(echo $(basename "${i}") | sed 's|_ATAC_peaks.final.narrowPeak.promoverlap.orth5||g' | awk -F'_' '{print $1}')
  spID=$(echo $(basename "${i}") | sed 's|_ATAC_peaks.final.narrowPeak.promoverlap.orth5||g' | awk -F'_' '{print $1}' | sed 's/1a//g' | sed 's/1b//g' | sed 's/2a//g' | sed 's/2b//g' | sed 's/3a//g' | sed 's/3b//g')
  tissue=$(echo $(basename "${i}") | sed 's|_ATAC_peaks.final.narrowPeak.promoverlap.orth5||g' | awk -F'_' '{print $2}')
  totalpeaks=$(awk '$15!="NULL" && $15!="IDR"' ${i} | wc -l | awk '{print $1}')
  tssnear=$(awk '$15!="NULL" && $15!="IDR"' ${i} | awk '$22<=5000 && $22>=-5000' | wc -l | awk '{print $1}')
  tssnearidrt=$(awk '$15!="NULL" && $15!="IDR"' ${i} | awk '$22<=5000 && $22>=-5000 && $15=="T"' | wc -l | awk '{print $1}')
  tssnearidrf=$(awk '$15!="NULL" && $15!="IDR"' ${i} | awk '$22<=5000 && $22>=-5000 && $15=="F"' | wc -l | awk '{print $1}')
  echo -e ${spID}'_'${tissue}'\t'${sp}'\t'${tissue}'\t'${totalpeaks}'\t'${tssnear}'\t'${tssnearidrt}'\t'${tssnearidrf} >> 1b-ii.5kb_stats.out2
done


# then calculate the average (based on number of samples per tissue) of each column, per species to use for plotting
for j in Ac; do
  grep $j 1b-ii.5kb_stats.out2 | awk '$3=="3dpf"' > 1b-ii.5kb_stats.out2.$j.3dpf
  grep $j 1b-ii.5kb_stats.out2 | awk '$3=="7dpf"' > 1b-ii.5kb_stats.out2.$j.7dpf
  grep $j 1b-ii.5kb_stats.out2 | awk '$3=="12dpf"' > 1b-ii.5kb_stats.out2.$j.12dpf
done

printf 'species\ttissue\ttotalpeaksavg\ttss5kbavg\ttss5kbidrtavg\ttss5kbidrfavg\n' > 1b-ii.5kb_stats.out3

for i in 1b-ii.5kb_stats.out2.*.*; do
  sp=$(echo $(basename "${i}") | sed 's|1b-ii.5kb_stats.out2.||g' | awk -F'.' '{print $1}')
  tissue=$(echo $(basename "${i}") | sed 's|1b-ii.5kb_stats.out2.||g' | awk -F'.' '{print $2}')
  totalpeaksavg=$(awk '{ sum += $4 } END { if (NR > 0) print sum / NR }' ${i})
  tssnearavg=$(awk '{ sum += $5 } END { if (NR > 0) print sum / NR }' ${i})
  tssnearidrtavg=$(awk '{ sum += $6 } END { if (NR > 0) print sum / NR }' ${i})
  tssnearidrfavg=$(awk '{ sum += $7 } END { if (NR > 0) print sum / NR }' ${i})
  echo -e ${sp}'\t'${tissue}'\t'${totalpeaksavg}'\t'${tssnearavg}'\t'${tssnearidrtavg}'\t'${tssnearidrfavg} >> 1b-ii.5kb_stats.out3
done


##**** NOTE: 'collatepeakfeature.sh' BELOW HAS ALREADY BEEN RAN SO NO NEED TO RUN AGAIN
##**** Code included here so aware of how it was generated!

# d. test enrichment of all tissue-specific in all peaks IDR T and F > (for 2.3 above that will be Fig. 1c)

# Create a FINAL set of A) collated (by stage) narrow peaks (IDR T and F), and B) collated (by tissue) IDR T and F narrow peaks according to feature
peakannotdir=/home/hsmidhuk/hc-storage/ATAC_Acalliptera/4b.peak_phylop/0.phylop_run

cd $peakannotdir

nano collatepeakfeature.sh

#!/bin/bash -e
#SBATCH -p low # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 28000 # memory pool for all cores
#SBATCH -t 0-05:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR

peakannotdir=/home/hsmidhuk/hc-storage/ATAC_Acalliptera/4b.peak_phylop/0.phylop_run

cd $peakannotdir/Astatotilapia_calliptera
cat 1aAc_3dpf_ATAC_peaks.final.narrowPeak 1bAc_3dpf_ATAC_peaks.final.narrowPeak | sort -k1,1 -k2,2n > Ac_3dpf_ATAC_peaks.final.allIDR.narrowPeak
cat 1aAc_3dpf_ATAC_peaks.final.narrowPeak.features.gff 1bAc_3dpf_ATAC_peaks.final.narrowPeak.features.gff | sort -k1,1 -k2,2n > Ac_3dpf_ATAC_peaks.final.allIDR.narrowPeak.features.gff
cat 2aAc_7dpf_ATAC_peaks.final.narrowPeak 2bAc_7dpf_ATAC_peaks.final.narrowPeak | sort -k1,1 -k2,2n > Ac_7dpf_ATAC_peaks.final.allIDR.narrowPeak
cat 2aAc_7dpf_ATAC_peaks.final.narrowPeak.features.gff 2bAc_7dpf_ATAC_peaks.final.narrowPeak.features.gff | sort -k1,1 -k2,2n > Ac_7dpf_ATAC_peaks.final.allIDR.narrowPeak.features.gff
cat 3aAc_12dpf_ATAC_peaks.final.narrowPeak 3bAc_12dpf_ATAC_peaks.final.narrowPeak | sort -k1,1 -k2,2n > Ac_12dpf_ATAC_peaks.final.allIDR.narrowPeak
cat 3aAc_12dpf_ATAC_peaks.final.narrowPeak.features.gff 3bAc_12dpf_ATAC_peaks.final.narrowPeak.features.gff | sort -k1,1 -k2,2n > Ac_12dpf_ATAC_peaks.final.allIDR.narrowPeak.features.gff


### FINAL IDR true AND false peak and peaks overlapping features files are (use THESE all peaks for the paper and then select IDR true for examples!):
## /home/hsmidhuk/hc-storage/ATAC_Acalliptera/4b.peak_phylop/0.phylop_run/*/*_*_ATAC_peaks.final.allIDR.narrowPeak
## /home/hsmidhuk/hc-storage/ATAC_Acalliptera/4b.peak_phylop/0.phylop_run/*/*_*_ATAC_peaks.final.allIDR.narrowPeak.features.gff

# run the above
sbatch collatepeakfeature.sh # DONE - ALREADY RAN SO DO NOT RE-RUN!


##**** NOTE: NEED TO INSTALL ON HPC USING SINGULARITY - Genomic Association Tester (gat): https://gat.readthedocs.io/en/latest/installation.html
# Usage instructions here - https://gat.readthedocs.io/en/latest/usage.html

nano 1.peakenrich_Ab.allIDR.sh

#!/bin/bash
#SBATCH -p ei-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 1-23:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

# Run the Genomic Association Tester (gat): https://gat.readthedocs.io/en/latest/usage.html

peakannotdir=/home/hsmidhuk/hc-storage/ATAC_Acalliptera/4b.peak_phylop/0.phylop_run

# This takes, as input:
# A set of intervals S with segments of interest to test – peaks
for i in 3dpf 7dpf 12dpf; do
  cut -f1-6,11 $peakannotdir/Astatotilapia_calliptera/Ac_${i}_ATAC_peaks.final.allIDR.narrowPeak | sort -V -k1,1 -k2,2n > Ac_${i}_ATAC_peaks.collated.allIDR.narrowPeak
done

# A set of intervals A with annotations to test against – bed file of Annotations:
cut -f1-6 /home/hsmidhuk/hc-storage/ATAC_Acalliptera/4b.peak_phylop/final_annot/Astatotilapia_calliptera.fAstCal1.2.annot.bed > Astatotilapia_calliptera.fAstCal1.2.annot.basic.bed

# A set of intervals W describing a workspace - Annotation of whole length of chromosomes (just lengths)
awk '{print $1,"0",$2}' OFS='\t' /home/hsmidhuk/hc-storage/ATAC_Acalliptera/4b.peak_phylop/Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.and.nonchromosomal.fa.chrom.sizes > Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.and.nonchromosomal.fa.chrom.sizes

gat-run.py --verbose=5 \
           --log=gatnormed_Ac_3dpf_Peak_AllAnnot_allIDR_ChrBgd.tsv.log \
            --segments=Ac_3dpf_ATAC_peaks.collated.allIDR.narrowPeak \
            --annotations=Astatotilapia_calliptera.fAstCal1.2.annot.basic.bed \
            --workspace=Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.and.nonchromosomal.fa.chrom.sizes \
            --counter=segment-overlap \
            --ignore-segment-tracks \
            --qvalue-method=BH \
            --pvalue-method=norm \
            >& gatnormed_Ac_3dpf_Peak_AllAnnot_allIDR_ChrBgd.tsv

gat-run.py --verbose=5 \
           --log=gatnormed_Ac_7dpf_Peak_AllAnnot_allIDR_ChrBgd.tsv.log \
            --segments=Ac_7dpf_ATAC_peaks.collated.allIDR.narrowPeak \
            --annotations=Astatotilapia_calliptera.fAstCal1.2.annot.basic.bed \
            --workspace=Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.and.nonchromosomal.fa.chrom.sizes \
            --counter=segment-overlap \
            --ignore-segment-tracks \
            --qvalue-method=BH \
            --pvalue-method=norm \
            >& gatnormed_Ac_7dpf_Peak_AllAnnot_allIDR_ChrBgd.tsv

gat-run.py --verbose=5 \
           --log=gatnormed_Ac_12dpf_Peak_AllAnnot_allIDR_ChrBgd.tsv.log \
            --segments=Ac_12dpf_ATAC_peaks.collated.allIDR.narrowPeak \
            --annotations=Astatotilapia_calliptera.fAstCal1.2.annot.basic.bed \
            --workspace=Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.and.nonchromosomal.fa.chrom.sizes \
            --counter=segment-overlap \
            --ignore-segment-tracks \
            --qvalue-method=BH \
            --pvalue-method=norm \
            >& gatnormed_Ac_12dpf_Peak_AllAnnot_allIDR_ChrBgd.tsv


# run the above
sbatch 1.peakenrich_Ac.allIDR.sh

# collate the enrichments, adding species and stage
for i in gatnormed_*_*_Peak_AllAnnot_allIDR_ChrBgd.tsv; do
  awk '{print FILENAME, $0}' OFS='\t' ${i} | sed '1,1d' | sed 's|gatnormed_||g' | sed 's|_Peak_AllAnnot_allIDR_ChrBgd.tsv||g' | sed 's|_|\t|' | cut -f1-2,4-27 > ${i}.tmp1
done

printf 'species\tstage\tannotation\tobserved\texpected\tCI95low\tCI95high\tstddev\tfold\tl2fold\tpvalue\tqvalue\ttrack_nsegments\ttrack_size\ttrack_density\tannotation_nsegments\tannotation_size\tannotation_density\toverlap_nsegments\toverlap_size\toverlap_density\tpercent_overlap_nsegments_track\tpercent_overlap_size_track\tpercent_overlap_nsegments_annotation\tpercent_overlap_size_annotation\n' > gatnormed_header

cat gatnormed_header gatnormed_*_*_Peak_AllAnnot_allIDR_ChrBgd.tsv.tmp1 > gatnormed_collated_Peak_AllAnnot_allIDR_ChrBgd.tsv

# then, plot the enrichments 
# Tarang has example plots in (~/Documents/TGAC/Technology_Development/ATAC-Seq/ATAC_manuscript/ATAC_Bioinf_ManuscriptPlots.R) but have a go yourself!


# 3.	Functional enrichment analysis
# •	Perform Gene Ontology enrichment analysis of narrow peaks overlapping gene promoter regions using the ‘g:GOst’ module of g:Profiler (https://biit.cs.ut.ee/gprofiler/gost), and plot output using ggplot.  

# Files to use: /home/hsmidhuk/hc-storage/ATAC_Acalliptera/4b.peak_phylop/0.phylop_run/Astatotilapia_calliptera/Ac_*dpf_ATAC_peaks.final.allIDR.narrowPeak.features.gff
# To do: extract the 'ENSACLG' ids as a list for each of the three stages (so three files) to upload into g:GOst


# 4.	Differential peak analysis
# •	Differential chromatin accessibility analysis using edgeR (see https://ivanek.github.io/analysisOfGenomicsDataWithR/12_ATACSeq_html.html#differential-accessibility-analysis).

# Files to use:


# 5.	Biological interpretation
# •	Define active TF footprints in gene promoter regions of key developmental genes (based on genes in ontology terms of point 2 above), and perform hierarchical clustering, based on TF bit-score of motif match in R, to identify core or developmental stage specific binding patterns.     
 
# Files to use:


# 6.	TF binding dynamics
# •	Extra: Use principles from the TOBIAS framework (https://github.com/loosolab/TOBIAS) 16 to further interrogate TF binding dynamics during embryonic development e.g., rewiring of binding networks. 
# •	Extra: Research literature on how you could model protein-DNA complexes and compare them between developmental stages in homologous loci. 
