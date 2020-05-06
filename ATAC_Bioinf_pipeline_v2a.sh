#!/bin/sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 12000 # memory pool for all cores
#SBATCH -t 0-23:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

################################################################################################################

# ATAC-seq pipeline - Part 1
# March 2020: Tarang K. Mehta, Earlham Institute, Norwich, UK

################################################################################################################

### Script usage

# 1. Create a topmost directory separately, placing this script in there and running from that directory e.g. /tgac/workarea/group-vh/Tarang/ATACseq/2.run2
# 2. sbatch ATAC_Bioinf_pipeline_v2a.sh

################################################################################################################

# ~ This pipeline should be proceeded with 'ATAC_Bioinf_pipeline_v2b_gDNA.sh' and then 'ATAC_Bioinf_pipeline_v2b.sh' that does:

# 1-2: Trimming and alignment of gDNA reads
# 1-8: Trimming and alignment of ATAC reads, and then filtering, calling peaks, and annotation

################################################################################################################

# Add variables here:
scripts=(/tgac/workarea/group-vh/Tarang/ATACseq/2.run2) # place all scripts in the topmost directory - create this separately and place this script in there too

libids=($scripts/libids.txt) # 2-col space-delimited file where col1 is the *_{R1,R2}.fastq.merged.gz and col2 is the species_tissue_experiment_barcode_{R1,R2}.fastq.merged.gz e.g. Mz_L_ATAC/gDNA

rawreadfolder1=(/tgac/data/reads/TarangMehta_EI_EI_TM_ENQ-1771_A_03/170214_D00507_0252_AHG5KGBCXY)
rawreadfolder2=(/tgac/data/reads/TarangMehta_EI_EI_TM_ENQ-1771_A_03/170320_D00507_0270_AHHYLLBCXY)
rawreadfolder3=(/tgac/data/reads/)

WD=(/tgac/workarea/group-vh/Tarang/ATACseq/2.run2) # insert the working directory
rawreaddir=($WD/0.rawreads)

################################################################################################################

### 0. Merge files if sequenced over multiple lanes - Add the species ID and tissue to create read1 and read2 e.g. Mz_L_ATAC_read1 and Mz_L_ATAC_read2

# 0a. Create folder linked to reads
mkdir -p $rawreaddir
cd $rawreaddir

# 0b. Merge files of same libraries sequenced on different lanes

echo '# -- 0. File merging started -- #'

# Merge from test lane first
for lane1 in $rawreadfolder1/*.fastq.gz ; do lane2=$(echo $lane1| sed 's|/tgac/data/reads/TarangMehta_EI_EI_TM_ENQ-1771_A_03/170214_D00507_0252_AHG5KGBCXY/|/tgac/data/reads/TarangMehta_EI_EI_TM_ENQ-1771_A_03/170320_D00507_0270_AHHYLLBCXY/|g'); cat $lane1 $lane2 > $rawreaddir/"$(basename "$lane1" .gz).merged.gz" ; done

# Then, merge with new sequencing files {TO ADD}




################################################################################################################

### 1. Create the appropriate directory structure and file paths for downstream analysis (and to run following scripts)

# 1a. Create approporiate directory structure - this uses the space delimited file of 'merged fasta filename' and 'species_tissue_experiment_FASTA_filename' e.g. Mz_L_ATAC/gDNA for creating the symbolic links as per species structure

# Example is (note, these are made up examples)
# echo 'PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R1.fastq.merged.gz Pn1_T_ATAC_TCGCCTGC-AACCGCCA_L001_R1.fastq.merged.gz' > libids.txt
# echo 'PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R2.fastq.merged.gz Pn1_T_ATAC_TCGCCTGC-AACCGCCA_L001_R2.fastq.merged.gz' >> libids.txt

# A. read in the space delimited file and prepare working directories for all col2 entries
# B. within each species_tissue_experiment working directory, create a raw reads directory too

cd $WD

prefix=($scripts/prefix.txt)

echo '# -- 0. File merging complete -- #'

awk -F' ' '{print $2}' $libids | awk -F'_' '{print $1"_"$2"_"$3}' | sort -u > $prefix # create a prefix file to iterate
while IFS= read -r i; do
  # echo $i
  mkdir $i
  cd $i
  mkdir 0.rawreads
  cd ../
done < $prefix

# 1b. Create symbolic links to raw reads in per species_tissue e.g. Mz_L and experiment (ATAC/gDNA) - no file renaming required at this stage; will be done after trimming

cd $WD

libids2=($scripts/libids2.txt) # 2-col space-delimited file where col1 is the *_{R1,R2}.fastq.merged.gz and col2 is the species_tissue_experiment only e.g. Mz_L_ATAC/gDNA
sed 's/ATAC_.*/ATAC/g' $libids | sed 's/gDNA_.*/gDNA/g' > $libids2 # create the above

# while read, if statement to create symbolic link to raw reads for each species_tissue_experiment

while IFS= read -r i; do
  # echo $i
  spdir=$(echo $i | awk -F' ' '{print $2}') # pull out the species_tissue_experiment ID of each line
  # echo $spdir
  cd $WD/$spdir # change to each subdirectory
  dir=$(pwd) # create variable of the full path
  # echo $dir
  dir2=$(echo $dir | sed "s|$WD/||g") # pull out the species_tissue_experiment ID of the working path
  # echo $dir2
  read=$(echo $i | awk -F' ' '{print $1}') # pull out the corresponding read filename of each line
  # echo $read
  if [[ $spdir == $dir2 ]]; then # if the species_tissue_experiment ID in the list and the working directory match, then
    ln -s $rawreaddir/$read 0.rawreads/ # create symbolic link in raw reads dir to corresponding raw read
  fi
done < $libids2

################################################################################################################
