#!/bin/sh

#!/bin/bash -e
#SBATCH -p ei-medium # partition (queue)
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
## NOTE: all pre-requisites to run all pipelines e.g. creating paths/files for merging/renaming etc. is ran in ATAC_Bioinf_pipeline_v2-abc_RunningData.sh

# 1. Create a topmost directory separately, placing this script in there and run from that directory e.g. /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/data/ATACseq/3.run2
# 2. Create a tab-delimited file for all files that require merging - each line is for each replicate, number of columns on each row does not matter. Place in topmost directory ($WD).
# 3. Create a 2-column space-delimited file where col1 is the *_{R1,R2}.fastq.merged.gz and col2 is the species_tissue_experiment_barcode_{R1,R2}.fastq.merged.gz e.g. Mz_L_ATAC/gDNA. Place in topmost directory ($WD).
# 4. sbatch ATAC_Bioinf_pipeline_v2a.sh

################################################################################################################

# ~ This pipeline should be proceeded with 'ATAC_Bioinf_pipeline_v2b_gDNA.sh' and then 'ATAC_Bioinf_pipeline_v2b.sh' that does:

# 1-2: Trimming and alignment of gDNA reads
# 1-8: Trimming and alignment of ATAC reads, and then filtering, calling peaks, and annotation

################################################################################################################

# Add variables here:
scripts=(/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/data/ATACseq/3.run2) # place all scripts in the topmost directory - create this separately and place this script in there too

libids=($scripts/libids.txt) # 2-col space-delimited file where col1 is the *_{R1,R2}.fastq.merged.gz and col2 is the species_tissue_experiment_barcode_{R1,R2}.fastq.merged.gz e.g. Mz_L_ATAC/gDNA

rawreadfolder1=(/tgac/data/reads/TarangMehta_EI_EI_TM_ENQ-1771_A_03/170214_D00507_0252_AHG5KGBCXY)
rawreadfolder2=(/tgac/data/reads/TarangMehta_EI_EI_TM_ENQ-1771_A_03/170320_D00507_0270_AHHYLLBCXY)
rawreadfolder3=(/ei/data/reads/PIP-2582/200717_A00478_0125_BHTF7MDMXX)

WD=(/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/data/ATACseq/3.run2) # insert the working directory
rawreaddir=($WD/0.rawreads)

prefix=($scripts/prefix.txt)
prefixgDNA=($scripts/prefixgDNA.txt)
prefixATAC=($scripts/prefixATAC.txt)

gDNAscr=(ATAC_Bioinf_pipeline_v2b_gDNA.sh)
ATACscr=(ATAC_Bioinf_pipeline_v2b.sh)

mergefiles1=($scripts/mergefiles_1.txt) # create a tab-delimited file where col1 is the sampleID e.g. Pnm1_L_gDNA, and all proceeding columns of each line are all the files that you want to merge e.g. replicates sequenced on different lanes. This will be automatically merged below so column numbers for each row can be different.
mergefiles2=($scripts/mergefiles_2.txt) # this file removes the sample ID - each line is all the files that you want to merge e.g. replicates sequenced on different lanes. This will be automatically merged below so column numbers for each row can be different.

################################################################################################################

### 0. Merge files if sequenced over multiple lanes - Add the species ID and tissue to create read1 and read2 e.g. Mz_L_ATAC_read1 and Mz_L_ATAC_read2

# 0a. Create folder linked to reads
mkdir -p $rawreaddir
cd $rawreaddir

# 0b. Merge files of same libraries sequenced on different lanes

echo '# -- 0. File merging started -- #'


# Using a a tab-delimited file where each line is all the files that you want to merge e.g. replicates sequenced on different lanes:
# a. Count the total number of columns (store as variable using awk)
# b. Generate alphanumeric variables for while loop according to number of columns.
# c. Use the column variables in a while IFS='\t' read -r A B C D XX loop, to merge the files and output using basename of col1

# this will get the filenames of run1 and run2
# for lane1 in $rawreadfolder1/*.fastq.gz; do
#   lane2=$(echo $lane1| sed 's|/tgac/data/reads/TarangMehta_EI_EI_TM_ENQ-1771_A_03/170214_D00507_0252_AHG5KGBCXY/|/tgac/data/reads/TarangMehta_EI_EI_TM_ENQ-1771_A_03/170320_D00507_0270_AHHYLLBCXY/|g')
#   echo -e $lane1' '$lane2
# done
#
# echo -e 'XX\tXX\tXX' > $mergefiles
# echo -e 'XX\tXX' >> $mergefiles

# This part (of creating paths/files for merging/renaming etc.) is ran in ATAC_Bioinf_pipeline_v2-abc_RunningData.sh so change there for specific files whilst using this pipeline
# 1. Creating the $mergefiles file in command line/excel and then port over (change line endings!)
# 2. For files that could be 'mixed' samples due to index hopping - create duplicates of the file but with different sample prefixes e.g. 63ATACNb4L_On2T and On2T_63ATACNb4L
# 3. Create the $libids file using $mergefiles1 - this is specific for your files so do change

cut -f 2- $mergefiles1 > $mergefiles2 # remove the first column
mergecols=$(awk -F'\t' '{print NF; exit}' $mergefiles2)
varstr=(variablestring.txt)
echo 'A B C D E F G H I J K L M N O' | tr ' ' '\n' > $varstr # create a sequence of letters to use as variables
mergecolsvar1=$(head -$mergecols $varstr | awk -F'\n' '{if(NR == 1) {printf $0} else {printf " "$0}}') # this will automatically create variable IDs e.g. 1,2..XX according to number of cols in merge file
mergecolsvar2=$(head -$mergecols $varstr | sed 's/^/\$/g' | awk -F'\n' '{if(NR == 1) {printf $0} else {printf " "$0}}') # this will automatically create variables e.g. $1,$2..$XX according to number of cols in merge file

while IFS=$'\t' read -r $(echo $mergecolsvar1); do
  # echo $A $B $C
  # eval echo ${mergecolsvar2}
  cat $(eval echo ${mergecolsvar2}) > "$(basename "$A" .gz).merged.gz"
done < $mergefiles2 ## Ensure final output is *_R1.fastq.merged.gz/*_R2.fastq.merged.gz

# # # Merge from test lane first
# for lane1 in $rawreadfolder1/*.fastq.gz ; do lane2=$(echo $lane1| sed 's|/tgac/data/reads/TarangMehta_EI_EI_TM_ENQ-1771_A_03/170214_D00507_0252_AHG5KGBCXY/|/tgac/data/reads/TarangMehta_EI_EI_TM_ENQ-1771_A_03/170320_D00507_0270_AHHYLLBCXY/|g'); cat $lane1 $lane2 > $rawreaddir/"$(basename "$lane1" .gz).merged.gz" ; done

################################################################################################################

### 1. Create the appropriate directory structure and file paths for downstream analysis (and to run following scripts)

# 1a. Create appropriate directory structure - this uses the space delimited file of 'merged fasta filename' and 'species_tissue_experiment_FASTA_filename' e.g. Mz_L_ATAC/gDNA for creating the symbolic links as per species structure

# Example is (note, these are made up examples)
# echo 'PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R1.fastq.merged.gz Pn1_T_ATAC_TCGCCTGC-AACCGCCA_L001_R1.fastq.merged.gz' > libids.txt
# echo 'PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R2.fastq.merged.gz Pn1_T_ATAC_TCGCCTGC-AACCGCCA_L001_R2.fastq.merged.gz' >> libids.txt

# A. read in the space delimited file and prepare working directories for all col2 entries
# B. within each species_tissue_experiment working directory, create a raw reads directory too

cd $WD

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
sed 's/ATAC_.*/ATAC/g' $libids | sed 's/_gDNA_.*/_gDNA/g' > $libids2 # create the above

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

### 2. Add relevant scripts to each directory

# Add each ATAC_Bioinf_pipeline_v2b_gDNA.sh script to each gDNA folder and ATAC_Bioinf_pipeline_v2b.sh to each ATAC folder

cd $WD

grep 'gDNA' $prefix > $prefixgDNA
grep 'ATAC' $prefix > $prefixATAC

while IFS= read -r i; do
  cp $gDNAscr $i
done < $prefixgDNA

while IFS= read -r i; do
  cp $ATACscr $i
done < $prefixATAC
