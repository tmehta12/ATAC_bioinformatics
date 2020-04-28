#!/bin/sh


################################################################################################################

# ATAC-seq pipeline - Part 1
# March 2020: Tarang K. Mehta, Earlham Institute, Norwich, UK

################################################################################################################

# Script usage: ./ATAC_Bioinf_pipeline_v2a.sh -s "spID" -g "spG" -c "String C"
# After this, run

################################################################################################################

# ~ This pipeline can (should) be proceeded 'ATAC_Bioinf_pipeline_v2b.sh', that does:

# 1-8: Trimming, alignment, filtering, calling peaks, and annotation

################################################################################################################

# This may or may not be required - Setting parameters for command line input (AMEND BELOW IF REQUIRED!)

helpFunction()
{
   echo ""
   echo "Usage: $0 -s spID -g spG -c parameterC"
   echo -e "\t-s spID = Species ID, preferably two short letters, tissue and experiment e.g. Metriaclima zebra Liver = Mz_L_ATAC/gDNA"
   echo -e "\t-g spG = Species genome ID e.g. hg19 or M_zebra_UMD1"
   echo -e "\t-c Description of what is parameterC"
   exit 1 # Exit script after printing help
}

while getopts "s:g:c:" opt
do
   case "$opt" in
      a ) spID="$OPTARG" ;;
      b ) spG="$OPTARG" ;;
      c ) parameterC="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$spID" ] || [ -z "$spG" ] || [ -z "$parameterC" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "$spID"
echo "$spG"
echo "$parameterC"

################################################################################################################

# Add variables here:
scripts=(/tgac/workarea/group-vh/Tarang/ATACseq/2.run2) # place all scripts in the topmost directory - create this separately

rawreadfolder1=(/tgac/data/reads/TarangMehta_EI_EI_TM_ENQ-1771_A_03/170214_D00507_0252_AHG5KGBCXY)
rawreadfolder2=(/tgac/data/reads/TarangMehta_EI_EI_TM_ENQ-1771_A_03/170320_D00507_0270_AHHYLLBCXY)
rawreadfolder3=(/tgac/data/reads/)

WD=(/tgac/workarea/group-vh/Tarang/ATACseq/2.run2) # insert the working directory

rawreaddir=($WD/0.rawreads) # assign raw reads dir

libids=($scripts/libids.txt) # 2-col space-delimited file where col1 is the *_{R1,R2}.fastq.merged.gz and col2 is the species_tissue_experiment_barcode_{R1,R2}.fastq.merged.gz e.g. Mz_L_ATAC/gDNA

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



echo '# -- 0. File merging complete -- #'

################################################################################################################

### 1. Create the appropriate directory structure and file paths for downstream analysis (and to run following scripts)

# 1a. Create approporiate directory structure - this uses the space delimited file of 'merged fasta filename' and 'species_tissue_experiment_FASTA_filename' e.g. Mz_L_ATAC/gDNA for creating the symbolic links as per species structure

# Example is (note, these are made up examples)
# echo 'PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R1.fastq.merged.gz Pn1_T_ATAC_TCGCCTGC-AACCGCCA_L001_R1.fastq.merged.gz' > libids.txt
# echo 'PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R2.fastq.merged.gz Pn1_T_ATAC_TCGCCTGC-AACCGCCA_L001_R2.fastq.merged.gz' >> libids.txt

# A. read in the space delimited file and prepare working directories for all col2 entries

prefix=($scripts/prefix.txt)
reads=(reads.txt)


awk -F' ' '{print $2}' $libids | awk -F'_' '{print $1"_"$2"_"$3}' > $prefix # create a prefix file to iterate
mapfile -t prefixmap < $prefix # assign prefixes to $prefixmap

awk -F' ' '{print $2}' $libids > $reads
mapfile -t reads < $reads # ${reads[0]} calls read1 AND ${reads[1]} calls read2'

# B. within each species_tissue_experiment working directory, create a raw reads directory too

# 1b. Create symbolic links to raw reads in per species_tissue e.g. Mz_L and experiment (ATAC/gDNA) - no file renaming required at this stage; will be done after trimming


# 1a.

spWD=(/tgac/workarea/group-vh/Tarang/ATACseq/2.run2/$i) # insert the working directory
sprawreaddir=($spWD/0.rawreads) # assign raw reads dir

mkdir -p $sprawreaddir
