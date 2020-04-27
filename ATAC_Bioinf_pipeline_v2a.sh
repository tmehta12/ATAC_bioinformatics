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
   echo -e "\t-s spID = Species ID, preferably two short letters and tissue e.g. Metriaclima zebra Liver = Mz_L"
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

################################################################################################################

### 0. Merge files if sequenced over multiple lanes - Add the species ID and tissue to create read1 and read2 e.g. Mz_L_read1 and Mz_L_read2

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
