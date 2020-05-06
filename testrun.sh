#!/bin/sh

mkdir -p /tgac/workarea/group-vh/Tarang/ATACseq/2.run2

scripts=(/tgac/workarea/group-vh/Tarang/ATACseq/2.run2) # place all scripts in the topmost directory - create this separately and place this script in there too
libids=($scripts/libids.txt) # 2-col space-delimited file where col1 is the *_{R1,R2}.fastq.merged.gz and col2 is the species_tissue_experiment_barcode_{R1,R2}.fastq.merged.gz e.g. Mz_L_ATAC/gDNA
WD=(/tgac/workarea/group-vh/Tarang/ATACseq/2.run2) # insert the working directory
# rawreaddir=($WD/0.rawreads)
rawreaddir=(/tgac/workarea/group-vh/Tarang/ATACseq/1.run1jan2017_twolanes/0.rawreads)
libids=($scripts/libids.txt) # 2-col space-delimited file where col1 is the *_{R1,R2}.fastq.merged.gz and col2 is the species_tissue_experiment_barcode_{R1,R2}.fastq.merged.gz e.g. Mz_L_ATAC/gDNA


cd $scripts

nano libids.txt # # 2-col space-delimited file where col1 is the *_{R1,R2}.fastq.merged.gz and col2 is the species_tissue_experiment_barcode_{R1,R2}.fastq.merged.gz e.g. Mz_L_ATAC/gDNA

PRO1563_S1_lib_ACCACTGT-TCGTCGGC_L001_R1.fastq.merged.gz Pnm1_L_gDNA_PRO1563_S1_lib_ACCACTGT-TCGTCGGC_L001_R1.fastq.merged.gz
PRO1563_S1_lib_CAGAATGC-GAACTGAG_L001_R1.fastq.merged.gz Ab5_L_ATAC_PRO1563_S1_lib_CAGAATGC-GAACTGAG_L001_R1.fastq.merged.gz
PRO1563_S1_lib_CAGAATGC-TAACTCTA_L001_R1.fastq.merged.gz Ab5_L_gDNA_PRO1563_S1_lib_CAGAATGC-TAACTCTA_L001_R1.fastq.merged.gz
PRO1563_S1_lib_CAGAGAGG-TCGTCGGC_L001_R1.fastq.merged.gz Pnm1_T_gDNA_PRO1563_S1_lib_CAGAGAGG-TCGTCGGC_L001_R1.fastq.merged.gz
PRO1563_S1_lib_CAGAGGTA-GAAGGACG_L001_R1.fastq.merged.gz Ab5_T_ATAC_PRO1563_S1_lib_CAGAGGTA-GAAGGACG_L001_R1.fastq.merged.gz
PRO1563_S1_lib_CAGAGGTA-TATTAGAG_L001_R1.fastq.merged.gz Ab5_T_gDNA_PRO1563_S1_lib_CAGAGGTA-TATTAGAG_L001_R1.fastq.merged.gz
PRO1563_S1_lib_CTCCGGTA-GGTTCCTC_L001_R1.fastq.merged.gz Ab6_L_ATAC_PRO1563_S1_lib_CTCCGGTA-GGTTCCTC_L001_R1.fastq.merged.gz
PRO1563_S1_lib_CTCCGGTA-TGAGCAAC_L001_R1.fastq.merged.gz Ab6_L_gDNA_PRO1563_S1_lib_CTCCGGTA-TGAGCAAC_L001_R1.fastq.merged.gz
PRO1563_S1_lib_CTCTAGGT-GTTAGTAC_L001_R1.fastq.merged.gz Ab6_T_ATAC_PRO1563_S1_lib_CTCTAGGT-GTTAGTAC_L001_R1.fastq.merged.gz
PRO1563_S1_lib_CTCTAGGT-TGAGCAAC_L001_R1.fastq.merged.gz Ab6_T_gDNA_PRO1563_S1_lib_CTCTAGGT-TGAGCAAC_L001_R1.fastq.merged.gz
PRO1563_S1_lib_GCTACGCT-TCGTCGGC_L001_R1.fastq.merged.gz Pnm1_L_ATAC_PRO1563_S1_lib_GCTACGCT-TCGTCGGC_L001_R1.fastq.merged.gz
PRO1563_S1_lib_GTATTATC-AGCGTCGA_L001_R1.fastq.merged.gz Nb4_L_gDNA_PRO1563_S1_lib_GTATTATC-AGCGTCGA_L001_R1.fastq.merged.gz
PRO1563_S1_lib_GTCAACCG-ATATAACC_L001_R1.fastq.merged.gz Nb4_T_gDNA_PRO1563_S1_lib_GTCAACCG-ATATAACC_L001_R1.fastq.merged.gz
PRO1563_S1_lib_GTCAACCG-TCCTTCTT_L001_R1.fastq.merged.gz Nb4_L_ATAC_PRO1563_S1_lib_GTCAACCG-TCCTTCTT_L001_R1.fastq.merged.gz
PRO1563_S1_lib_GTCCTTGA-TCTATCCT_L001_R1.fastq.merged.gz Nb4_T_ATAC_PRO1563_S1_lib_GTCCTTGA-TCTATCCT_L001_R1.fastq.merged.gz
PRO1563_S1_lib_GTCGTGAT-TCGTCGGC_L001_R1.fastq.merged.gz Pnm1_T_ATAC_PRO1563_S1_lib_GTCGTGAT-TCGTCGGC_L001_R1.fastq.merged.gz
PRO1563_S1_lib_TAGGTTAG-AACCATGG_L001_R1.fastq.merged.gz Nb5_L_ATAC_PRO1563_S1_lib_TAGGTTAG-AACCATGG_L001_R1.fastq.merged.gz
PRO1563_S1_lib_TAGGTTAG-CCGTCAGC_L001_R1.fastq.merged.gz Nb5_L_gDNA_PRO1563_S1_lib_TAGGTTAG-CCGTCAGC_L001_R1.fastq.merged.gz
PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R1.fastq.merged.gz Nb5_T_ATAC_PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R1.fastq.merged.gz
PRO1563_S1_lib_TCGCCTGC-CCTTGCCG_L001_R1.fastq.merged.gz Nb5_T_gDNA_PRO1563_S1_lib_TCGCCTGC-CCTTGCCG_L001_R1.fastq.merged.gz
PRO1563_S1_lib_ACCACTGT-TCGTCGGC_L001_R2.fastq.merged.gz Pnm1_L_gDNA_PRO1563_S1_lib_ACCACTGT-TCGTCGGC_L001_R2.fastq.merged.gz
PRO1563_S1_lib_CAGAATGC-GAACTGAG_L001_R2.fastq.merged.gz Ab5_L_ATAC_PRO1563_S1_lib_CAGAATGC-GAACTGAG_L001_R2.fastq.merged.gz
PRO1563_S1_lib_CAGAATGC-TAACTCTA_L001_R2.fastq.merged.gz Ab5_L_gDNA_PRO1563_S1_lib_CAGAATGC-TAACTCTA_L001_R2.fastq.merged.gz
PRO1563_S1_lib_CAGAGAGG-TCGTCGGC_L001_R2.fastq.merged.gz Pnm1_T_gDNA_PRO1563_S1_lib_CAGAGAGG-TCGTCGGC_L001_R2.fastq.merged.gz
PRO1563_S1_lib_CAGAGGTA-GAAGGACG_L001_R2.fastq.merged.gz Ab5_T_ATAC_PRO1563_S1_lib_CAGAGGTA-GAAGGACG_L001_R2.fastq.merged.gz
PRO1563_S1_lib_CAGAGGTA-TATTAGAG_L001_R2.fastq.merged.gz Ab5_T_gDNA_PRO1563_S1_lib_CAGAGGTA-TATTAGAG_L001_R2.fastq.merged.gz
PRO1563_S1_lib_CTCCGGTA-GGTTCCTC_L001_R2.fastq.merged.gz Ab6_L_ATAC_PRO1563_S1_lib_CTCCGGTA-GGTTCCTC_L001_R2.fastq.merged.gz
PRO1563_S1_lib_CTCCGGTA-TGAGCAAC_L001_R2.fastq.merged.gz Ab6_L_gDNA_PRO1563_S1_lib_CTCCGGTA-TGAGCAAC_L001_R2.fastq.merged.gz
PRO1563_S1_lib_CTCTAGGT-GTTAGTAC_L001_R2.fastq.merged.gz Ab6_T_ATAC_PRO1563_S1_lib_CTCTAGGT-GTTAGTAC_L001_R2.fastq.merged.gz
PRO1563_S1_lib_CTCTAGGT-TGAGCAAC_L001_R2.fastq.merged.gz Ab6_T_gDNA_PRO1563_S1_lib_CTCTAGGT-TGAGCAAC_L001_R2.fastq.merged.gz
PRO1563_S1_lib_GCTACGCT-TCGTCGGC_L001_R2.fastq.merged.gz Pnm1_L_ATAC_PRO1563_S1_lib_GCTACGCT-TCGTCGGC_L001_R2.fastq.merged.gz
PRO1563_S1_lib_GTATTATC-AGCGTCGA_L001_R2.fastq.merged.gz Nb4_L_gDNA_PRO1563_S1_lib_GTATTATC-AGCGTCGA_L001_R2.fastq.merged.gz
PRO1563_S1_lib_GTCAACCG-ATATAACC_L001_R2.fastq.merged.gz Nb4_T_gDNA_PRO1563_S1_lib_GTCAACCG-ATATAACC_L001_R2.fastq.merged.gz
PRO1563_S1_lib_GTCAACCG-TCCTTCTT_L001_R2.fastq.merged.gz Nb4_L_ATAC_PRO1563_S1_lib_GTCAACCG-TCCTTCTT_L001_R2.fastq.merged.gz
PRO1563_S1_lib_GTCCTTGA-TCTATCCT_L001_R2.fastq.merged.gz Nb4_T_ATAC_PRO1563_S1_lib_GTCCTTGA-TCTATCCT_L001_R2.fastq.merged.gz
PRO1563_S1_lib_GTCGTGAT-TCGTCGGC_L001_R2.fastq.merged.gz Pnm1_T_ATAC_PRO1563_S1_lib_GTCGTGAT-TCGTCGGC_L001_R2.fastq.merged.gz
PRO1563_S1_lib_TAGGTTAG-AACCATGG_L001_R2.fastq.merged.gz Nb5_L_ATAC_PRO1563_S1_lib_TAGGTTAG-AACCATGG_L001_R2.fastq.merged.gz
PRO1563_S1_lib_TAGGTTAG-CCGTCAGC_L001_R2.fastq.merged.gz Nb5_L_gDNA_PRO1563_S1_lib_TAGGTTAG-CCGTCAGC_L001_R2.fastq.merged.gz
PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R2.fastq.merged.gz Nb5_T_ATAC_PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R2.fastq.merged.gz
PRO1563_S1_lib_TCGCCTGC-CCTTGCCG_L001_R2.fastq.merged.gz Nb5_T_gDNA_PRO1563_S1_lib_TCGCCTGC-CCTTGCCG_L001_R2.fastq.merged.gz

### 1. Create the appropriate directory structure and file paths for downstream analysis (and to run following scripts)

cd $WD

prefix=($scripts/prefix.txt)

echo '# -- 0. File merging complete -- #'

awk -F' ' '{print $2}' $libids | awk -F'_' '{print $1"_"$2"_"$3}' | sort -u > $prefix # create a prefix file to iterate

while IFS= read -r i
do
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
