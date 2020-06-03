#!/bin/sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 8000 # memory pool for all cores
#SBATCH -t 0-15:00 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

################################################################################################################

# ATAC-seq pipeline - Part 3
# June 2020: Tarang K. Mehta, Earlham Institute, Norwich, UK

################################################################################################################

# Script usage: ./ATAC_Bioinf_pipeline_v2c.sh

## Place this script and the following files in $WD (created in first script)
# 1. As used in './ATAC_Bioinf_pipeline_v2a.sh': a 2-column space-delimited table where col1='R1/R2 filename's col2='desired species renamed filename: species_tissue_experiment e.g. Mz_L_ATAC/gDNA'
# 2. Scripts:
  # ATAC_Bioinf_pipeline_v2b_part5bD-a.py
  # ATAC_Bioinf_pipeline_v2b_part5bD.py
# 3. Run as an sbatch script with 8Gb memory and ~3 days runtime - will spawn off other jobs

################################################################################################################

# ~ This pipeline is ran for all ATAC narrow peak files generated by 'ATAC_Bioinf_pipeline_v2b.sh', and contains the following components:

# 1. IDR on all pairs of replicates, self-pseudoreplicates and pooled pseudoreplicates - IDR is optional. The IDR peaks are a subset of the naive overlap peaks that pass a specific IDR threshold of 10%.
# 	1a. IDR of true replicates
# 	1b. Compute Fraction of Reads in Peaks (FRiP)
# 2. Create signal tracks - bedtools
# 3. Annotation:
# 	3a. TSS enrichment
# 	3b. Fraction of Reads in annotated regions

################################################################################################################

# All variables are added (and can be amended) here

scripts=(/tgac/workarea/group-vh/Tarang/ATACseq/2.run2) # place all scripts in the topmost directory - create this separately
# WD=(/tgac/workarea/group-vh/Tarang/ATACseq/2.run2/$spID) # insert the working directory
email=Tarang.Mehta@earlham.ac.uk # SBATCH out and err send to address

### 1. IDR
idrdir=($scripts/1.IDR) # assign raw reads dir
prefixATAC=($scripts/prefixATAC.txt)
prefixpairs=($idrdir/prefixpairs.txt)
idrpairs=($idrdir/idrpairpaths.txt)
idrpair1=($idrdir/idrpairpaths_rep1.txt)
idrpair2=($idrdir/idrpairpaths_rep2.txt)
idrprefix=($idrdir/idrprefix.txt)
idrpair1a=($idrdir/idrpairpaths_rep1unzipped.txt) # unzipped narrow peaks file path for rep1
idrpair2a=($idrdir/idrpairpaths_rep2unzipped.txt) # unzipped narrow peaks file path for rep2
IDR_THRESH=0.1 # consider changing to 0.05
IDR_THRESH_TRANSFORMED=$(awk -v p=${IDR_THRESH} 'BEGIN{print -log(p)/log(10)}')
npeaks_idr=($idrdir/idr_10thresh_peakcount.txt) # Number of peaks passing IDR thresholds of 10%
npeaks_idrrep1=($idrdir/idr_10thresh_peakcount_rep1.txt) # Number of peaks passing IDR thresholds of 10% for replicate 1
npeaks_idrrep2=($idrdir/idr_10thresh_peakcount_rep2.txt) # Number of peaks passing IDR thresholds of 10% for replicate 2


################################################################################################################

### 1. Irreproducible Discovery Rate (IDR) on MACS2 narrow peaks

## It is worth following the ENCODE project’s “ATAC-seq Data Standards and Prototype Processing Pipeline” for replicated data on MACS2 peak calling:
  # Peak call with MACS2 > narrow peaks file > IDR on true replicates e.g. Ab5_L and Ab6_L: This is the dataset to use
    # IDR is A statistical procedure called the Irreproducible Discovery Rate (IDR) operates on the replicated peak set and compares consistency of ranks of these peaks in individual replicate/pseudoreplicate peak sets.
      # checks the reproducibility information from the duplicates using the IDR statistic
      # The basic idea is that if two replicates measure the same underlying biology, the most significant peaks, which are likely to be genuine signals, are expected to have high consistency between replicates, whereas peaks with low significance, which are more likely to be noise, are expected to have low consistency.
      # If the consistency between a pair of rank lists (peaks) that contains both significant and insignificant findings is plotted, a transition in consistency is expected
      # This consistency transition provides an internal indicator of the change from signal to noise and suggests how many peaks have been reliably detected - red being false, black being true
      # By fitting a bivariate rank distribuion, IDR finds a threshold to separate real peaks from noise
    # install idr locally: https://github.com/kundajelab/idr
      # wget https://github.com/kundajelab/idr/archive/2.0.4.zip
      # unzip 2.0.4.zip
      # cd idr-2.0.4/
      # ml gcc
      # ml zlib
      # python3 setup.py install
    # follow IDR details: https://docs.google.com/document/d/1f0Cm4vRyDQDu0bMehHD7P7KOMxTOP-HiNoIvL1VcBt8/edit#
      # set IDR_THRESH=0.1
    # If you have more than 2 true replicates select the longest peak list from all pairs that passes the IDR threshold.

## NOTE:
# 1. run IDR as an array for all replicate pairs
# 2. In cases where there are more than 2 true replicates select the longest peak list from all pairs that passes the IDR threshold.

mkdir -p $idrdir
cd $idrdir

# 1a. create a 2-column tab delimited file that has all narrowpeak file paths of pairs to compare
# 1b. split the above file into two files, one for each column
# 1c. using an array by iterating through each line of the two 1-column files from above, do the following for each pair:
  # 1cA. create a pooled-replicate narrowPeak file
  # 1cB. run IDR on each pair and get peaks passing threshold
  # 1cC. for cases where there are more than 2 true replicates, select the longest peak list from all pairs that passes the IDR threshold.

# =============================
# 1a. Create IDR comparison pairs
# create a 2-column tab delimited file that has all narrowpeak file paths of pairs to compare
# =============================

for a in $(awk '{print $1}' $prefixATAC); do
  for b in $(awk '{print $1}' $prefixATAC); do
    # echo -e "$a\t$b"
      echo "$a\t$b" | awk '{if($1 != $2) print $1, $2;}' OFS='\t' | awk -F"_" '{if($2==$4) print $0}' |
      awk '{if (substr($1,1,2)==substr($2,1,2)) {print $0, "YES"} else if (substr($1,1,2)!=substr($2,1,2)) {print $0, "NO"}}' OFS='\t' |
      grep -v 'NO' | cut -f1,2 >> $prefixpairs.temp
    done
done

awk -F'\t' '!seen[$1>$2 ? $1 FS $2 : $2 FS $1]++' $prefixpairs.temp > $prefixpairs # this removes duplicate pairs that are simply in a different order
rm $prefixpairs.temp

awk -F'\t' '{print "/tgac/workarea/group-vh/Tarang/ATACseq/2.run2/"$1"/5.peak_calling/"$1"_peaks.narrowPeak.gz","\t","/tgac/workarea/group-vh/Tarang/ATACseq/2.run2/"$2"/5.peak_calling/"$2"_peaks.narrowPeak.gz"}' $prefixpairs > $idrpairs

## AT THIS POINT, CONSIDER DOING A FILE CHECK OF ALL PAIRS AND ONLY CONTINUE WITH THOSE THAT ARE FOUND?


# =============================
# 1b. Separate the pairwise comparisons for an array
# split the above file into two files, one for each column
# =============================

cut -f1 $idrpairs > $idrpair1
cut -f2 $idrpairs > $idrpair2
paircount=$(wc -l $idrpairs | awk -F' ' '{print $1}') # assign variable for total number of pairs
IDRarray=0-$(expr $paircount - 1) # number of pairs for the array starting from 0

# =============================
# 1c. Run the IDR analyses by iterating in an array
# 1cA. unzip files and create a pooled-replicate narrowPeak file
# 1cB. Perform IDR analysis.
  # Generate a plot and IDR output with additional columns including IDR scores.
# 1cC. Get peaks passing IDR threshold of 10%

while IFS= read -r i; do
  gunzip $i
done < $idrpair1

while IFS= read -r i; do
  gunzip $i
done < $idrpair2

sed 's/.gz//g' $idrpair1 > $idrpair1a
sed 's/.gz//g' $idrpair2 > $idrpair2a


echo '#!/bin/bash -e' > 1c.IDR.sh
echo '#SBATCH -p tgac-short # partition (queue)' >> 1c.IDR.sh
echo '#SBATCH -N 1 # number of nodes' >> 1c.IDR.sh
echo '#SBATCH -n 1 # number of tasks' >> 1c.IDR.sh
echo "#SBATCH --array=$IDRarray" >> 1c.IDR.sh
echo '#SBATCH --mem-per-cpu 8000' >> 1c.IDR.sh
echo '#SBATCH -t 0-00:45' >> 1c.IDR.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 1c.IDR.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 1c.IDR.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 1c.IDR.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 1c.IDR.sh
printf '\n' >> 1c.IDR.sh
echo '# 1c. Run the IDR analyses by iterating in an array' >> 1c.IDR.sh
echo '# 1cA. unzip files and create a pooled-replicate narrowPeak file' >> 1c.IDR.sh
echo "mapfile -t idrrep1 < $idrpair1a" >> 1c.IDR.sh
echo "mapfile -t idrrep2 < $idrpair2a" >> 1c.IDR.sh
echo "awk '{print"' $1"_"$2}'"' $prefixpairs > $idrprefix" >> 1c.IDR.sh
echo "mapfile -t idrprefix < $idrprefix" >> 1c.IDR.sh
printf '\n' >> 1c.IDR.sh
echo 'cat ${idrrep1[${SLURM_ARRAY_TASK_ID}]} ${idrrep2[${SLURM_ARRAY_TASK_ID}]} > ${idrprefix[${SLURM_ARRAY_TASK_ID}]}_peaks.narrowPeak'
echo '# 1cB. Perform IDR analysis.' >> 1c.IDR.sh
echo '# Generate a plot and IDR output with additional columns including IDR scores.' >> 1c.IDR.sh
echo 'srun idr --samples ${idrrep1[${SLURM_ARRAY_TASK_ID}]} ${idrrep2[${SLURM_ARRAY_TASK_ID}]} --peak-list ${idrprefix[${SLURM_ARRAY_TASK_ID}]}_peaks.narrowPeak --input-file-type narrowPeak --output-file ${idrprefix[${SLURM_ARRAY_TASK_ID}]}.IDR0.1output --rank p.value --soft-idr-threshold '"${IDR_THRESH} --plot --use-best-multisummit-IDR" >> 1c.IDR.sh
echo '# 1cC. Get peaks passing IDR threshold of 10%' >> 1c.IDR.sh
echo '# IDR QC to report and using the IDR output' >> 1c.IDR.sh
echo '# 1. For each biological replicate pair, filter the IDR peaks based on the ${IDR_THRESH_TRANSFORMED} and sort descending based on signal.value (col7)' >> 1c.IDR.sh
printf '\n' >> 1c.IDR.sh
echo "awk 'BEGIN{OFS="'"\t"} $12>='"'"'"${IDR_THRESH_TRANSFORMED}"'"' {print "'$1,$2,$3,$4,$5,$6,$7,$8,$9,$10}'"' "'${idrprefix[${SLURM_ARRAY_TASK_ID}]}.IDR0.1output | sort | uniq | sort -k7,7rn > ${idrprefix[${SLURM_ARRAY_TASK_ID}]}.IDR0.1Transf.narrowPeak' >> 1c.IDR.sh
printf '\n' >> 1c.IDR.sh
echo '# 2. For each biological replicate pair, count the number of lines passing ${IDR_THRESH_TRANSFORMED} - this is where IDR finds the threshold to separate real peaks from noise' >> 1c.IDR.sh
echo 'wc -l ${idrprefix[${SLURM_ARRAY_TASK_ID}]}.IDR0.1Trans.narrowPeak >>'" $npeaks_idr # Number of peaks passing IDR thresholds of 10%" >> 1c.IDR.sh
echo '# 3. For each biological replicate in a pair, assign the wc -l as max_numPeaks_Rep peaks' >> 1c.IDR.sh
awk -F'.' '{print $1}' $npeaks_idr | awk -F'_' '{print $1"_"$2"_"$3" "$4"_"$5"_"$6}' | sed 's/ /\t/g' > $npeaks_idr.tmp
echo "cut -f1,2 $npeaks_idr.tmp > $npeaks_idrrep1" >> 1c.IDR.sh
echo "cut -f1,3 $npeaks_idr.tmp > $npeaks_idrrep2" >> 1c.IDR.sh
echo '# cat $npeaks_idr.tmp2 $npeaks_idr.tmp3 | sort -k2,2 > $npeaks_idr2' >> 1c.IDR.sh
echo "rm *tmp*" >> 1c.IDR.sh
printf '\n' >> 1c.IDR.sh
echo "# 4. For each biological replicate pair, sort the MACS2 narrowPeak file based on signal.value column (col7) and add another column on end 'max_numPeaks_Rep', adding 'T' for True and 'F' for False of peaks that pass IDR, based on wc -l" >> 1c.IDR.sh
printf '\n' >> 1c.IDR.sh
echo '# mapfile -t idrrep1 < $idrpair1a' >> 1c.IDR.sh
echo '# mapfile -t idrrep2 < $idrpair2a' >> 1c.IDR.sh
echo '# mapfile -t idrrep1peak < $npeaks_idrrep1' >> 1c.IDR.sh
echo '# mapfile -t idrrep2peak < $npeaks_idrrep2' >> 1c.IDR.sh


echo '#!/bin/bash -e' > 1d.IDR.sh
echo '#SBATCH -p tgac-medium # partition (queue)' >> 1d.IDR.sh
echo '#SBATCH -N 1 # number of nodes' >> 1d.IDR.sh
echo '#SBATCH -n 1 # number of tasks' >> 1d.IDR.sh
echo '#SBATCH --mem-per-cpu 8000' >> 1d.IDR.sh
echo '#SBATCH -t 0-02:45' >> 1d.IDR.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 1d.IDR.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 1d.IDR.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 1d.IDR.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 1d.IDR.sh
printf '\n' >> 1d.IDR.sh
echo '# Double while read array:' >> 1d.IDR.sh
echo '# 1. Read in the peak numbers and narrowPeak files and do an IF statement to match the corresponding files' >> 1d.IDR.sh
echo '# 2. sort peaks file based on col7 (signal value - fold-change at peak summit)' >> 1d.IDR.sh
echo '# 3. read from file paths, open; read from peaks' >> 1d.IDR.sh
echo '# 4. use number in col1 to add T to that line and all preceding lines and add F for all lines thereafter' >> 1d.IDR.sh
printf '\n' >> 1d.IDR.sh
echo 'while IFS= read -u 3 -r reppeak && IFS= read -u 4 -r npeak; do' >> 1d.IDR.sh
echo -e '\t# echo ${npeak}' >> 1d.IDR.sh
echo -e '\t# echo ${reppeak}' >> 1d.IDR.sh
echo -e '\tpeakline=$(echo ${npeak} | awk -F'"' ' '{print "'$1}'"') # get the total count of IDR passed peaks" >> 1d.IDR.sh
echo -e '\tpeakprefix=$(echo ${npeak} | awk -F'"' ' '{print "'$2}'"') # get the file prefix this corresponds to" >> 1d.IDR.sh
echo -e '\t# echo $peakline' >> 1d.IDR.sh
echo -e '\t# echo $peakprefix' >> 1d.IDR.sh
echo -e '\tout=$(echo $(basename ${reppeak}) | sed -e '"'s/.narrowPeak/.final.narrowPeak/') # create an outfile" >> 1d.IDR.sh
echo -e '\tinprefix=$(echo $(basename ${reppeak}) | sed -e '"'s/_peaks.narrowPeak//') # get the prefix of the narrow peak in file" >> 1d.IDR.sh
echo -e '\t# echo $out' >> 1d.IDR.sh
echo -e '\t# echo $inprefix' >> 1d.IDR.sh
echo -e '\tif [[ $peakprefix = $inprefix ]]; then # check that the peaks passing IDR and narrowPeak file match' >> 1d.IDR.sh
echo -e '\t\techo "peakcount>> $peakprefix = $inprefix <<narrowPeaks_file: npeaks and peaks file ARE matched"'
echo -e '\t\tsort -k7,7rn ${reppeak} | awk -v peakline="$peakline" '"'{if(NR>=1 && NR<=peakline)print "'$0,"T";else print $0,"F";}'"' OFS='\t' > "'${out} # check with awk '"'FNR>=103003 && FNR<=103006'" >> 1d.IDR.sh
echo -e '\telse' >> 1d.IDR.sh
echo -e '\t\techo "peakcount>> $peakprefix != $inprefix <<narrowPeaks_file: npeaks and peaks file NOT matched"' >> 1d.IDR.sh
echo -e '\tfi' >> 1d.IDR.sh
echo "done 3<$idrpair1a 4<$npeaks_idrrep1" >> 1d.IDR.sh


# IDR output
# Broad peak output files are the same except that they do not include the the summit columns (e.g. columns 10, 18, and 22 for samples with 2 replicates)
#
#     1. chrom string
#     Name of the chromosome for common peaks
#
#     2. chromStart int
#     The starting position of the feature in the chromosome or scaffold for common peaks, shifted based on offset. The first base in a chromosome is numbered 0.
#
#     3. chromEnd int
#     The ending position of the feature in the chromosome or scaffold for common peaks. The chromEnd base is not included in the display of the feature.
#
#     4. name string
#     Name given to a region (preferably unique) for common peaks. Use '.' if no name is assigned.
#
#     5. score int
#     Contains the scaled IDR value, min(int(log2(-125IDR), 1000). e.g. peaks with an IDR of 0 have a score of 1000, idr 0.05 have a score of int(-125log2(0.05)) = 540, and idr 1.0 has a score of 0.
#
#     6. strand [+-.] Use '.' if no strand is assigned.
#
#     7. signalValue float
#     Measurement of enrichment for the region for merged peaks. When a peak list is provided this is the value from the peak list.
#
#     8. p-value float
#     Merged peak p-value. When a peak list is provided this is the value from the peak list.
#
#     9. q-value float
#     Merged peak q-value. When a peak list is provided this is the value from the peak list.
#
#     10. summit int
#     Merged peak summit
#
#     11. localIDR float -log10(Local IDR value)
#
#     12. globalIDR float -log10(Global IDR value)
#
#     13. rep1_chromStart int
#     The starting position of the feature in the chromosome or scaffold for common replicate 1 peaks, shifted based on offset. The first base in a chromosome is numbered 0.
#
#     14. rep1_chromEnd int
#     The ending position of the feature in the chromosome or scaffold for common replicate 1 peaks. The chromEnd base is not included in the display of the feature.
#
#     15. rep1_signalValue float
#     Signal measure from replicate 1. Note that this is determined by the --rank option. e.g. if --rank is set to signal.value, this corresponds to the 7th column of the narrowPeak, whereas if it is set to p.value it corresponds to the 8th column.
#
#     16. rep1_summit int
#     The summit of this peak in replicate 1.
#
# [rep 2 data]
#
# ...
#
# [rep N data]


# The plot (*.IDR0.1output.png) for each quadrant is described below:
# Upper Left: Replicate 1 peak ranks versus Replicate 2 peak ranks - peaks that do not pass the specified idr threshold are colored red.
# Upper Right: Replicate 1 log10 peak scores versus Replicate 2 log10 peak scores - peaks that do not pass the specified idr threshold are colored red.
# Bottom Row: Peak rank versus IDR scores are plotted in black. The overlayed boxplots display the distribution of idr values in each 5% quantile. The IDR values are thresholded at the optimization precision - 1e-6 by default.

## FINAL PEAK FILES based on IDR QC are *.final.narrowPeak - the last column of T (True) or F (False) show peaks passing IDR

# =============================

echo '# -- 1. IDR started -- #'

JOBID1=$( sbatch -W --array=$IDRarray 1c.IDR.sh | awk '{print $4}' ) # Run the first job and then store the first job to variable JOBID1 (taken by awk once run)
JOBID2=$( sbatch -W --dependency=afterok:${JOBID1} 1d.IDR.sh | awk '{print $4}' ) # JOB2 depends on JOB1 completing successfully

################################################################################################################

### 2. Create signal tracks - normally bedtools

# Signal tracks are generated from BAM file (Raw) and bias corrected by HINT-ATAC (Bias corrected)
# This could be rolled in with TF footprinting using HINT-ATAC, see this: https://www.regulatory-genomics.org/hint/tutorial/

# install HINT-ATAC by installing RGT - Regulatory Genomics Toolbox: https://github.com/CostaLab/reg-gen

source rgt-0.12.3

# CHECK THAT THE INPUT OF OUTPUT FROM ABOVE IS OK

# 1. Customise RGT data folder and data.config file for own genome files etc.: http://www.regulatory-genomics.org/rgt/rgt-data-folder/
  # Please run download-db.sh to download all required RGT database files to /opt/software/conda_env/share/rgt-0.12.3/db/\n\n
  # need to find the rgt-data folder!!
# 2. follow scripts here to run: https://www.regulatory-genomics.org/hint/tutorial/


## ~ INSERT CODE HERE ~ ##

echo '# -- 1. IDR has completed -- #'

echo '# -- 2. Creating signal tracks has started -- #'

JOBID3=$( sbatch -W --dependency=afterok:${JOBID2} XX.sh | awk '{print $4}' ) # JOB2 depends on JOB1 completing successfully


################################################################################################################

### 3. Annotation:
# 	3a. TSS enrichment - plot
  # The TSS enrichment calculation is a signal to noise calculation.
  # The reads around a reference set of TSSs are collected to form an aggregate distribution of reads centered on the TSSs and extending to 1000 bp in either direction (for a total of 2000bp).
  # This distribution is then normalized by taking the average read depth in the 100 bps at each of the end flanks of the distribution (for a total of 200bp of averaged data) and calculating a fold change at each position over that average read depth.
  # This means that the flanks should start at 1, and if there is high read signal at transcription start sites (highly open regions of the genome) there should be an increase in signal up to a peak in the middle.
  # We take the signal value at the center of the distribution after this normalization as our TSS enrichment metric.

# 	3b. Fraction of Reads in annotated regions

## ~ INSERT CODE HERE ~ ##

echo '# 3a. TSS enrichment plotting' >> 5.peakcall.sh
echo '# This calls two python scripts - make sure they are in $scripts' >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo '# create a 2kb window around TSS (+/- 1kb) bed file e.g.' >> 5.peakcall.sh
echo '# chr1	134210701	134214701	+' >> 5.peakcall.sh
echo '# chr1	33724603	33728603	-' >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo '# 3aA. Protein-coding genes GTF > BED: variable $annot of gzipped GTF file (*gtf.gz); or 2) longest protein-coding gene annotations as 6-column BED: col1-scaff,col2-start,col3-end,col4-geneID,col5-XX,col6-strand' >> 5.peakcall.sh
echo '# output is genebed=($peakcall/$spG'_refGene.bed') defined as variable at top' >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo 'source bedops-2.4.28' >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo '# the below is for either one of your own gtf or bed files' >> 5.peakcall.sh
echo "case $annot in" >> 5.peakcall.sh
echo -e "\t*gtf.gz) gunzip -c $annot | awk 'OFS="'"\\t" {if ($3=="gene") {print $1,$4-1,$5,$10,0,$7,$18}}'"' | tr -d '"'";'"' | grep 'protein_coding' | awk '{print "'$1,$2,$3,$4,$5,$6}'"' OFS="'"\\t" > '"$genebed ;; # GTF > 0-based BED of protein_coding" >> 5.peakcall.sh
echo -e "\t*.bed) awk '{"'print $1,$2,$3,$4,0,$6}'"' OFS="'"\\t"'" $annot > $genebed ;; # BED format ONLY for my files that have six cols (where 6th col is strand)" >> 5.peakcall.sh
echo "esac" >> 5.peakcall.sh
# below is for the old annotations that are a little odd
# echo "case $annot in" >> 5.peakcall.sh
# echo -e "\t*gtf.gz) gunzip -c $annot | awk 'OFS="'"\\t" {if ($3=="gene" || $3=="exon") {print $1,$4-1,$5,$10,0,$7,$18}}'"' | tr -d '"'";'"' | awk '{print "'$1,$2,$3,$4,$5,$6}'"' OFS="'"\\t" > '"$genebed ;; # GTF > 0-based BED of protein_coding" >> 5.peakcall.sh
# echo -e "\t*.bed) awk '{"'print $1,$2,$3,$4,0,$6}'"' OFS="'"\\t"'" $annot > $genebed ;; # BED format ONLY for my files that have six cols (where 6th col is strand)" >> 5.peakcall.sh
# echo "esac" >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
# # below is the same as what is echoed above
# case $annot in
#   *gtf.gz) gunzip -c $annot | awk 'OFS="\t" {if ($3=="gene" || $3=="exon") {print $1,$4-1,$5,$10,0,$7,$18}}' | tr -d '";' | awk '{print $1,$2,$3,$4,$5,$6}' OFS='\t' > $genebed ;; # GTF > 0-based BED of protein_coding
#   *.bed) awk '{print $1,$2,$3,$4,0,$6}' OFS='\t' $annot > $genebed ;; # BED format ONLY for my files that have six cols (where 6th col is strand)
# esac
# # the below is for either a gtf file from ensembl etc. and your own bed file
# case $annot in
#   *gtf.gz) gunzip -c $annot | awk 'OFS="\t" {if ($3=="gene" || $3=="exon") {print $1,$4-1,$5,$10,0,$7,$18}}' | tr -d '";' | grep -wiF 'protein_coding' | awk '{print $1,$2,$3,$4,$5,$6}' OFS='\t' > $genebed ;; # GTF > 0-based BED of protein_coding
#   *.bed) awk '{print $1,$2,$3,$4,0,$6}' OFS='\t' $annot > $genebed ;; # BED format ONLY for my files that have six cols (where 6th col is strand)
# esac
echo '# 3aB. Then split them by strand and pad around the stranded-start position of the annotation (taking TSS +/- 1000=1kb)' >> 5.peakcall.sh
echo "awk '("'$6 == "+") { print $0 }'"' $genebed | awk 'BEGIN{ OFS="'"\t" }($2 > 1000){ print $1, ($2 - 1000), ($2 + 1000), $4, $5, $6  }'"' > $genebed.tss.for.padded.bed" >> 5.peakcall.sh
echo "awk '("'$6 == "-") { print $0 }'"' $genebed | awk 'BEGIN{ OFS="'"\t" }($3 > 1000){ print $1, ($3 - 1000), ($3 + 1000), $4, $5, $6  }'"' > $genebed.tss.rev.padded.bed" >> 5.peakcall.sh
echo "bedops --everything $genebed.tss.for.padded.bed $genebed.tss.rev.padded.bed > $genebedtss" >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo '# 3aC. Keep only TSS regions within chromosomal bounds - prep scaffold sizes file (col1=scaffoldID; col2=0; col3=length) from genome fasta' >> 5.peakcall.sh
echo 'bioawk -c fastx'" '{ print "'$name, length($seq) }'"' < $gFA | awk '{print "'$1,"0",$2}'"' OFS="'"\t" > '"$scafflen" >> 5.peakcall.sh
echo "bedops --element-of 100% $genebedtss $scafflen > $genebedtss2" >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
# echo '# 5bD. Use final TSS (+/- 1kb) bed file as input to calculate TSS enrichment and plot with python script ATAC_Bioinf_pipeline_v2b_part5bD.py' >> 5.peakcall.sh
# echo "python3 $scripts/ATAC_Bioinf_pipeline_v2b_part5bD-a.py $fastqr1 # input fastq can be native or gzipped" >> 5.peakcall.sh
# echo "python3 $scripts/ATAC_Bioinf_pipeline_v2b_part5bD.py $Test1 $genebedtss2 $spID $read_len $scafflen"' # usage: python3 $scripts/ATAC_Bioinf_pipeline_v2b_part5bD.py'" 'FINAL_BAM' 'TSS' 'OUTPUT_PREFIX' 'read_len' 'CHROMSIZES'" >> 5.peakcall.sh
# printf '\n' >> 5.peakcall.sh


echo '# -- 2. Creating signal tracks has completed -- #'

echo '# -- 3. Peak annotation has started -- #'

JOBID3=$( sbatch -W --dependency=afterok:${JOBID2} XX.sh | awk '{print $4}' ) # JOB3 depends on JOB2 completing successfully


################################################################################################################

### 4. Differential analysis of peaks

### NOTE: be careful when considering differential peaks as some may be only offset by a few bases. In this secnario, consider the average number of mapped reads over a window.

## FOR SIMPLE DIFFERENTIAL ANALYSIS OF PEAKS, USE HOMER; SOME CODE HERE: https://dtc-coding-dojo.github.io/main//blog/Analysing_ATAC_and_CHIPseq_data/
## Then use DiffBind to:
  # A. Determine tissue-specific peaks in each species
    # Tissue-specificity of ATAC-seq peaks was determined using DiffBind (https://www.bioconductor.org/packages/release/bioc/ html/DiffBind.html),
    # This provided the peak coordinates for each of the biological replicates of all tissues profiled as input, plus the mapped and shifted sequencing reads (parameters ‘method = DBA_EDGER, bFullLibrarySize = FALSE, bSubControl = FALSE, bTagwise = FALSE’).
    # All peaks identified with a log2 fold change equal or greater than 1 in one tissue compared to all others were selected as tissue-specific.
  # B. Determine tissue-specific peaks between species same tissues e.g. Ab5_L vs Nb5_L
    # you need to find a way to compare different species - use association to orthologous genes?


## ~ INSERT CODE HERE ~ ##

echo '# -- 3. Peak annotation has completed -- #'

echo '# -- 4. Differential analysis of peaks has started -- #'

JOBID4=$( sbatch -W --dependency=afterok:${JOBID3} XX.sh | awk '{print $4}' ) # JOB4 depends on JOB3 completing successfully

################################################################################################################

echo '# -- 4. Differential analysis of peaks has completed -- #'

### Finish the script
exit 0
