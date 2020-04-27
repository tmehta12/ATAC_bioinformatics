#!/bin/sh

################################################################################################################

# ATAC-seq pipeline - Part 2
# March 2020: Tarang K. Mehta, Earlham Institute, Norwich, UK

################################################################################################################

# Script usage: ./ATAC_Bioinf_pipeline_v2b.sh -s "spID" -g "spG" -f "gFA" -m "mtID" -u "Usr" -a "annot"
# e.g. ./ATAC_Bioinf_pipeline_v2b.sh -s Mz1_L_ATAC -g M_zebra_UMD1 -f /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -m KT221043 -u mehtat
# Note: Script is adapted for SBATCH usage

## Place the following files in $scripts - to be created prior to running
# 1. Prepare a 2-column space-delimited table where col1='R1/R2 filename's col2='desired species renamed filename'
# 2. Scripts:
  # ATAC_Bioinf_pipeline_v2b_part3a.py
  # ATAC_Bioinf_pipeline_v2b_part3b.R
  # ATAC_Bioinf_pipeline_v2b_part5bD.py
# 3. Run as an sbatch script with 8Gb memory and >15 day runtime

################################################################################################################

# ~ This pipeline follows on from 'ATAC_Bioinf_pipeline_v2a.sh', that:

# 0. Merges files sequenced over multiple lanes

#### NOTE: TO DO (as of 08/04/2020) - create a separate script that is only part 1 and part2 for gDNA mapping
## Think about how the gDNA control is called:
  # 1. There needs to be a separate pipeline for gDNA controls that runs until peak calling e.g no mitochondrial reads removed - this script can be adapted at the end


################################################################################################################

# ~ This pipeline is ran as species-specific, and contains the following components:

# 1. Trim adaptors - trimgalore
# 2. Read alignment - bowtie2 > samtools sorted bam
# 3. Remove mitochondrial mapped - samtools > new bam
# 	3a. Fragment/insert size distribution
# 4. Post alignment filtering - Sort, Map and remove duplicates, Remove reads unmapped, not primary alignment, reads failing platform, duplicates: samtools > final bam
# 	4a. Convert PE Bam to tagalign - bedtools
# 	4b. Calculate Cross-correlation QC scores - phantompeakqualtools (v1.2.1)
# 	4c. Generate self-pseudoreplicates for each replicate - probably not required
# 	4d. Generate pooled dataset and pooled-pseudoreplicates
# 	4e. TN5 shifting of tagaligns for ATAC Seq - gawk
# 	4f. Calculate Jensen-Shannon distance (JSD) - deeptools (v3.3.0) plotFingerprint
# 	4g. Calculate GC bias
# 	4h. Fragment length statistics (PE only)
# 5. Call peaks on replicates, self-pseudoreplicates, pooled data and pooled-pseudoreplicates
# 	5a. peak calling - macs2
# 	5b. Blacklist filtering for peaks - bedtools
# 	5c. Bed to bigbed conversion for narrowpeaks
# 	5d. Naive overlap thresholding for MACS2 peak calls
# 6. IDR on all pairs of replicates, self-pseudoreplicates and pooled pseudoreplicates - IDR is optional. The IDR peaks are a subset of the naive overlap peaks that pass a specific IDR threshold of 10%.
# 	6a. IDR of true replicates
# 	6b. Compute Fraction of Reads in Peaks (FRiP)
# 7. Create signal tracks - bedtools
# 8. Annotation:
# 	8a. TSS enrichment
# 	8b. Fraction of Reads in annotated regions

## This script will spawn off several other jobs and thus, can be ran with 12Gb memory (to create directories etc.) but with a long run time

################################################################################################################

# Setting parameters for command line input

helpFunction()
{
   echo ""
   echo "Usage: $0 -s spID -g spG -f gFA -m mtID -u Usr"
   echo -e "\t-s spID = Species ID, preferably two short letters, individual and tissue e.g. Metriaclima zebra Individual 1 Liver ATAC/gDNA = Mz1_L_ATAC/Mz1_L_gDNA"
   echo -e "\t-g spG = Species genome ID e.g. hg19 or M_zebra_UMD1"
   echo -e "\t-f gFA = Full path to genome assembly in FASTA format"
   echo -e "\t-m mtID = NCBI accession number to complete mitochondrial genome FASTA"
   echo -e "\t-u Usr = ssh username e.g. mehtat"
   echo -e "\t-a annot = Full path to genome annotations: either 1) gzipped GTF file (*gtf.gz); or 2) longest protein-coding gene annotations as 6-column BED: col1-scaff,col2-start,col3-end,col4-geneID,col5-XX,col6-strand"
   exit 1 # Exit script after printing help
}

while getopts "s:g:f:m:u:a:" opt
do
   case "$opt" in
      s ) spID="$OPTARG" ;;
      g ) spG="$OPTARG" ;;
      f ) gFA="$OPTARG" ;;
      m ) mtID="$OPTARG" ;;
      u ) Usr="$OPTARG" ;;
      a ) annot="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$spID" ] || [ -z "$spG" ] || [ -z "$gFA" ] || [ -z "$mtID" ] || [ -z "$Usr" ] || [ -z "$annot" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "$spID"
echo "$spG"
echo "$gFA"
echo "$mtID"
echo "$Usr"
echo "$annot"

################################################################################################################

# All variables are added (and can be amended) here

scripts=(/tgac/workarea/group-vh/Tarang/ATACseq/2.run2) # place all scripts in the topmost directory - create this separately
WD=(/tgac/workarea/group-vh/Tarang/ATACseq/2.run2/$spID) # insert the working directory
email=Tarang.Mehta@earlham.ac.uk # SBATCH out and err send to address

### 1. Trim adaptors
rawreaddir=($WD/0.rawreads) # assign raw reads dir
trimdir=($WD/1.adaptor_trimming) # assign trimmed reads dir
trimarray=X-X # INSERT the number range of paired *fastq.merged.gz to trim in zero base e.g. 10 pairs = 0-9
libids=($trimdir/libids.txt) # 2-col space-delimited file where col1 is the *_{R1,R2}.fastq.merged.gz and col2 is the species_tissueID_barcode_{R1,R2}.fastq.merged.gz
libids2=$libids.sh

### 2. Read alignment
readalign=($WD/2.read_alignment)
RAarray=X-X # INSERT the number range of paired Mz *fastq.merged.gz to align in zero base e.g. 10 pairs = 0-9
idx=$spG # species bowtie index
bam='.bam' # output BAM to add prefix beforehand
log='.align.log' # output bowtie alignment log to add prefix beforehand
read1=($readalign/$spID'read1.txt')
read2=($readalign/$spID'read2.txt')
prefix=($readalign/$spID'prefix.txt')
reads=(reads.txt)
multimapping=4 # bowtie alignments reported for multimapping  (4 by default for ENCODE3)
bwt_thread=32 # bowtie parallel search threads and as SBATCH cpus-per task (if running array)
fgQC1=(_flagstat_qc1.txt) # mapping stats from flagstat (SAMstats)

### 3. Remove mitochondrial mapped and plot fragment distribution
mtfilt=($WD/3.Mtfilt_fragcnt) # mitochondrial filtering and frag count directory
mtgen=($mtfilt/$mtID'.fasta') # species mitochondrial genome
blstout=($spID.genome_mt.blast) # blast outfile
mtscaff=($spID.genome_mt.filtered.blast) # filtered BLAST to mitochondrial genome output

### 4. Post alignment filtering
filtdir=($WD/4.postalign_filt) # post alignment filtering directory

### 5. ATAC peak calling (test vs control)
peakcall=($WD/5.peak_calling) # ATAC peak calling directory
genebed=($peakcall/$spG'_refGene.bed') # GTF (or abnormal BED) > BED file of protein coding genes only
genebedtss=($peakcall/$spG'_refGene.tss.padded.bed') # BED file of protein coding genes with TSS +/- 1kb
scafflen=($peakcall/$spG'_scaffbounds.bed') # scaffold length boundaries of input genome assembly (col1=scaff; col2=0; col3=scaffold_length)
genebedtss2=($peakcall/$spG'_refGene.tss.padded.filt.bed') # BED file of protein coding genes with TSS +/- 1kb and filtered for any out of bound genes
fastqr1=($rawreaddir/*_R1.fastq.merged.gz)
read_len=($peakcall/*_read_length.txt)

################################################################################################################

### 1. Trim adaptors
mkdir -p $trimdir
cd $trimdir

# run as an array in sbatch script

echo '#!/bin/bash -e' > 1a.trimadaptors.sh
echo '#SBATCH -p tgac-medium # partition (queue)' >> 1a.trimadaptors.sh
echo '#SBATCH -N 1 # number of nodes' >> 1a.trimadaptors.sh
echo '#SBATCH -n 1 # number of tasks' >> 1a.trimadaptors.sh
echo "#SBATCH --array=$trimarray" >> 1a.trimadaptors.sh
echo '#SBATCH --mem-per-cpu 24000' >> 1a.trimadaptors.sh
echo '#SBATCH -t 0-14:59' >> 1a.trimadaptors.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 1a.trimadaptors.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 1a.trimadaptors.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 1a.trimadaptors.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 1a.trimadaptors.sh
printf '\n' >> 1a.trimadaptors.sh
echo 'source trim_galore-0.5.0' >> 1a.trimadaptors.sh
printf '\n' >> 1a.trimadaptors.sh
echo "for read1 in $rawreaddir/*_R1.fastq.merged.gz; do read2="'$(echo $read1 | sed'" 's/_R1.fastq.merged.gz/_R2.fastq.merged.gz/'); "'echo $read1 >> trimarrayread1; echo $read2 >> trimarrayread2; done' >> 1a.trimadaptors.sh
printf '\n' >> 1a.trimadaptors.sh
echo 'mapfile -t read1 < trimarrayread1 # assign read1 as elements to $read1 variable' >> 1a.trimadaptors.sh
echo 'mapfile -t read2 < trimarrayread2 # assign read2 as elements to $read2 variable' >> 1a.trimadaptors.sh
printf '\n' >> 1a.trimadaptors.sh
echo 'srun trim_galore --output_dir $trimdir --paired --fastqc ${read1[${SLURM_ARRAY_TASK_ID}]} ${read2[${SLURM_ARRAY_TASK_ID}]}' >> 1a.trimadaptors.sh

# assign '1a.trimadaptors.sh' to variable: JOBID1 and run
JOBID1=$( sbatch --hold --job-nam=START --array=$trimarray 1a.trimadaptors.sh | awk '{print $4}' ) # Run the first job and then store the first job to variable JOBID1 (taken by awk once run)

echo '# -- 1a. Adaptor trimming started -- #'

# rename the files according to species and tissue - provide this as a 2-column SPACE-delimited table that will be assigned as a variable above, placed in the trimdir
# only provide for the two paired files you are working on and not all files

# Example is (note, these are made up examples)
# echo 'PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R1.fastq.merged.gz Pn1T_TCGCCTGC-AACCGCCA_L001_R1.fastq.merged.gz' > libids.txt
# echo 'PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R2.fastq.merged.gz Pn1T_TCGCCTGC-AACCGCCA_L001_R2.fastq.merged.gz' >> libids.txt

echo '#!/bin/bash -e' > 1b.renamefiles.sh
echo '#SBATCH -p tgac-short # partition (queue)' >> 1b.renamefiles.sh
echo '#SBATCH -N 1 # number of nodes' >> 1b.renamefiles.sh
echo '#SBATCH -n 1 # number of tasks' >> 1b.renamefiles.sh
echo '#SBATCH --mem 8000' >> 1b.renamefiles.sh
echo '#SBATCH -t 0-00:20' >> 1b.renamefiles.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 1b.renamefiles.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 1b.renamefiles.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 1b.renamefiles.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 1b.renamefiles.sh
printf '\n' >> 1b.renamefiles.sh
echo "sed 's/^/mv /g'" $scripts"/"$libids ">" $libids2 >> 1b.renamefiles.sh
echo "sed -i '1 i\\n'" $libids2 >> 1b.renamefiles.sh
echo "sed -i '1 i\#!/bin/sh'" $libids2 >> 1b.renamefiles.sh
printf '\n' >> 1b.renamefiles.sh
echo "sh" $libids2 >> 1b.renamefiles.sh

# assign '1b.renamefiles.sh' to variable: JOBID2 and run
JOBID2=$( sbatch --dependency=afterok:${JOBID1} 1b.renamefiles.sh | awk '{print $4}' ) # JOB2 depends on JOB1 completing successfully

echo '# -- 1a. Adaptor trimming completed -- #'

echo '# -- 1b. Renaming files started -- #'

################################################################################################################

### 2. Read alignment

mkdir -p $readalign
cd $readalign

# 2a. Build genome indexes for bowtie2 alignment - variables assigned on command line input

echo '#!/bin/bash -e' > 2a.genomeindex.sh
echo '#SBATCH -p tgac-medium # partition (queue)' >> 2a.genomeindex.sh
echo '#SBATCH -N 1 # number of nodes' >> 2a.genomeindex.sh
echo '#SBATCH -n 1 # number of tasks' >> 2a.genomeindex.sh
echo '#SBATCH --mem 48000' >> 2a.genomeindex.sh
echo '#SBATCH -t 0-23:59' >> 2a.genomeindex.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 2a.genomeindex.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 2a.genomeindex.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 2a.genomeindex.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 2a.genomeindex.sh
printf '\n' >> 2a.genomeindex.sh
echo 'ml bowtie2/2.2.6' >> 2a.genomeindex.sh
printf '\n' >> 2a.genomeindex.sh
echo '#build indexes and assign to variables' >> 2a.genomeindex.sh
echo 'bowtie2-build $gFA $spG' >> 2a.genomeindex.sh
echo '#bowtie2-build /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa M_zebra_UMD1' >> 2a.genomeindex.sh
echo '#bowtie2-build /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta P_nyererei_v1' >> 2a.genomeindex.sh
echo '#bowtie2-build /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta A_burtoni_v1' >> 2a.genomeindex.sh
echo '#bowtie2-build /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta N_brichardi_v1' >> 2a.genomeindex.sh
echo '#bowtie2-build /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta O_niloticus_UMD1' >> 2a.genomeindex.sh

JOBID3=$( sbatch --dependency=afterok:${JOBID2} 2a.genomeindex.sh | awk '{print $4}' ) # JOB3 depends on JOB2 completing successfully

echo '# -- 1b. Renaming files completed -- #'

echo '# -- 2a.'$spID' genome index building started -- #'


# 2b. bowtie2 read alignments to genome indexes

echo '#!/bin/bash -e' > 2b.readalign.sh
echo '#SBATCH -p ei-largemem # partition (queue)' >> 2b.readalign.sh
echo '#SBATCH -N 1 # number of nodes' >> 2b.readalign.sh
echo '#SBATCH -c '$bwt_thread '# number of cores' >> 2b.readalign.sh
echo '#SBATCH --mem 512GB' >> 2b.readalign.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 2b.readalign.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 2b.readalign.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 2b.readalign.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 2b.readalign.sh
printf '\n' >> 2b.readalign.sh
echo 'ml bowtie2/2.2.6' >> 2b.readalign.sh
echo 'ml samtools/1.7' >> 2b.readalign.sh
printf '\n' >> 2b.readalign.sh
echo "awk -F' ' '{print \$2}' "$scripts"/"$libids " > " $reads >> 2b.readalign.sh
echo 'mapfile -t reads < '$reads' # ${reads[0]} calls read1 AND ${reads[1]} calls read2' >> 2b.readalign.sh
echo "awk -F' ' '{print \$2}' " $scripts"/"$libids " | awk -F'_' '{print \$1}' > "$prefix "# create a prefix file to iterate" >> 2b.readalign.sh
echo 'mapfile -t prefixmap < '$prefix '# assign prefixes to $prefixmap' >> 2b.readalign.sh
echo '# run bowtie2 with multimapping and threading, then output sorted BAM file' >> 2b.readalign.sh
echo 'srun bowtie2 -k' $multimapping '-X2000 --mm --threads' $bwt_thread '-x' $idx '-1 ${reads[0]} -2 ${reads[1]} 2>'$prefixmap$log '| samtools view -Su /dev/stdin | samtools sort -o $prefixmap'$bam >> 2b.readalign.sh
printf '\n' >> 2b.readalign.sh
echo 'samtools flagstat' $prefixmap$bam '>' $prefixmap$fgQC1 '# output alignment stats' >> 2b.readalign.sh

# ## this is if you are running several libraries (>2) in an array
# echo '#!/bin/bash -e' > 2b.readalign.sh
# echo '#SBATCH -p ei-largemem # partition (queue)' >> 2b.readalign.sh
# echo '#SBATCH -N 1 # number of nodes' >> 2b.readalign.sh
# echo '#SBATCH --cpus-per-task='$bwt_thread '# number of cores' >> 2b.readalign.sh
# echo '#SBATCH --array='$RAarray >> 2b.readalign.sh
# echo '#SBATCH --mem 512GB' >> 2b.readalign.sh
# echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 2b.readalign.sh
# echo "#SBATCH --mail-user=$email # send-to address" >> 2b.readalign.sh
# echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 2b.readalign.sh
# echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 2b.readalign.sh
# printf '\n' >> 2b.readalign.sh
# echo 'ml bowtie2/2.2.6' >> 2b.readalign.sh
# echo 'ml samtools/1.7' >> 2b.readalign.sh
# printf '\n' >> 2b.readalign.sh
# echo "grep '_L001_R1.fastq.merged.gz' " $scripts"/"$libids " | awk -F' '{print \$2}' | grep ^"$spID" > "$read1 "# create a list file of read1 for the array" >> 2b.readalign.sh
# echo "grep '_L001_R2.fastq.merged.gz' " $scripts"/"$libids " | awk -F' '{print \$2}' | grep ^"$spID" > "$read2 "# create a list file of read2 for the array" >> 2b.readalign.sh
# echo "awk -F' ' '{print \$2}' " $scripts"/"$libids " | awk -F'_' '{print \$1}' | grep "$spID" > "$prefix "# create a prefix file to iterate on the array" >> 2b.readalign.sh
# printf '\n' >> 2b.readalign.sh
# echo 'mapfile -t read1map < '$read1 '# assign read1 as elements to $read1map' >> 2b.readalign.sh
# echo 'mapfile -t read2map < '$read2 '# assign read2 as elements to $read2map' >> 2b.readalign.sh
# echo 'mapfile -t prefixmap < '$prefix '# assign prefixes to $prefixmap' >> 2b.readalign.sh
# printf '\n' >> 2b.readalign.sh
# echo '# run bowtie2 with multimapping and threading, then output sorted BAM file' >> 2b.readalign.sh
# echo 'srun bowtie2 -k' $multimapping '-X2000 --mm --threads' $bwt_thread '-x' $idx '-1 ${read1map[${SLURM_ARRAY_TASK_ID}]} -2 ${read2map[${SLURM_ARRAY_TASK_ID}]} 2>${prefixmap[${SLURM_ARRAY_TASK_ID}]}'$log '| samtools view -Su /dev/stdin | samtools sort -o ${prefixmap[${SLURM_ARRAY_TASK_ID}]}'$bam >> 2b.readalign.sh
# printf '\n' >> 2b.readalign.sh
# echo 'samtools flagstat ${prefixmap[${SLURM_ARRAY_TASK_ID}]}'$bam '>' $readalign'/${prefixmap[${SLURM_ARRAY_TASK_ID}]}'$fgQC1 '# output alignment stats' >> 2b.readalign.sh

# JOBID4=$( sbatch --dependency=afterok:${JOBID3} --array=$RAarray 2b.readalign.sh | awk '{print $4}' ) # JOB4 depends on JOB3 completing successfully

JOBID4=$( sbatch --dependency=afterok:${JOBID3} 2b.readalign.sh | awk '{print $4}' ) # JOB4 depends on JOB3 completing successfully

echo '# -- 2a.'$spID' genome index building completed -- #'

echo '# -- 2b.'$spID' read alignment started -- #'

################################################################################################################

### 3. Remove mitochondrial mapped and plot fragment distribution

mkdir -p $mtfilt
cd $mtfilt

## 3a. download mitochondrial genome as FASTA from NCBI
echo '#!/bin/bash -e' > 3.mtfilt_fragcount_A.sh
echo '#SBATCH -p ei-short # partition (queue)' >> 3.mtfilt_fragcount_A.sh
echo '#SBATCH -N 1 # number of nodes' >> 3.mtfilt_fragcount_A.sh
echo '#SBATCH -c 1 # number of cores' >> 3.mtfilt_fragcount_A.sh
echo '#SBATCH --mem 8000 # memory pool for all cores' >> 3.mtfilt_fragcount_A.sh
echo '#SBATCH -t 0-0:45 # time (D-HH:MM)' >> 3.mtfilt_fragcount_A.sh
echo '#SBATCH -o slurm.%j.out # STDOUT' >> 3.mtfilt_fragcount_A.sh
echo '#SBATCH -e slurm.%j.err # STDERR' >> 3.mtfilt_fragcount_A.sh
echo '#SBATCH --mail-type=END,FAIL,TIME_LIMIT_75 # notifications for job done & fail' >> 3.mtfilt_fragcount_A.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to addressUSERNAME=mehtat' >> 3.mtfilt_fragcount_A.sh
printf '\n' >> 3.mtfilt_fragcount_A.sh
echo "# this script will access the software node to download the FASTA file from the internet" >> 3.mtfilt_fragcount_A.sh
echo "USERNAME=$Usr" >> 3.mtfilt_fragcount_A.sh
echo 'HOSTNAME="software"' >> 3.mtfilt_fragcount_A.sh
echo "PWD=$(pwd)" >> 3.mtfilt_fragcount_A.sh
echo 'SCRIPT="cd '${PWD}'; wget -O '$mtID'.fasta https://www.ncbi.nlm.nih.gov/search/api/sequence/'$mtID'/?report=fasta; exit"' >> 3.mtfilt_fragcount_A.sh
echo 'ssh  -oStrictHostKeyChecking=no -l ${USERNAME} ${HOSTNAME} "${SCRIPT}"' >> 3.mtfilt_fragcount_A.sh

JOBID5=$( sbatch --dependency=afterok:${JOBID4} 3.mtfilt_fragcount_A.sh | awk '{print $4}' ) # JOB5 depends on JOB4 completing successfully

echo '# -- 2b.'$spID' read alignment completed -- #'

echo '# -- 3a.'$spID' mitochondrial genome download started: '$mtID' -- #'

## 3b. remove mitochondrial mapped reads from ATAC data and plot fragment lengths
echo '#!/bin/bash -e' > 3.mtfilt_fragcount_B.sh
echo '#SBATCH -p tgac-medium # partition (queue)' >> 3.mtfilt_fragcount_B.sh
echo '#SBATCH -N 1 # number of nodes' >> 3.mtfilt_fragcount_B.sh
echo '#SBATCH -n 1 # number of tasks' >> 3.mtfilt_fragcount_B.sh
echo '#SBATCH --mem 24000' >> 3.mtfilt_fragcount_B.sh
echo '#SBATCH -t 0-23:59' >> 3.mtfilt_fragcount_B.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 3.mtfilt_fragcount_B.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 3.mtfilt_fragcount_B.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 3.mtfilt_fragcount_B.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 3.mtfilt_fragcount_B.sh
printf '\n' >> 3.mtfilt_fragcount_B.sh
echo 'ml samtools/1.7' >> 3.mtfilt_fragcount_B.sh
printf '\n' >> 3.mtfilt_fragcount_B.sh
echo '# 2. create a blast database of each genome and store under variables' >> 3.mtfilt_fragcount_B.sh
echo 'ml blast/2.3.0' >> 3.mtfilt_fragcount_B.sh
echo "makeblastdb -in $gFA -parse_seqids -dbtype nucl # create a blast database of the genome assembly" >> 3.mtfilt_fragcount_B.sh
printf '\n' >> 3.mtfilt_fragcount_B.sh
echo '# 3. blast the mitochondrial genomes against the assembly - output tabular format' >> 3.mtfilt_fragcount_B.sh
echo "blastn -db $gFA -outfmt 6 -evalue 1e-3 -word_size 11 -show_gis -num_alignments 10 -max_hsps 20 -num_threads 5 -out $spID.genome_mt.blast -query $mtgen # blast the mitochondrial genome against the input genome assembly, output tabular format" >> 3.mtfilt_fragcount_B.sh
printf '\n' >> 3.mtfilt_fragcount_B.sh
echo '# 4. Python script: filter BLAST hits based on pident>=93 and evalue<=1e-10; then pident>75% according to BLAST hit and alignment length to mitochondrial genome' >> 3.mtfilt_fragcount_B.sh
echo 'ml python/3.5' >> 3.mtfilt_fragcount_B.sh
echo "python $scripts/ATAC_Bioinf_pipeline_v2b_part3a.py $blstout $mtgen # output is stored as variable "'$mtscaff at top' >> 3.mtfilt_fragcount_B.sh
printf '\n' >> 3.mtfilt_fragcount_B.sh
echo '# 5. Store the genome scaffolds names that are mitochondrial genome under variable $scaffarray and create grep command $grepscaff' >> 3.mtfilt_fragcount_B.sh
echo "IFS=$'\\\n' scaffarray=("'$(cut -f2 '$mtscaff" | awk ""'"'!x[$0]++'"'))" '# this takes the $mtscaff as input, takes unique genome scaffolds (col2) that match mtDNA (using awk instead of sort -u so top hit ordering retained) and then assigns to the variable $scaffarray. Accessed using ${scaffarray[0]}..${scaffarray[2]}' >> 3.mtfilt_fragcount_B.sh
echo 'echo Scaffold/s matching mitochondrial genome are/is: ${scaffarray[@]} # this will echo each scaffold that matches mtDNA' >> 3.mtfilt_fragcount_B.sh
echo 'echo Total number of scaffolds matching mitochondrial genome is: ${#scaffarray[@]} # this will echo the number of scaffolds that match mtDNA' >> 3.mtfilt_fragcount_B.sh
echo '# for sc in ${scaffarray[@]}; do echo $sc ; done # this will loop through each element (Scaff) of the array and echo on new line' >> 3.mtfilt_fragcount_B.sh
echo 'grepscaff=$(echo ${scaffarray[@]} | sed '"'s/ / | grep -v /g' | sed 's/^/grep -v /') # Assign grep commands to variable: input grep -v and then pipe in front of each element of the array e.g. UNK2407 UNK3173 UNK1343 > grep -v UNK2407 | grep -v UNK3173 | grep -v UNK1343" >> 3.mtfilt_fragcount_B.sh
printf '\n' >> 3.mtfilt_fragcount_B.sh >> 3.mtfilt_fragcount_B.sh
echo '# 6. Loop through variable of scaffold/s matching the mitochondrial genome and filter reads assigned to these scaffolds.' >> 3.mtfilt_fragcount_B.sh
echo 'ml samtools/1.3' >> 3.mtfilt_fragcount_B.sh
echo 'ml perl_activeperl/5.18' >> 3.mtfilt_fragcount_B.sh
echo 'ml zlib/1.2.8' >> 3.mtfilt_fragcount_B.sh
printf '\n' >> 3.mtfilt_fragcount_B.sh
echo '# 1. assign mtDNA filtered bam to new file with extension .nochrM.bam' >> 3.mtfilt_fragcount_B.sh
echo '# 2. remove all reads mapping to mitochondrial scaffold/s and index' >> 3.mtfilt_fragcount_B.sh
echo "for bam_file in $readalign/*.bam; do bam_file_nochrM="'$(echo $bam_file | sed -e '"'s/.bam/.nochrM.bam/'); samtools idxstats "'$bam_file | cut -f1 | $grepscaff | xargs samtools view -b $bam_file > $bam_file_nochrM; samtools index $bam_file_nochrM; done #samtools view -h $bam_file | grep ${scaffarray[0]} | wc -l # you can type this to test' >> 3.mtfilt_fragcount_B.sh
printf '\n' >> 3.mtfilt_fragcount_B.sh
echo '# 7. Calculate fragment length count for each' >> 3.mtfilt_fragcount_B.sh
echo 'frag_length=$(echo $bam_file_nochrM | sed -e '"'s/.bam/_frag_length_count.txt/')" >> 3.mtfilt_fragcount_B.sh
echo 'samtools view $bam_file_nochrM | awk '"'"'$9>0'"' | cut -f9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > "'$frag_length' >> 3.mtfilt_fragcount_B.sh
printf '\n' >> 3.mtfilt_fragcount_B.sh
echo '# 8. Plot fragment length count in R' >> 3.mtfilt_fragcount_B.sh
echo 'source R-3.5.2' >> 3.mtfilt_fragcount_B.sh
echo "R CMD BATCH --no-save --no-restore --args "'$bam_file_nochrM '"$scripts/ATAC_Bioinf_pipeline_v2b_part3b.R ATAC_Bioinf_pipeline_v2b_part3b.Rout # this creates two files - Rplots.pdf (which has the image!) and another (empty) image file with the actual filename. Simply rename Rplots.pdf" >> 3.mtfilt_fragcount_B.sh
echo 'mv Rplots.pdf "$(basename "$bam_file_nochrM" .bam).fraglength.pdf" # rename Rplots.pdf to *.fraglength.pdf' >> 3.mtfilt_fragcount_B.sh

JOBID6=$( sbatch --dependency=afterok:${JOBID5} 3.mtfilt_fragcount_B.sh | awk '{print $4}' ) # JOB6 depends on JOB5 completing successfully

echo '# -- 3a.'$spID' mitochondrial genome downloaded: '$mtID' -- #'

echo '# -- 3b.'$spID' mitochondrial removed and fragment dist plot started -- #'

################################################################################################################

### 4. Post alignment filtering
#   4a. Filter reads (Sort, Map and remove duplicates, Remove reads unmapped, not primary alignment, reads failing platform, duplicates) - samtools > final bam

mkdir -p $filtdir
cd $filtdir

echo '#!/bin/bash -e' > 4.postalign_filt.sh
echo '#SBATCH -p tgac-medium # partition (queue)' >> 4.postalign_filt.sh
echo '#SBATCH -N 1 # number of nodes' >> 4.postalign_filt.sh
echo '#SBATCH -C 1 # number of cores' >> 4.postalign_filt.sh
echo '#SBATCH -T 2 # number of threads per core' >> 4.postalign_filt.sh
echo '#SBATCH --mem 24000' >> 4.postalign_filt.sh
echo '#SBATCH -t 0-05:59' >> 4.postalign_filt.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 4.postalign_filt.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 4.postalign_filt.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 4.postalign_filt.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 4.postalign_filt.sh
printf '\n' >> 4.postalign_filt.sh
echo 'ml samtools/1.3' >> 4.postalign_filt.sh
echo 'ml sambamba/0.6.5' >> 4.postalign_filt.sh
echo 'ml python/3.5' >> 4.postalign_filt.sh
echo 'ml zlib/1.2.8' >> 4.postalign_filt.sh
echo 'ml glib/2.40' >> 4.postalign_filt.sh
echo 'ml java/8.45' >> 4.postalign_filt.sh
echo 'ml picard/1.140' >> 4.postalign_filt.sh
echo 'ml bedtools/2.26.0' >> 4.postalign_filt.sh
echo 'ml zlib/1.2.8' >> 4.postalign_filt.sh
echo 'ml glib/2.40' >> 4.postalign_filt.sh
echo 'ml R/3.2.3' >> 4.postalign_filt.sh
printf '\n' >> 4.postalign_filt.sh
echo "for bam_file in ${mtfilt}/*.nochrM.bam; do" >> 4.postalign_filt.sh
echo '\t# variables for output files' >> 4.postalign_filt.sh
echo '\tbam_file_sorted=$(echo $bam_file | sed -e '"'s/.bam/.sorted.bam/')" >> 4.postalign_filt.sh
echo '\tbam_file_dup=$(echo $bam_file | sed -e '"'s/.bam/.sorted.dup.bam/')" >> 4.postalign_filt.sh
echo '\tnodup_filt_bam_file=$(echo $bam_file | sed -e '"'s/.bam/.nodup.filt.bam/') # final bam file" >> 4.postalign_filt.sh
echo '\tnodup_filt_bam_index_file=$(echo $bam_file | sed -e '"'s/.bam/.nodup.filt.bam.bai/') # index file" >> 4.postalign_filt.sh
echo '\tnodup_filt_bam_file_mapstats=$(echo $bam_file | sed -e '"'s/.bam/.flagstat.qc/') # QC file" >> 4.postalign_filt.sh
echo '\tpbc_file_qc=$(echo $bam_file | sed -e '"'s/.bam/.pbc.qc/') # library complexity" >> 4.postalign_filt.sh
echo '\tnodup_filt_bam_file_sorted=$(echo $bam_file | sed -e '"'s/.bam/.srt.bam/') # final bam file, sorted (temp)" >> 4.postalign_filt.sh
echo '\t# Filter reads' >> 4.postalign_filt.sh
echo "\tsambamba sort -m 24G -t 2 -o $filtdir/"'$bam_file_sorted -u $bam_file' >> 4.postalign_filt.sh
echo "\tsambamba markdup -l 0 -t 10 $filtdir/"'$bam_file_sorted '"$filtdir/"'$bam_file_dup' >> 4.postalign_filt.sh
echo "\tsamtools view -F 1804 -f 2 -q 30 -b $filtdir/"'$bam_file_dup > '"$filtdir/"'$nodup_filt_bam_file' >> 4.postalign_filt.sh
echo "\tsamtools index $filtdir/"'$nodup_filt_bam_file '"$filtdir/"'$nodup_filt_bam_index_file' >> 4.postalign_filt.sh
echo "\tsamtools flagstat $filtdir/"'$nodup_filt_bam_file > '"$filtdir/"'$nodup_filt_bam_file_mapstats' >> 4.postalign_filt.sh
echo '\t# Plotting the fragment length distribution' >> 4.postalign_filt.sh
echo "\tjava -Xmx2g -jar /tgac/software/production/picardtools/1.84/x86_64/bin/CollectInsertSizeMetrics.jar \R=$gFA \I=$filtdir/"'$nodup_filt_bam_file \O="'"$filtdir/"'${nodup_filt_bam_file}_PicardInsertMetrics.jar.txt" \H="'"$filtdir/"'${nodup_filt_bam_file}_insert_size_histogram.pdf" \M=0.5' >> 4.postalign_filt.sh
echo 'done' >> 4.postalign_filt.sh

# for bam_file in ${mtfilt}/*.nochrM.bam; do
#   # variables for output files
#   bam_file_sorted=$(echo $bam_file | sed -e 's/.bam/.sorted.bam/')
#   bam_file_dup=$(echo $bam_file | sed -e 's/.bam/.sorted.dup.bam/')
#   nodup_filt_bam_file=$(echo $bam_file | sed -e 's/.bam/.nodup.filt.bam/') # final bam file
#   nodup_filt_bam_index_file=$(echo $bam_file | sed -e 's/.bam/.nodup.filt.bam.bai/') # index file
#   nodup_filt_bam_file_mapstats=$(echo $bam_file | sed -e 's/.bam/.flagstat.qc/') # QC file
#   pbc_file_qc=$(echo $bam_file | sed -e 's/.bam/.pbc.qc/') # library complexity
#   nodup_filt_bam_file_sorted=$(echo $bam_file | sed -e 's/.bam/.srt.bam/') # final bam file, sorted (temp)
#   # Filter reads
#   sambamba sort -m 24G -t 2 -o $filtdir/$bam_file_sorted -u $bam_file
#   sambamba markdup -l 0 -t 10 $filtdir/$bam_file_sorted $filtdir/$bam_file_dup
#   samtools view -F 1804 -f 2 -q 30 -b $filtdir/$bam_file_dup > $filtdir/$nodup_filt_bam_file
#   samtools index $filtdir/$nodup_filt_bam_file $filtdir/$nodup_filt_bam_index_file
#   samtools flagstat $filtdir/$nodup_filt_bam_file > $filtdir/$nodup_filt_bam_file_mapstats
#   # Plotting the fragment length distribution
#   java -Xmx2g -jar /tgac/software/production/picardtools/1.84/x86_64/bin/CollectInsertSizeMetrics.jar \R=$gFA \I=$filtdir/$nodup_filt_bam_file \O="$filtdir/${nodup_filt_bam_file}_PicardInsertMetrics.jar.txt" \H="$filtdir/${nodup_filt_bam_file}_insert_size_histogram.pdf" \M=0.5
# done

JOBID7=$( sbatch --dependency=afterok:${JOBID6} 4.postalign_filt.sh | awk '{print $4}' ) # JOB7 depends on JOB6 completing successfully

echo '# -- 3b.'$spID' mitochondrial removed and fragment dist plot completed -- #'

echo '# -- 4.'$spID' Post alignment filtering started -- #'

################################################################################################################

### 5. ATAC peak calling (test vs control)
# 	5a. Convert PE Bam to tagalign (start/end positions of each read) - bedtools
#   5b. TSS enrichment - plot
  # The TSS enrichment calculation is a signal to noise calculation.
  # The reads around a reference set of TSSs are collected to form an aggregate distribution of reads centered on the TSSs and extending to 1000 bp in either direction (for a total of 2000bp).
  # This distribution is then normalized by taking the average read depth in the 100 bps at each of the end flanks of the distribution (for a total of 200bp of averaged data) and calculating a fold change at each position over that average read depth.
  # This means that the flanks should start at 1, and if there is high read signal at transcription start sites (highly open regions of the genome) there should be an increase in signal up to a peak in the middle.
  # We take the signal value at the center of the distribution after this normalization as our TSS enrichment metric.
# 	5c. TN5 shifting of tagaligns - shift reads +4 bp for the +strand and -5 bp for the -strand
# 	5d. peak calling - macs2 NOTE: consider running another peak-calling program and take the intersection. Also, for analysis only consider open-chromatin so filter based on that?


## ~ COMPLETED THIS CODE AS OF 08/04/2020 - NEEDS TESTING USING DATA AND THEN TURN INTO ECHO ~ ##

## Think about how the gDNA control is called:
  # 1. There needs to be a separate pipeline for gDNA controls that runs until peak calling e.g no mitochondrial reads removed - this script can be adapted at the end

mkdir -p $peakcall
cd $peakcall

echo '#!/bin/bash -e' > 5.peakcall.sh
echo '#SBATCH -p tgac-medium # partition (queue)' >> 5.peakcall.sh
echo '#SBATCH -N 1 # number of nodes' >> 5.peakcall.sh
echo '#SBATCH -C 1 # number of cores' >> 5.peakcall.sh
echo '#SBATCH -T 1 # number of threads per core' >> 5.peakcall.sh
echo '#SBATCH --mem 24000' >> 5.peakcall.sh
echo '#SBATCH -t 0-05:59' >> 5.peakcall.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 5.peakcall.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 5.peakcall.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 5.peakcall.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo 'ml MACS' >> 5.peakcall.sh
echo 'ml python/3.5' >> 5.peakcall.sh
echo 'ml bedtools/2.25.0' >> 5.peakcall.sh
echo 'ml GCC' >> 5.peakcall.sh
echo 'ml zlib' >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo 'ulimit -Sn 10000 # set the open file limit to 10,000' >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh

## ONCE CHECKED, MOVE VARIABLES TO THE TOP
Test1=$filtdir/$spID'.nodup.filt.bam' # the filename needs checking!!
Control1=$(echo $Test1 | sed -e 's/ATAC/gDNA/g') # the filename needs checking - does it pick up the corresponding gDNA dir AND file!

# Sum total length of all scaffolds in assembly - needed as input for Macs2
source bioawk-1.0
Gensz=$(bioawk -c fastx '{ print $name, length($seq) }' < $gFA | awk -F"\t" '{print;x+=$2}END{print "Total " x}' | tail -1 | sed 's/Total //g') # this will output the sum length of scaffolds and assign to variable $Gensz

# MACS2 parameters
Macs2PvalThresh="0.05"  # The p-value threshold for calling peaks
Macs2SmoothWindow=150  # The window size to smooth alignment signal over
Macs2ShiftSize=$(python -c "print(int(${Macs2SmoothWindow}/2))") # This uses --extsize 75; DNA wrapped around the nucleosome is ~147bp, however in ATAC-seq, we cut (using transposase) at 'any' open chromatin region so average fragment length can be < 150
output_prefix="$peakcall/$spID"

# 5a. Convert each bam to a .tagAlign - contains the start/end positions of each read:
tagalign_test1=$peakcall/$(echo $(basename $Test1) | sed -e 's/.bam/.tagAlign.gz/')
tagalign_control1=$peakcall/$(echo $(basename $Control1) | sed -e 's/.bam/.tagAlign.gz/')
bedtools bamtobed -i $Test1 | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -c  > $tagalign_test1
bedtools bamtobed -i $Control1 | awk 'BEGIN{OFS="\t"}{$4="N";$5="1000";print $0}' | gzip -c > $tagalign_control1

# 5b. TSS enrichment plotting
# This calls one python script 'ATAC_Bioinf_pipeline_v2b_part5bD.py' and requires encode_lib_genomic.py and encode_lib_common.py are in the same scripts folder

# create a 2kb window around TSS (+/- 1kb) bed file e.g.
# chr1	134210701	134214701	+
# chr1	33724603	33728603	-

# 5bA. Protein-coding genes GTF > BED: variable $annot of gzipped GTF file (*gtf.gz); or 2) longest protein-coding gene annotations as 6-column BED: col1-scaff,col2-start,col3-end,col4-geneID,col5-XX,col6-strand"
# output is genebed=($peakcall/$spG'_refGene.bed') defined as variable at top

source bedops-2.4.28

case "$annot" in
*gtf.gz)
        zcat $annot | awk 'OFS="\t" {if ($3=="gene") {print $1,$4-1,$5,$10,$7,$18}}' | tr -d '";' | grep -wiF 'protein_coding' | awk '{print $1,$2,$3,$4,$5}' OFS='\t' > $genebed # GTF > 0-based BED of protein_coding
        ;;
*.bed)
        awk '{print $1,$2,$3,$4,$6}' OFS='\t' $annot > $genebed # BED format ONLY for my files that have six cols (where 6th col is strand)
        ;;
esac

# 5bB. Then split them by strand and pad around the stranded-start position of the annotation (taking TSS +/- 1000=1kb)
awk '($5 == "+") { print $0 }' $genebed | awk 'BEGIN{ OFS="\t" }($2 > 1000){ print $1, ($2 - 1000), ($2 + 1000), $4, $5  }' > $genebed'.tss.for.padded.bed'
awk '($5 == "-") { print $0 }' $genebed | awk 'BEGIN{ OFS="\t" }($3 > 1000){ print $1, ($3 - 1000), ($3 + 1000), $4, $5  }' > $genebed'.tss.rev.padded.bed'
bedops --everything $genebed'.tss.for.padded.bed' $genebed'.tss.rev.padded.bed' > $genebedtss

# 5bC. Keep only TSS regions within chromosomal bounds - prep scaffold sizes file (col1=scaffoldID; col2=0; col3=length) from genome fasta
bioawk -c fastx '{ print $name, length($seq) }' < $gFA | awk '{print $1,"0",$2}' OFS='\t' > $scafflen
bedops --element-of 100% $genebedtss $scafflen > $genebedtss2

# 5bD. Use final TSS (+/- 1kb) bed file as input to calculate TSS enrichment and plot with python script 'ATAC_Bioinf_pipeline_v2b_part5bD.py'
ml python/3.5
python3 ATAC_Bioinf_pipeline_v2b_part5bD-a.py $fastqr1 # input fastq can be native or gzipped
python3 ATAC_Bioinf_pipeline_v2b_part5bD.py $Test1 $genebedtss2 $spID $read_len $scafflen # usage: python3 ATAC_Bioinf_pipeline_v2b_part5bD.py 'FINAL_BAM' 'TSS' 'OUTPUT_PREFIX' 'read_len' 'CHROMSIZES'

# 5c. Tn5 shifting of tagaligns
shifted_tag="$peakcall/$spID"'.tn5.tagAlign.gz'
zcat $tagalign_test1 | awk -F $'\t' 'BEGIN {OFS = FS}{ if ($6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}' | gzip -c > $shifted_tag

# 5d. Run MACS2:
macs2 callpeak -t $shifted_tag -c $tagalign_control1 -f BED -n $output_prefix -g $GenSz -p $Macs2PvalThresh --nomodel --shift -$Macs2ShiftSize --extsize $Macs2SmoothWindow -B --SPMR --keep-dup all --call-summits
# --nomodel and --extsize 150 tells MACS2 to use 150bp as fragment size to pileup sequencing reads.
# -g XX lets MACS2 consider a genome size as background.
# -B --SPMR ask MACS2 to generate pileup signal file of 'fragment pileup per million reads' in bedGraph format.

# Generate a fold change file comparing the sample to the control and logLR track
macs2 bdgcmp -t $output_prefix\_treat_pileup.bdg -c $output_prefix\_control_lambda.bdg -o $output_prefix\_FE.bdg -m FE
macs2 bdgcmp -t $output_prefix\_treat_pileup.bdg -c $output_prefix\_control_lambda.bdg -o $output_prefix\_logLR.bdg -m logLR -p 0.00001


JOBID8=$( sbatch --dependency=afterok:${JOBID7} 5.peakcall.sh | awk '{print $4}' ) # JOB8 depends on JOB7 completing successfully

echo '# -- 4.'$spID' Post alignment filtering completed -- #'

echo '# -- 5.'$spID' Peak calling started -- #'


################################################################################################################

### 6. IDR on all pairs of replicates, self-pseudoreplicates and pooled pseudoreplicates - IDR is optional. The IDR peaks are a subset of the naive overlap peaks that pass a specific IDR threshold of 10%.
# 	6a. IDR of true replicates
# 	6b. Compute Fraction of Reads in Peaks (FRiP)

## ~ INSERT CODE HERE ~ ##

JOBID9=$( sbatch --dependency=afterok:${JOBID8} XX.sh | awk '{print $4}' ) # JOB9 depends on JOB8 completing successfully

echo '# -- 5.'$spID' Peak calling completed -- #'

echo '# -- 6.'$spID' IDR started -- #'

################################################################################################################

### 7. Create signal tracks - bedtools

## ~ INSERT CODE HERE ~ ##

JOBID10=$( sbatch --dependency=afterok:${JOBID9} XX.sh | awk '{print $4}' ) # JOB10 depends on JOB9 completing successfully

echo '# -- 6.'$spID' IDR completed -- #'

echo '# -- 7.'$spID' Creating signal tracks has started -- #'

################################################################################################################

### 8. Annotation:
# 	8a. TSS enrichment
# 	8b. Fraction of Reads in annotated regions

## ~ INSERT CODE HERE ~ ##

JOBID11=$( sbatch --dependency=afterok:${JOBID10} XX.sh | awk '{print $4}' ) # JOB11 depends on JOB10 completing successfully

echo '# -- 7.'$spID' Creating signal tracks has completed -- #'

echo '# -- 8.'$spID' Peak annotation has started -- #'

################################################################################################################

### NOTE: be careful when considering differential peaks as some may be only offset by a few bases. In this secnario, consider the average number of mapped reads over a window.

## ~ INSERT CODE HERE ~ ##

JOBID12=$( sbatch --dependency=afterok:${JOBID11} XX.sh | awk '{print $4}' ) # JOB12 depends on JOB11 completing successfully

echo '# -- 8.'$spID' Peak annotation has completed -- #'
