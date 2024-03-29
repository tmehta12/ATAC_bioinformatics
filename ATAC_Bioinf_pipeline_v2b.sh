#!/bin/sh

#!/bin/bash -e
#SBATCH -p ei-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 8000 # memory pool for all cores
#SBATCH -t 4-23:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

################################################################################################################

# ATAC-seq pipeline - Part 2
# March 2020: Tarang K. Mehta, Earlham Institute, Norwich, UK

################################################################################################################

# Script usage: ./ATAC_Bioinf_pipeline_v2b.sh -s "spID" -g "spG" -f "gFA" -m "mtID" -u "Usr" -a "annot"
# e.g. ./ATAC_Bioinf_pipeline_v2b.sh -s Mz1_L_ATAC -g M_zebra_UMD1 -f /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -m KT221043 -u mehtat -a M_zebra_UMD1.gtf.gz
# Note: Script is adapted for SBATCH usage and uses gDNA controls (can be removed from peak calling step)

## Place this script and the following files in $WD (created in first script)
# 1. As used in './ATAC_Bioinf_pipeline_v2a.sh': a 2-column space-delimited table where col1='R1/R2 filename's col2='desired species renamed filename: species_tissue_experiment e.g. Mz_L_ATAC/gDNA'
# 2. Scripts:
  # ATAC_Bioinf_pipeline_v2b_part3a.py
  # ATAC_Bioinf_pipeline_v2b_part3b.R
# 3. Run as an sbatch script with 8Gb memory and ~6 days runtime - will spawn off other jobs

################################################################################################################

# ~ This pipeline follows on from 'ATAC_Bioinf_pipeline_v2a.sh' and 'ATAC_Bioinf_pipeline_v2b_gDNA.sh', that:

# 0. Merges files sequenced over multiple lanes
# 1-2: Trimming and alignment of gDNA reads

################################################################################################################

# ~ This pipeline is ran as species-specific, and contains the following components:

# 1. Trim adaptors - trimgalore
# 2. Read alignment - bowtie2 > samtools sorted bam
# 3. Remove mitochondrial mapped - samtools > new bam
# 	3a. Fragment/insert size distribution
# 4. Post alignment filtering
#   4a. Filter reads (Sort, Map and remove duplicates, Remove reads unmapped, not primary alignment, reads failing platform, duplicates) - samtools > final bam
# 5. ATAC peak-calling
# 	5a. Convert PE Bam to tagalign (start/end positions of each read) - bedtools
# 	5b. TN5 shifting of tagaligns - shift reads +4 bp for the +strand and -5 bp for the -strand
# 	5c. count-based peak calling using Poisson distribution - macs2
#   5d. count-based peak calling using a program that considers properly paired, unpaired and secondary alignments (unlike MACS2) - Genrich
#   5e. markov model based peak calling specific for ATAC-seq data - HMMRATAC (this is currently surpressed to run)
# 6. bed to bigbed conversion for narrowpeaks - bedClip and bedToBigbed from ucsc_tools

################################################################################################################

# Setting parameters for command line input

helpFunction()
{
   echo ""
   echo "Usage: $0 -s spID -g spG -f gFA -m mtID -u Usr -a annot"
   echo -e "\t-s spID = Species ID, preferably two short letters, individual and tissue e.g. Metriaclima zebra Individual 1 Liver ATAC/gDNA = Mz1_L_ATAC/Mz1_L_gDNA Note: this naming convention needs to be the same as renaming in space delimited file"
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

scripts=(/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/data/ATACseq/3.run2) # place all scripts in the topmost directory - create this separately
WD=(/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/data/ATACseq/3.run2/$spID) # insert the working directory
email=Tarang.Mehta@earlham.ac.uk # SBATCH out and err send to address

### 1. Trim adaptors and renaming
rawreaddir=($WD/0.rawreads) # assign raw reads dir
trimdir=($WD/1.adaptor_trimming) # assign trimmed reads dir
trimarray=0-0 # INSERT the number range of paired *fastq.merged.gz to trim in zero base e.g. 10 pairs = 0-9
libids=($scripts/libids.txt) # 2-col space-delimited file where col1 is the *_{R1,R2}.fastq.merged.gz and col2 is the species_tissueID_experiment_barcode_{R1,R2}.fastq.merged.gz
libids1=($trimdir/libids1.txt) # this has the new trimmed ID gzipped fastq file paths
libids2=($trimdir/libids.sh)

### 2. Read alignment
readalign=($WD/2.read_alignment)
RAarray=0-0 # INSERT the number range of paired Mz *fastq.merged.gz to align in zero base e.g. 10 pairs = 0-9
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
bam_file_nochrM=($spID.nochrM.bam)

### 4. Post alignment filtering
filtdir=($WD/4.postalign_filt) # post alignment filtering directory

### 5. ATAC peak calling (test vs control)
peakcall=($WD/5.peak_calling) # ATAC peak calling directory
genebed=($peakcall/$spG'_refGene.bed') # GTF (or abnormal BED) > BED file of protein coding genes only
genebedtss=($peakcall/$spG'_refGene.tss.padded.bed') # BED file of protein coding genes with TSS +/- 1kb
scafflen=($peakcall/$spG'_scaffbounds.bed') # scaffold length boundaries of input genome assembly (col1=scaff; col2=0; col3=scaffold_length)
genebedtss2=($peakcall/$spG'_refGene.tss.padded.filt.bed') # BED file of protein coding genes with TSS +/- 1kb and filtered for any out of bound genes
fastqr1=($rawreaddir/*_R1.fastq.merged.gz)
read_len=($rawreaddir/*_read_length.txt)
Test1=$filtdir/$spID'.nochrM.nodup.filt.bam'
Control1=$(echo $Test1 | sed -e 's/_ATAC/_gDNA/g' | sed 's/.nochrM.nodup.filt//' | sed 's/4.postalign_filt/2.read_alignment/')
# MACS2 parameters
Macs2PvalThresh="0.05"  # The p-value threshold for calling peaks
Macs2SmoothWindow=150  # The window size to smooth alignment signal over
Macs2ShiftSize=$(python -c "print(int(${Macs2SmoothWindow}/2))") # This uses --extsize 75; DNA wrapped around the nucleosome is ~147bp, however in ATAC-seq, we cut (using transposase) at 'any' open chromatin region so average fragment length can be < 150
output_prefix="$peakcall/$spID"
tagalign_test1=$peakcall/$(echo $(basename $Test1) | sed -e 's/.bam/.tagAlign.gz/')
tagalign_control1=$(echo $readalign | sed -e 's/_ATAC/_gDNA/g')/$(echo $(basename $Control1) | sed -e 's/.bam/.tagAlign.gz/')
shifted_tag="$peakcall/$spID"'.tn5.tagAlign.gz'
Test2=$spID'.nochrM.nodup.filt.querysorted.bam'
Control2=$(echo $spID'.bam' | sed -e 's/_ATAC/_gDNA.querysorted/g')
Test1index=$filtdir/$spID'.nochrM.nodup.filt.bam.bai'
scafflen2=($peakcall/$spG'_scaffbounds.txt')

### 6. bed to bigbed conversion
bigbed=($peakcall/${spID}_peaks.narrowPeak.bb)
peak=($peakcall/${spID}_peaks.narrowPeak)
peakgz=($peakcall/${spID}_peaks.narrowPeak.gz)


################################################################################################################

### 1. Trim adaptors
mkdir -p $trimdir
cd $trimdir

# run as an array in sbatch script

echo '#!/bin/bash -e' > 1a.trimadaptors.sh
echo '#SBATCH -p ei-medium # partition (queue)' >> 1a.trimadaptors.sh
echo '#SBATCH -N 1 # number of nodes' >> 1a.trimadaptors.sh
echo '#SBATCH -n 1 # number of tasks' >> 1a.trimadaptors.sh
echo "#SBATCH --array=$trimarray" >> 1a.trimadaptors.sh
echo '#SBATCH --mem-per-cpu 24000' >> 1a.trimadaptors.sh
echo '#SBATCH -t 1-23:59' >> 1a.trimadaptors.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 1a.trimadaptors.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 1a.trimadaptors.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 1a.trimadaptors.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 1a.trimadaptors.sh
printf '\n' >> 1a.trimadaptors.sh
echo 'source trim_galore-0.5.0' >> 1a.trimadaptors.sh
echo 'ml java' >> 1a.trimadaptors.sh
echo 'source fastqc-0.11.9' >> 1a.trimadaptors.sh
printf '\n' >> 1a.trimadaptors.sh
echo "for read1 in $rawreaddir/*_R1.fastq.merged.gz; do read2="'$(echo $read1 | sed'" 's/_R1.fastq.merged.gz/_R2.fastq.merged.gz/'); "'echo $read1 >> trimarrayread1; echo $read2 >> trimarrayread2; done' >> 1a.trimadaptors.sh
printf '\n' >> 1a.trimadaptors.sh
echo 'mapfile -t read1 < trimarrayread1 # assign read1 as elements to $read1 variable' >> 1a.trimadaptors.sh
echo 'mapfile -t read2 < trimarrayread2 # assign read2 as elements to $read2 variable' >> 1a.trimadaptors.sh
printf '\n' >> 1a.trimadaptors.sh
echo "srun trim_galore --output_dir $trimdir --paired"' --fastqc ${read1[${SLURM_ARRAY_TASK_ID}]} ${read2[${SLURM_ARRAY_TASK_ID}]}' >> 1a.trimadaptors.sh

echo '# -- 1a. Adaptor trimming started -- #'

# assign '1a.trimadaptors.sh' to variable: JOBID1 and run
JOBID1=$( sbatch -W --array=$trimarray 1a.trimadaptors.sh | awk '{print $4}' ) # Run the first job and then store the first job to variable JOBID1 (taken by awk once run)

# rename the files according to species, tissue and experiment - provide this as a 2-column SPACE-delimited table that will be assigned as a variable above, placed in the trimdir
# only provide for the two paired files you are working on and not all files

# Example is (note, these are made up examples)
# echo 'PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R1.fastq.merged.gz Pn1_T_ATAC_TCGCCTGC-AACCGCCA_L001_R1.fastq.merged.gz' > libids.txt
# echo 'PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R2.fastq.merged.gz Pn1_T_ATAC_TCGCCTGC-AACCGCCA_L001_R2.fastq.merged.gz' >> libids.txt

echo '#!/bin/bash -e' > 1b.renamefiles.sh
echo '#SBATCH -p ei-short # partition (queue)' >> 1b.renamefiles.sh
echo '#SBATCH -N 1 # number of nodes' >> 1b.renamefiles.sh
echo '#SBATCH -n 1 # number of tasks' >> 1b.renamefiles.sh
echo '#SBATCH --mem 8000' >> 1b.renamefiles.sh
echo '#SBATCH -t 0-00:20' >> 1b.renamefiles.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 1b.renamefiles.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 1b.renamefiles.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 1b.renamefiles.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 1b.renamefiles.sh
printf '\n' >> 1b.renamefiles.sh
echo "grep $spID $libids | sed 's/R1.fastq.merged.gz/R1.fastq.merged.gz_val_1.fq.gz/g' | sed 's/R2.fastq.merged.gz/R2.fastq.merged.gz_val_2.fq.gz/g' > $libids1 # grep the relevant species files from the long list" >> 1b.renamefiles.sh
echo "sed 's/^/mv /g' $libids1 > $libids2" >> 1b.renamefiles.sh
echo "sed -i '1 i\#!/bin/sh' $libids2" >> 1b.renamefiles.sh
printf '\n' >> 1b.renamefiles.sh
echo "sh $libids2" >> 1b.renamefiles.sh

echo '# -- 1a. Adaptor trimming completed -- #'

echo '# -- 1b. Renaming files started -- #'

# assign '1b.renamefiles.sh' to variable: JOBID2 and run
JOBID2=$( sbatch -W --dependency=afterok:${JOBID1} 1b.renamefiles.sh | awk '{print $4}' ) # JOB2 depends on JOB1 completing successfully


################################################################################################################

### 2. Read alignment

mkdir -p $readalign
cd $readalign

# 2a. Build genome indexes for bowtie2 alignment - variables assigned on command line input

echo '#!/bin/bash -e' > 2a.genomeindex.sh
echo '#SBATCH -p ei-medium # partition (queue)' >> 2a.genomeindex.sh
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
echo "bowtie2-build $gFA $spG" >> 2a.genomeindex.sh
echo '#bowtie2-build /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa M_zebra_UMD1' >> 2a.genomeindex.sh
echo '#bowtie2-build /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/PunNye1.0/P_nyererei_v1.assembly.fasta P_nyererei_v1' >> 2a.genomeindex.sh
echo '#bowtie2-build /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/HapBur1.0/H_burtoni_v1.assembly.fasta A_burtoni_v1' >> 2a.genomeindex.sh
echo '#bowtie2-build /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/NeoBri1.0/N_brichardi_v1.assembly.fasta N_brichardi_v1' >> 2a.genomeindex.sh
echo '#bowtie2-build /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/O_niloticus_UMD1/O_niloticus_UMD1.fasta O_niloticus_UMD1' >> 2a.genomeindex.sh

echo '# -- 1b. Renaming files completed -- #'

echo '# -- 2a.'$spID' genome index building started -- #'

JOBID3=$( sbatch -W --dependency=afterok:${JOBID2} 2a.genomeindex.sh | awk '{print $4}' ) # JOB3 depends on JOB2 completing successfully


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
echo "awk -F' ' '{print \$2}' $libids1 | sed -e"' "s|^|'$trimdir'/|g" > '"$reads" >> 2b.readalign.sh
echo "mapfile -t reads < $reads"' # ${reads[0]} calls read1 AND ${reads[1]} calls read2' >> 2b.readalign.sh
echo "awk -F' ' '{print \$2}' " $libids1 " | awk -F'_' '{print \$1\"_\"\$2\"_\"\$3}' > "$prefix "# create a prefix file to iterate" >> 2b.readalign.sh
echo "mapfile -t prefixmap < $prefix"' # assign prefixes to $prefixmap' >> 2b.readalign.sh
echo '# run bowtie2 with multimapping and threading, then output sorted BAM file' >> 2b.readalign.sh
echo 'srun bowtie2 -k '$multimapping' -X2000 --mm --threads '$bwt_thread' -x '$idx' -1 ${reads[0]} -2 ${reads[1]} 2>${prefixmap[0]}'$log '| samtools view -Su /dev/stdin | samtools sort -o ${prefixmap[0]}'$bam >> 2b.readalign.sh
printf '\n' >> 2b.readalign.sh
echo 'samtools flagstat ${prefixmap[0]}'"$bam > "'${prefixmap[0]}'"$fgQC1 # output alignment stats" >> 2b.readalign.sh


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
# echo "grep '_L001_R1.fastq.merged.gz' " $libids " | awk -F' '{print \$2}' | grep ^"$spID" > "$read1 "# create a list file of read1 for the array" >> 2b.readalign.sh
# echo "grep '_L001_R2.fastq.merged.gz' " $libids " | awk -F' '{print \$2}' | grep ^"$spID" > "$read2 "# create a list file of read2 for the array" >> 2b.readalign.sh
# echo "awk -F' ' '{print \$2}' " $libids " | awk -F'_' '{print \$1}' | grep "$spID" > "$prefix "# create a prefix file to iterate on the array" >> 2b.readalign.sh
# printf '\n' >> 2b.readalign.sh
# echo 'mapfile -t read1map < '$read1 '# assign read1 as elements to $read1map' >> 2b.readalign.sh
# echo 'mapfile -t read2map < '$read2 '# assign read2 as elements to $read2map' >> 2b.readalign.sh
# echo 'mapfile -t prefixmap < '$prefix '# assign prefixes to $prefixmap' >> 2b.readalign.sh
# printf '\n' >> 2b.readalign.sh
# echo '# run bowtie2 with multimapping and threading, then output sorted BAM file' >> 2b.readalign.sh
# echo 'srun bowtie2 -k' $multimapping '-X2000 --mm --threads' $bwt_thread '-x' $idx '-1 ${read1map[${SLURM_ARRAY_TASK_ID}]} -2 ${read2map[${SLURM_ARRAY_TASK_ID}]} 2>${prefixmap[${SLURM_ARRAY_TASK_ID}]}'$log '| samtools view -Su /dev/stdin | samtools sort -o ${prefixmap[${SLURM_ARRAY_TASK_ID}]}'$bam >> 2b.readalign.sh
# printf '\n' >> 2b.readalign.sh
# echo 'samtools flagstat ${prefixmap[${SLURM_ARRAY_TASK_ID}]}'$bam '>' $readalign'/${prefixmap[${SLURM_ARRAY_TASK_ID}]}'$fgQC1 '# output alignment stats' >> 2b.readalign.sh

# JOBID4=$( sbatch -W --dependency=afterok:${JOBID3} --array=$RAarray 2b.readalign.sh | awk '{print $4}' ) # JOB4 depends on JOB3 completing successfully

echo '# -- 2a.'$spID' genome index building completed -- #'

echo '# -- 2b.'$spID' read alignment started -- #'

JOBID4=$( sbatch -W --dependency=afterok:${JOBID3} 2b.readalign.sh | awk '{print $4}' ) # JOB4 depends on JOB3 completing successfully



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
echo 'ssh -o StrictHostKeyChecking=no -l ${USERNAME} ${HOSTNAME} "${SCRIPT}"' >> 3.mtfilt_fragcount_A.sh

echo '# -- 2b.'$spID' read alignment completed -- #'

echo '# -- 3a.'$spID' mitochondrial genome download started: '$mtID' -- #'

JOBID5=$( sbatch -W --dependency=afterok:${JOBID4} 3.mtfilt_fragcount_A.sh | awk '{print $4}' ) # JOB5 depends on JOB4 completing successfully

## 3b. remove mitochondrial mapped reads from ATAC data and plot fragment lengths
echo '#!/bin/bash -e' > 3.mtfilt_fragcount_B.sh
echo '#SBATCH -p ei-medium # partition (queue)' >> 3.mtfilt_fragcount_B.sh
echo '#SBATCH -N 1 # number of nodes' >> 3.mtfilt_fragcount_B.sh
echo '#SBATCH -n 1 # number of tasks' >> 3.mtfilt_fragcount_B.sh
echo '#SBATCH --mem 48000' >> 3.mtfilt_fragcount_B.sh
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
echo '# if-else statement - if $mtscaff exists then continue with normal steps (filtering etc.), else if it does not exist then create a new file with same name $bam_file_nochrM' >> 3.mtfilt_fragcount_B.sh
echo "if [ -f $mtscaff ]; then" >> 3.mtfilt_fragcount_B.sh
echo -e '\techo '"$mtscaff" ' DOES exist - filtering chrM mapped reads, outputing new BAM and plotting fragment length count' >> 3.mtfilt_fragcount_B.sh
echo -e '\t# 5. Store the genome scaffolds names that are mitochondrial genome under variable $scaffarray and create grep command $grepscaff' >> 3.mtfilt_fragcount_B.sh
echo -e "\tIFS=$'\\\n' scaffarray=("'$(cut -f2 '$mtscaff" | awk ""'"'!x[$0]++'"'))" '# this takes the $mtscaff as input, takes unique genome scaffolds (col2) that match mtDNA (using awk instead of sort -u so top hit ordering retained) and then assigns to the variable $scaffarray. Accessed using ${scaffarray[0]}..${scaffarray[2]}' >> 3.mtfilt_fragcount_B.sh
echo -e '\techo Scaffold/s matching mitochondrial genome are/is: ${scaffarray[@]} # this will echo each scaffold that matches mtDNA' >> 3.mtfilt_fragcount_B.sh
echo -e '\techo Total number of scaffolds matching mitochondrial genome is: ${#scaffarray[@]} # this will echo the number of scaffolds that match mtDNA' >> 3.mtfilt_fragcount_B.sh
echo -e '\t# for sc in ${scaffarray[@]}; do echo $sc ; done # this will loop through each element (Scaff) of the array and echo on new line' >> 3.mtfilt_fragcount_B.sh
echo -e '\tgrepscaff=$(echo ${scaffarray[@]} | sed '"'s/ / | grep -v /g' | sed 's/^/grep -v /') # Assign grep commands to variable: input grep -v and then pipe in front of each element of the array e.g. UNK2407 UNK3173 UNK1343 > grep -v UNK2407 | grep -v UNK3173 | grep -v UNK1343" >> 3.mtfilt_fragcount_B.sh
echo -e '\t# 6. Loop through variable of scaffold/s matching the mitochondrial genome and filter reads assigned to these scaffolds.' >> 3.mtfilt_fragcount_B.sh
echo -e '\tml samtools/1.3' >> 3.mtfilt_fragcount_B.sh
echo -e '\tml perl_activeperl/5.18' >> 3.mtfilt_fragcount_B.sh
echo -e '\tml zlib/1.2.8' >> 3.mtfilt_fragcount_B.sh
echo -e '\t# A. assign mtDNA filtered bam to new file with extension .nochrM.bam' >> 3.mtfilt_fragcount_B.sh
echo -e '\t# B. remove all reads mapping to mitochondrial scaffold/s and index' >> 3.mtfilt_fragcount_B.sh
echo -e "\tfor bam_file in $readalign/*.bam; do bam_file_nochrM="'$(echo $bam_file | sed -e '"'s/.bam/.nochrM.bam/' | sed -e 's/2.read_alignment/3.Mtfilt_fragcnt/g'); samtools idxstats "'$bam_file | cut -f1 | $grepscaff | xargs samtools view -b $bam_file > $bam_file_nochrM; samtools index $bam_file_nochrM; done #samtools view -h $bam_file | grep ${scaffarray[0]} | wc -l # you can type this to test' >> 3.mtfilt_fragcount_B.sh
echo -e '\t# 7. Calculate fragment length count for each' >> 3.mtfilt_fragcount_B.sh
echo -e '\tfrag_length=$(echo $bam_file_nochrM | sed -e '"'s/.bam/_frag_length_count.txt/')" >> 3.mtfilt_fragcount_B.sh
echo -e '\tsamtools view $bam_file_nochrM | awk '"'"'$9>0'"' | cut -f9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > "'$frag_length' >> 3.mtfilt_fragcount_B.sh
echo 'else' >> 3.mtfilt_fragcount_B.sh
echo -e '\techo '"$mtscaff"' DOES NOT exist - creating a non-chrM_filtered BAM file and plotting fragment length count' >> 3.mtfilt_fragcount_B.sh
echo -e '\t# A. assign non-mtDNA filtered bam to new file with extension .nochrM.bam (for complete naming conventions)' >> 3.mtfilt_fragcount_B.sh
echo -e "\tfor bam_file in $readalign/*.bam; do bam_file_nochrM="'$(echo $bam_file | sed -e '"'s/.bam/.nochrM.bam/' | sed -e 's/2.read_alignment/3.Mtfilt_fragcnt/g'); xargs samtools view -b "'$bam_file > $bam_file_nochrM; samtools index $bam_file_nochrM; done' >> 3.mtfilt_fragcount_B.sh
echo 'fi' >> 3.mtfilt_fragcount_B.sh

# ml samtools/1.7
# # 2. create a blast database of each genome and store under variables
# ml blast/2.3.0
# makeblastdb -in $gFA -parse_seqids -dbtype nucl # create a blast database of the genome assembly
# # 3. blast the mitochondrial genomes against the assembly - output tabular format
# blastn -db $gFA -outfmt 6 -evalue 1e-3 -word_size 11 -show_gis -num_alignments 10 -max_hsps 20 -num_threads 5 -out Ab5_L_ATAC.genome_mt.blast -query /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/data/ATACseq/3.run2/Ab5_L_ATAC/3.Mtfilt_fragcnt/NC_027289.1.fasta # blast the mitochondrial genome against the input genome assembly, output tabular format
# # 4. Python script: filter BLAST hits based on pident>=93 and evalue<=1e-10; then pident>75% according to BLAST hit and alignment length to mitochondrial genome
# ml python/3.5
# python $scripts/ATAC_Bioinf_pipeline_v2b_part3a.py $blstout $mtgen # output is stored as variable $mtscaff at top
# # if-else statement - if $mtscaff exists then continue with normal steps (filtering etc.), else if it does not exist then create a new file with same name '$bam_file_nochrM'
# if [ -f $mtscaff ]; then
#     echo "$mtscaff DOES exist - filtering chrM mapped reads, outputing new BAM and plotting fragment length count"
#     # 5. Store the genome scaffolds names that are mitochondrial genome under variable $scaffarray and create grep command $grepscaff
#     IFS=$'\n' scaffarray=($(cut -f2 $mtscaff | awk '!x[$0]++')) # this takes the $mtscaff as input, takes unique genome scaffolds (col2) that match mtDNA (using awk instead of sort -u so top hit ordering retained) and then assigns to the variable $scaffarray. Accessed using ${scaffarray[0]}..${scaffarray[2]}
#     echo Scaffold/s matching mitochondrial genome are/is: ${scaffarray[@]} # this will echo each scaffold that matches mtDNA
#     echo Total number of scaffolds matching mitochondrial genome is: ${#scaffarray[@]} # this will echo the number of scaffolds that match mtDNA
#     # for sc in ${scaffarray[@]}; do echo $sc ; done # this will loop through each element (Scaff) of the array and echo on new line
#     grepscaff=$(echo ${scaffarray[@]} | sed 's/ / | grep -v /g' | sed 's/^/grep -v /') # Assign grep commands to variable: input grep -v and then pipe in front of each element of the array e.g. UNK2407 UNK3173 UNK1343 > grep -v UNK2407 | grep -v UNK3173 | grep -v UNK1343
#     # 6. Loop through variable of scaffold/s matching the mitochondrial genome and filter reads assigned to these scaffolds.
#     ml samtools/1.3
#     ml perl_activeperl/5.18
#     ml zlib/1.2.8
#     # A. assign mtDNA filtered bam to new file with extension .nochrM.bam
#     # B. remove all reads mapping to mitochondrial scaffold/s and index
#     for bam_file in $readalign/*.bam; do bam_file_nochrM=$(echo $bam_file | sed -e 's/.bam/.nochrM.bam/' | sed -e 's/2.read_alignment/3.Mtfilt_fragcnt/g'); samtools idxstats $bam_file | cut -f1 | $grepscaff | xargs samtools view -b $bam_file > $bam_file_nochrM; samtools index $bam_file_nochrM; done #samtools view -h $bam_file | grep ${scaffarray[0]} | wc -l # you can type this to test
#     # 7. Calculate fragment length count for each
#     frag_length=$(echo $bam_file_nochrM | sed -e 's/.bam/_frag_length_count.txt/')
#     samtools view $bam_file_nochrM | awk '$9>0' | cut -f9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > $frag_length
# else
#     echo "$mtscaff DOES NOT exist - creating a non-chrM_filtered BAM file and plotting fragment length count"
#     # A. assign non-mtDNA filtered bam to new file with extension .nochrM.bam (for complete naming conventions)
#     for bam_file in $readalign/*.bam; do bam_file_nochrM=$(echo $bam_file | sed -e 's/.bam/.nochrM.bam/' | sed -e 's/2.read_alignment/3.Mtfilt_fragcnt/g'); xargs samtools view -b $bam_file > $bam_file_nochrM; samtools index $bam_file_nochrM; done
# fi
# 8. Plot fragment length count in R
# source R-3.5.2
# bam_file_nochrM=($spID.nochrM.bam)
# R CMD BATCH --no-save --no-restore '--args $bam_file_nochrM' $scripts/ATAC_Bioinf_pipeline_v2b_part3b.R ATAC_Bioinf_pipeline_v2b_part3b.Rout # this creates two files - Rplots.pdf (which has the image!) and another (empty) image file with the actual filename. Simply rename Rplots.pdf
# mv Rplots.pdf "$(basename "$bam_file_nochrM" .bam).fraglength.pdf" # rename Rplots.pdf to *.fraglength.pdf

echo '# -- 3a.'$spID' mitochondrial genome downloaded: '$mtID' -- #'

echo '# -- 3b.'$spID' mitochondrial removed and fragment dist plot started -- #'

JOBID6=$( sbatch -W --dependency=afterok:${JOBID5} 3.mtfilt_fragcount_B.sh | awk '{print $4}' ) # JOB6 depends on JOB5 completing successfully

# 3b-A. Plot fragment lengths - note that this will produce a ggplot error and there is no STDOUT but the file will be Rplots.pdf
echo '#!/bin/bash -e' > 3.mtfilt_fragcount_B_a.sh
echo '#SBATCH -p ei-medium # partition (queue)' >> 3.mtfilt_fragcount_B_a.sh
echo '#SBATCH -N 1 # number of nodes' >> 3.mtfilt_fragcount_B_a.sh
echo '#SBATCH -n 1 # number of tasks' >> 3.mtfilt_fragcount_B_a.sh
echo '#SBATCH --mem 24000' >> 3.mtfilt_fragcount_B_a.sh
echo '#SBATCH -t 0-04:59' >> 3.mtfilt_fragcount_B_a.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 3.mtfilt_fragcount_B_a.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 3.mtfilt_fragcount_B_a.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 3.mtfilt_fragcount_B_a.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 3.mtfilt_fragcount_B_a.sh
printf '\n' >> 3.mtfilt_fragcount_B_a.sh
echo '# 8. Plot fragment length count in R' >> 3.mtfilt_fragcount_B_a.sh
echo 'source R-3.5.2' >> 3.mtfilt_fragcount_B_a.sh
echo "R CMD BATCH --no-save --no-restore '--args $bam_file_nochrM' $scripts/ATAC_Bioinf_pipeline_v2b_part3b.R ATAC_Bioinf_pipeline_v2b_part3b.Rout # this creates two files - Rplots.pdf (which has the image!) and another (empty) image file with the actual filename. Simply rename Rplots.pdf" >> 3.mtfilt_fragcount_B_a.sh

JOBID6a=$( sbatch -W --dependency=afterok:${JOBID6} 3.mtfilt_fragcount_B_a.sh | awk '{print $4}' ) # JOB6a depends on JOB6 completing successfully - NOTE: this will not be used as a dependency since it spits out a ggplot error despite running fine!

################################################################################################################

### 4. Post alignment filtering
#   4a. Filter reads (Sort, Map and remove duplicates, Remove reads unmapped, not primary alignment, reads failing platform, duplicates) - samtools > final bam


mkdir -p $filtdir
cd $filtdir
mkdir $filtdir/tmp

echo '#!/bin/bash -e' > 4.postalign_filt.sh
echo '#SBATCH -p ei-medium # partition (queue)' >> 4.postalign_filt.sh
echo '#SBATCH -N 1 # number of nodes' >> 4.postalign_filt.sh
echo '#SBATCH --mem 120000' >> 4.postalign_filt.sh
echo '#SBATCH -t 0-10:59' >> 4.postalign_filt.sh
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
echo 'ml bedtools/2.25.0' >> 4.postalign_filt.sh
echo 'ml zlib/1.2.8' >> 4.postalign_filt.sh
echo 'ml glib/2.40' >> 4.postalign_filt.sh
echo 'ml R/3.2.3' >> 4.postalign_filt.sh
printf '\n' >> 4.postalign_filt.sh
echo "for bam_file in ${mtfilt}/*.nochrM.bam; do" >> 4.postalign_filt.sh
echo -e '\t# variables for output files' >> 4.postalign_filt.sh
echo -e '\tbam_file_sorted=$(echo $bam_file | sed -e '"'s/.bam/.sorted.bam/' | sed -e 's/3.Mtfilt_fragcnt/4.postalign_filt/g')" >> 4.postalign_filt.sh
echo -e '\tbam_file_dup=$(echo $bam_file | sed -e '"'s/.bam/.sorted.dup.bam/' | sed -e 's/3.Mtfilt_fragcnt/4.postalign_filt/g')" >> 4.postalign_filt.sh
echo -e '\tnodup_filt_bam_file=$(echo $bam_file | sed -e '"'s/.bam/.nodup.filt.bam/' | sed -e 's/3.Mtfilt_fragcnt/4.postalign_filt/g') # final bam file" >> 4.postalign_filt.sh
echo -e '\tnodup_filt_bam_index_file=$(echo $bam_file | sed -e '"'s/.bam/.nodup.filt.bam.bai/' | sed -e 's/3.Mtfilt_fragcnt/4.postalign_filt/g') # index file" >> 4.postalign_filt.sh
echo -e '\tnodup_filt_bam_file_mapstats=$(echo $bam_file | sed -e '"'s/.bam/.flagstat.qc/' | sed -e 's/3.Mtfilt_fragcnt/4.postalign_filt/g') # QC file" >> 4.postalign_filt.sh
echo -e '\tpbc_file_qc=$(echo $bam_file | sed -e '"'s/.bam/.pbc.qc/' | sed -e 's/3.Mtfilt_fragcnt/4.postalign_filt/g') # library complexity" >> 4.postalign_filt.sh
echo -e '\tnodup_filt_bam_file_sorted=$(echo $bam_file | sed -e '"'s/.bam/.srt.bam/' | sed -e 's/3.Mtfilt_fragcnt/4.postalign_filt/g') # final bam file, sorted (temp)" >> 4.postalign_filt.sh
echo -e '\t# Filter reads' >> 4.postalign_filt.sh
echo -e "\tsambamba sort -m 88G -t 1 --tmpdir tmp -o "'$bam_file_sorted -u $bam_file' >> 4.postalign_filt.sh
echo -e "\tsambamba markdup -l 0 -t 1 "'$bam_file_sorted '""'$bam_file_dup' >> 4.postalign_filt.sh
echo -e "\tsamtools view -F 1804 -f 2 -q 30 -b "'$bam_file_dup > '""'$nodup_filt_bam_file' >> 4.postalign_filt.sh
echo -e "\tsamtools index "'$nodup_filt_bam_file '""'$nodup_filt_bam_index_file' >> 4.postalign_filt.sh
echo -e "\tsamtools flagstat "'$nodup_filt_bam_file > '""'$nodup_filt_bam_file_mapstats' >> 4.postalign_filt.sh
echo -e '\t# Plotting the fragment length distribution' >> 4.postalign_filt.sh
echo -e "\tjava -Xmx2g -jar /tgac/software/production/picardtools/1.84/x86_64/bin/CollectInsertSizeMetrics.jar \R=$gFA \I="'$nodup_filt_bam_file \O="'""'${nodup_filt_bam_file}_PicardInsertMetrics.jar.txt" \H="'""'${nodup_filt_bam_file}_insert_size_histogram.pdf" \M=0.5' >> 4.postalign_filt.sh
echo 'done' >> 4.postalign_filt.sh

# for bam_file in ${mtfilt}/*.nochrM.bam; do
#   # variables for output files
#   bam_file_sorted=$(echo $bam_file | sed -e 's/.bam/.sorted.bam/' | sed -e 's/3.Mtfilt_fragcnt/4.postalign_filt/g')
#   bam_file_dup=$(echo $bam_file | sed -e 's/.bam/.sorted.dup.bam/' | sed -e 's/3.Mtfilt_fragcnt/4.postalign_filt/g')
#   nodup_filt_bam_file=$(echo $bam_file | sed -e 's/.bam/.nodup.filt.bam/' | sed -e 's/3.Mtfilt_fragcnt/4.postalign_filt/g') # final bam file
#   nodup_filt_bam_index_file=$(echo $bam_file | sed -e 's/.bam/.nodup.filt.bam.bai/' | sed -e 's/3.Mtfilt_fragcnt/4.postalign_filt/g') # index file
#   nodup_filt_bam_file_mapstats=$(echo $bam_file | sed -e 's/.bam/.flagstat.qc/' | sed -e 's/3.Mtfilt_fragcnt/4.postalign_filt/g') # QC file
#   pbc_file_qc=$(echo $bam_file | sed -e 's/.bam/.pbc.qc/' | sed -e 's/3.Mtfilt_fragcnt/4.postalign_filt/g') # library complexity
#   nodup_filt_bam_file_sorted=$(echo $bam_file | sed -e 's/.bam/.srt.bam/' | sed -e 's/3.Mtfilt_fragcnt/4.postalign_filt/g') # final bam file, sorted (temp)
#   # Filter reads
#   sambamba sort -m 54G -t 2 -o $bam_file_sorted -u $bam_file
#   sambamba markdup -l 0 -t 2 $bam_file_sorted $bam_file_dup
#   samtools view -F 1804 -f 2 -q 30 -b $bam_file_dup > $nodup_filt_bam_file
#   samtools index $nodup_filt_bam_file $nodup_filt_bam_index_file
#   samtools flagstat $nodup_filt_bam_file > $nodup_filt_bam_file_mapstats
#   # Plotting the fragment length distribution
#   java -Xmx2g -jar /tgac/software/production/picardtools/1.84/x86_64/bin/CollectInsertSizeMetrics.jar \R=$gFA \I=$nodup_filt_bam_file \O="${nodup_filt_bam_file}_PicardInsertMetrics.jar.txt" \H="${nodup_filt_bam_file}_insert_size_histogram.pdf" \M=0.5
# done

echo '# -- 3b.'$spID' mitochondrial removed and fragment dist plot completed -- #'

mv $mtfilt/Rplots.pdf $mtfilt/"$(basename "$bam_file_nochrM" .bam).fraglength.pdf" # rename Rplots.pdf to *.fraglength.pdf

echo '# -- 4.'$spID' Post alignment filtering started -- #'

JOBID7=$( sbatch -W --dependency=afterok:${JOBID6} 4.postalign_filt.sh | awk '{print $4}' ) # JOB7 depends on JOB6 completing successfully

################################################################################################################

### 5. ATAC peak calling (test vs control)
# 	5a. Convert PE Bam to tagalign (start/end positions of each read) - bedtools
# 	5b. TN5 shifting of tagaligns - shift reads +4 bp for the +strand and -5 bp for the -strand
# 	5c. count-based peak calling using Poisson distribution - macs2
#   5d. count-based peak calling using another program that considers properly paired, unpaired and secondary alignments (unlike MACS2) - Genrich
#   5e. markov model based peak calling specific for ATAC-seq data - HMMRATAC

mkdir -p $peakcall
cd $peakcall

echo '#!/bin/bash -e' > 5.peakcall.sh
echo '#SBATCH -p ei-medium # partition (queue)' >> 5.peakcall.sh
echo '#SBATCH -N 1 # number of nodes' >> 5.peakcall.sh
echo '#SBATCH --mem 48000' >> 5.peakcall.sh
echo '#SBATCH -t 0-12:59' >> 5.peakcall.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 5.peakcall.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 5.peakcall.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 5.peakcall.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo 'ml MACS' >> 5.peakcall.sh
# echo 'ml python/3.5' >> 5.peakcall.sh
echo 'ml bedtools/2.25.0' >> 5.peakcall.sh
echo 'ml GCC' >> 5.peakcall.sh
echo 'ml zlib' >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo 'ulimit -Sn 10000 # set the open file limit to 10,000' >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo '# Sum total length of all scaffolds in assembly - needed as input for Macs2' >> 5.peakcall.sh
echo 'source bioawk-1.0' >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo 'GenSz=$(bioawk -c fastx'" '{ print "'$name, length($seq) }'"' < $gFA | awk "'-F"\t" '"'{print;x+="'$2}END{print "Total " x}'"' | tail -1 | sed 's/Total //g') # this will output the sum length of scaffolds and assign to variable "'$Gensz' >> 5.peakcall.sh
echo 'bioawk -c fastx'" '{ print "'$name, length($seq) }'"' < $gFA | awk '{print "'$1,"0",$2}'"' OFS="'"\t" > '"$scafflen" >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo "# 5a. Convert each bam to a .tagAlign - contains the start/end positions of each read:" >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo "bedtools bamtobed -i $Test1 | awk 'BEGIN{OFS="'"\t"}{$4="N";$5="1000";print $0}'"' | gzip -c  > $tagalign_test1" >> 5.peakcall.sh
echo "bedtools bamtobed -i $Control1 | awk 'BEGIN{OFS="'"\t"}{$4="N";$5="1000";print $0}'"' | gzip -c > $tagalign_control1" >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo '# 5c. Tn5 shifting of tagaligns' >> 5.peakcall.sh
echo "zcat $tagalign_test1 | awk -F "'$"\t" '"'BEGIN {OFS = FS}{ if ("'$6 == "+") {$2 = $2 + 4} else if ($6 == "-") {$3 = $3 - 5} print $0}'"' | gzip -c > $shifted_tag" >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo '# 5d. Run MACS2:' >> 5.peakcall.sh
echo "macs2 callpeak -t $shifted_tag -c $tagalign_control1 -f BED -n $output_prefix -g "'$GenSz -p '"$Macs2PvalThresh --nomodel --shift -$Macs2ShiftSize --extsize $Macs2SmoothWindow -B --SPMR --keep-dup all --call-summits" >> 5.peakcall.sh
echo '# --nomodel and --extsize 150 tells MACS2 to use 150bp as fragment size to pileup sequencing reads.' >> 5.peakcall.sh
echo '# -g XX lets MACS2 consider a genome size as background.' >> 5.peakcall.sh
echo '# -B --SPMR ask MACS2 to generate pileup signal file of fragment pileup per million reads in bedGraph format.' >> 5.peakcall.sh
printf '\n' >> 5.peakcall.sh
echo '# Generate a fold change file comparing the sample to the control and logLR track' >> 5.peakcall.sh
echo "macs2 bdgcmp -t $output_prefix\_treat_pileup.bdg -c $output_prefix\_control_lambda.bdg -o $output_prefix\_FE.bdg -m FE" >> 5.peakcall.sh
echo "macs2 bdgcmp -t $output_prefix\_treat_pileup.bdg -c $output_prefix\_control_lambda.bdg -o $output_prefix\_logLR.bdg -m logLR -p 0.00001" >> 5.peakcall.sh
echo '# 5e. peak calling using another program - Genrich (installed on HPC by CiS) and take the intersection' >> 5.peakcall.sh
echo 'ml samtools/1.3' >> 5.peakcall.sh
echo "samtools sort -o $Test2 -O bam -n $Test1" >> 5.peakcall.sh
echo "samtools sort -o $Control2 -O bam -n $Control1" >> 5.peakcall.sh
echo 'source package 8bf6d6cb-b9e9-4215-a3d8-b17a76fec816 # sources Genrich on HPC - guidance: https://informatics.fas.harvard.edu/atac-seq-guidelines.html#another-peak-caller-why; https://github.com/jsh58/Genrich' >> 5.peakcall.sh
echo "Genrich -t $Test2 -c $Control2 -o $output_prefix'_Genrich.peaks' -p $Macs2PvalThresh -j -y -r -v # output is ENCODE narrowPeak format" >> 5.peakcall.sh

# ENCODE narrowPeak format
# 1. chrom 	Name of the chromosome
# 2. chromStart 	Starting position of the peak (0-based)
# 3. chromEnd 	Ending position of the peak (not inclusive)
# 4. name 	peak_N, where N is the 0-based count
# 5. score 	Average AUC (total AUC / bp) × 1000, rounded to the nearest int (max. 1000)
# 6. strand 	. (no orientation)
# 7. signalValue 	Total area under the curve (AUC)
# 8. pValue 	Summit -log10(p-value)
# 9. qValue 	Summit -log10(q-value), or -1 if not available (e.g. without -q)
# 10. peak 	Summit position (0-based offset from chromStart): the midpoint of the peak interval with the highest significance (the longest interval in case of ties)

echo '#!/bin/bash -e' > 5f.peakcall.sh
echo '#SBATCH -p ei-medium # partition (queue)' >> 5f.peakcall.sh
echo '#SBATCH -N 1 # number of nodes' >> 5f.peakcall.sh
echo '#SBATCH --mem 48000' >> 5f.peakcall.sh
echo '#SBATCH -t 0-04:59' >> 5f.peakcall.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 5f.peakcall.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 5f.peakcall.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 5f.peakcall.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 5f.peakcall.sh
printf '\n' >> 5f.peakcall.sh
echo '# 5f. Markov model based peak caller: HMMRATAC - a three-state semi-supervised hidden Markov model (HMM) to simultaneously segment the genome into open chromatin regions with high signal, nucleosomal regions with moderate signals, and background regions with low signals, respectively' >> 5f.peakcall.sh
printf '\n' >> 5f.peakcall.sh
echo 'source package 3352c4a6-04ca-4bcd-a3e6-9ace9fa757b9' >> 5f.peakcall.sh
echo 'ml samtools/1.3' >> 5f.peakcall.sh
echo '# samtools view -H $Test1 | perl -ne'" 'if(/^@SQ.*?SN:(\w+)\s+LN:(\d+)/){print "'$1,"\t",$2,"\n"}'"' > "'$Geninfo # Make genome information (chromosome sizes) from the BAM file to get a genome.info file' >> 5f.peakcall.sh
echo "cut -f1,3 ${scafflen} > ${scafflen2}" >> 5f.peakcall.sh
echo "HMMRATAC -Xms20g -Xmx80g -XX:ParallelGCThreads=2 -b $Test1 -i $Test1index -g $scafflen2 -o $spID # Run HMMRATAC" >> 5f.peakcall.sh
echo 'awk -v OFS="'"\t"'" '"'"'$13>=10 {print}'"' $spID"'_peaks.gappedPeak'" > $spID.filteredPeaks.gappedPeak # Filter HMMRATAC output by the score, if desired. Score threshold will depend on dataset, score type and user preference. A threshold of 10 would be:" >> 5f.peakcall.sh
echo 'awk -v OFS="\t" '"'"'$5>=10 {print}'"' $spID"'_summits.bed > '"$spID.filteredSummits.bed # filter the summit file by the same threshold" >> 5f.peakcall.sh

echo '# -- 4.'$spID' Post alignment filtering completed -- #'

echo '# -- 5.'$spID' Peak calling started -- #'

JOBID8=$( sbatch -W --dependency=afterok:${JOBID7} 5.peakcall.sh | awk '{print $4}' ) # JOB8 depends on JOB7 completing successfully
# JOBID9=$( sbatch -W --dependency=afterok:${JOBID7} 5f.peakcall.sh | awk '{print $4}' ) # JOB9 depends on JOB7 completing successfully

## Ran three peakcallers for the following purpose:
# 1. We will ultimately use the MACS2 narrow peak file for peaks since this is what is recommended by:
  # A. ENCODE project’s “ATAC-seq Data Standards and Prototype Processing Pipeline”
  # B. https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1929-3
# 2. Will perform (in a separate script), IDR as per ENCODE project’s “ATAC-seq Data Standards and Prototype Processing Pipeline” for true replicated data
# 3. The three peak callers we used are two count-based (MACS2 and Genrich, developed by the ATAC-seq experts) and a Markov Model based caller developed specifically for ATAC-seq data
# 4. We will keep data of peaks from all three peak callers for test purposes and can call intersections for reviewers if required - maybe ensure that peaks that are focused on are at least present in 2/3 sets.

################################################################################################################

### 6. bed to bigbed conversion for narrowpeaks - bedClip and bedToBigbed from ucsc_tools

# Peaks have been called with MACS2 > narrow peaks file: This is the basic dataset

cd $peakcall

echo '# -- 5.'$spID' Peak calling completed -- #'

echo '# -- 6.'$spID' Bed to BigBed conversion started -- #'

gzip $peak # gzip compress the narrowPeak file

source ucsc_utils-v333

echo 'table narrowPeak' > narrowPeak.as
echo '"BED6+4 Peaks of signal enrichment based on pooled, normalized (interpreted) data."' >> narrowPeak.as
echo '(' >> narrowPeak.as
echo -e '\tstring chrom;        "Reference sequence chromosome or scaffold"' >> narrowPeak.as
echo -e '\tuint   chromStart;   "Start position in chromosome"' >> narrowPeak.as
echo -e '\tuint   chromEnd;     "End position in chromosome"' >> narrowPeak.as
echo -e '\tstring name;	 "Name given to a region (preferably unique). Use . if no name is assigned"' >> narrowPeak.as
echo -e '\tuint   score;        "Indicates how dark the peak will be displayed in the browser (0-1000) "' >> narrowPeak.as
echo -e '\tchar[1]  strand;     "+ or - or . for unknown"' >> narrowPeak.as
echo -e '\tfloat  signalValue;  "Measurement of average enrichment for the region"' >> narrowPeak.as
echo -e '\tfloat  pValue;       "Statistical significance of signal value (-log10). Set to -1 if not used."' >> narrowPeak.as
echo -e '\tfloat  qValue;       "Statistical significance with multiple-test correction applied (FDR -log10). Set to -1 if not used."' >> narrowPeak.as
echo -e '\tint   peak;         "Point-source called for this peak; 0-based offset from chromStart. Set to -1 if no point-source called."' >> narrowPeak.as
echo ')' >> narrowPeak.as

# cut -f1,3 ${scafflen} > ${scafflen2}
zcat ${peakgz} | sort -k1,1 -k2,2n > ${bigbed}.tmp
bedClip ${bigbed}.tmp ${scafflen2} ${bigbed}.tmp2

bedToBigBed -type=bed6+4 -as=narrowPeak.as ${bigbed}.tmp2 ${scafflen2} ${bigbed}
rm -f ${bigbed}.tmp ${bigbed}.tmp2

echo '# -- 6.'$spID' Bed to BigBed conversion completed -- #'

################################################################################################################

### Finish the script
echo -e '# --------------------\nEXITING SCRIPT - all completed\n# --------------------'

exit 0
