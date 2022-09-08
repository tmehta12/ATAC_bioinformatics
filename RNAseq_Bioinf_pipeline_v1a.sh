#!/bin/sh

#!/bin/bash -e
#SBATCH -p tgac-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 4000 # memory pool for all cores
#SBATCH -t 1-23:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

################################################################################################################

# RNA-seq pipeline - Part 1
# July 2021: Tarang K. Mehta, Earlham Institute, Norwich, UK

################################################################################################################

# ~ This pipeline will be ran as species-specific, and contains the following components:

# 1. Trimming and renaming - Trim Galore
# 2. Mapping to genome (can also map to transcriptome but opting for less bias this way) - HiSat2
# 3. Alignment QC - Qualimap2
# 4. Counts - HTseq
# 5. Differential gene expression - DESeq2
  # Outputs: normalized quantification tables, some important statistics for the whole gene or transcript list, and the list of significantly differentially expressed genes or transcripts (with default threshold of FDR<0.05)
  # The raw count is normalized based on Trimmed Mean of M values (TMM) (if edgeR is used) or the median-of-ratios method (if DESeq2 is used) when the reads are mapped to a genome.

################################################################################################################

# Script usage: ./RNAseq_Bioinf_pipeline_v1a.sh -s "spID" -g "spG" -f "gFA" -m "mtID" -u "Usr" -a "annot"
# e.g. ./RNAseq_Bioinf_pipeline_v1a.sh -s Mz_1_Brain -g M_zebra_UMD1 -f /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -u mehtat -a M_zebra_UMD1.gtf
# Note: Script is adapted for SBATCH usage

## Inputs:
# 1. Raw fastq.gz file locations
# 2. Target species genome FASTA
# 3. Target species genome annotation gtf

## Place this script and the following files in $wd
# 1. $wd/libids.txt: 2-col tab-delimited file where col1 is the *_{R1,R2}.fastq.gz and col2 is the species_tissueID_experiment_barcode_{R1,R2}.fastq.gz (for renaming purposes)
  # This was created for my files like so:
  # ls -1 $readsdir > col1.txt
  # # create col2 - species_individual_tissueID_barcode_{R1,R2}.fastq.gz
  # ls -1 $readsdir | awk -F'_' '{print $2}' | sed 's|Mz|_Mz_|g' | sed 's|Pnm|_Pn_|g' | sed 's|Pmn|_Pn_|g' | sed 's|Ab|_Ab_|g' | sed 's|Nb|_Nb_|g' | sed 's|On|_On_|g' | awk -F'_' '{print $2}' > sp.txt # species
  # ls -1 $readsdir | awk -F'_' '{print $2}' | sed 's|Mz|_Mz_|g' | sed 's|Pnm|_Pnm_|g' | sed 's|Pmn|_Pmn_|g' | sed 's|Ab|_Ab_|g' | sed 's|Nb|_Nb_|g' | sed 's|On|_On_|g' | awk -F'_' '{print $3}' | sed 's|B|Brain|g' | sed 's|E|Eye|g' | sed 's|L|Liver|g' | sed 's|T|Testis|g' | sed 's|G|Gill|g' | sed 's/[0-9]\{1,\}/_&_/' | sed 's/[0-9]\+$//' | sed 's|$|_|g' > ind_tissue.txt # _individual_tissueID_
  # ls -1 $readsdir | sed 's|_1_1_|_1_|g' | awk -F'_' '{print $6}' | sed 's|$|_|g' > bc.txt # barcode
  # ls -1 $readsdir | sed 's|_1_1_|_1_|g' | awk -F'_' '{print $8}' > ext.txt # read extension
  # paste -d '' sp.txt ind_tissue.txt bc.txt ext.txt > col2.txt
  # paste -d'\t' col1.txt col2.txt > libids.txt

################################################################################################################

# Setting parameters for command line input

helpFunction()
{
   echo ""
   echo "Usage: $0 -s spID -g spG -f gFA -u Usr -t transcr"
   echo -e "\t-s spID = Species ID, preferably two short letters '_' individual '_' tissue e.g. P. nyererei Individual 1 Brain = Pn_1_Brain Note: this naming convention needs to be the same as renaming in 2-column tab-delimited file"
   echo -e "\t-g spG = Species genome ID e.g. hg19 or M_zebra_UMD1"
   echo -e "\t-f gFA = Full path to genome FASTA file"
   echo -e "\t-u Usr = ssh username e.g. mehtat"
   echo -e "\t-t transcr = Full path to genome annotations: either 1) gzipped GTF file (*gtf.gz); or 2) longest protein-coding gene annotations as 6-column BED: col1-scaff,col2-start,col3-end,col4-geneID,col5-XX,col6-strand"
   exit 1 # Exit script after printing help
}

while getopts "s:g:f:u:t:" opt
do
   case "$opt" in
      s ) spID="$OPTARG" ;;
      g ) spG="$OPTARG" ;;
      f ) gFA="$OPTARG" ;;
      u ) Usr="$OPTARG" ;;
      t ) transcr="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$spID" ] || [ -z "$spG" ] || [ -z "$gFA" ] || [ -z "$Usr" ] || [ -z "$transcr" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "$spID"
echo "$spG"
echo "$gFA"
echo "$Usr"
echo "$transcr"

################################################################################################################

# Variables

topdir=(/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch) # place all scripts in the topmost directory - create this separately
WD=(${topdir}/RNAseq)
wd=(${topdir}/RNAseq/$spID) # insert the working directory
email=Tarang.Mehta@earlham.ac.uk # SBATCH out and err send to address

### 0. Reads
readsgp=(/ei/data/reads/PIP-2735/210607_A00478_0165_AHF53CDSX2) # manual lib prep reads - 8 samples
readsgp2=(/ei/data/reads/PIP-2734/210607_A00478_0165_AHF53CDSX2) # automated lib prep reads - 37 samples
readsdirmain=(${WD}/0.reads)
readsdir=(${wd}/0.rawreads)

### 1. Trimming and renaming - Trim Galore
trimdir=($wd/1.adaptor_trimming) # assign trimmed reads dir
trimarray=0-0 # INSERT the number range of paired *fastq.gz to trim in zero base e.g. 10 pairs = 0-9
libids=($WD/libids.txt) # 2-col space-delimited file where col1 is the *_{R1,R2}.fastq.gz and col2 is the species_tissueID_experiment_barcode_{R1,R2}.fastq.gz # NOTE:
libids1=($wd/${spID}'libids1.txt') # this has the new trimmed ID gzipped fastq file paths
libids2=($wd/${spID}'libids.sh')

### 2. Mapping to genome - HiSat2
readalign=($wd/2.mapping) # assign mapping to genome dir
hisat_thread=32 # get HISAT2 to launch a specified number of parallel search threads
read1=($readalign/${spID}'read1.txt')
read2=($readalign/${spID}'read2.txt')
prefix=($readalign/${spID}'prefix.txt')
reads=($readalign/${spID}'reads.txt')
bam='.bam' # output BAM to add prefix beforehand

### 3. Alignment QC - Qualimap2
qualimap=($wd/3.qualimap)

### 4. Counts - HTseq
htseqdir=($wd/4.htseqcounts)

### 5. Differential gene expression - DESeq2
deseqdir=($WD/5.deseq)
prefixes=($deseqdir/prefixes.txt)
expdes=($deseqdir/expdesign.txt)




################################################################################################################

### 0. Copying the reads

## i. copy the fastq files over
# mkdir -p ${readsdirmain} # DONE
# cp ${readsgp}/*.fastq.gz $readsdirmain # DONE
# cp ${readsgp2}/*.fastq.gz $readsdirmain # DONE

################################################################################################################

### 1. Trimming using Trim Galore (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) and subsequently, an additional QC report is generated.

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
echo 'ml java' >> 1a.trimadaptors.sh
echo 'source fastqc-0.11.9' >> 1a.trimadaptors.sh
printf '\n' >> 1a.trimadaptors.sh
echo "for read1 in $readsdir/*_R1.fastq.gz; do read2="'$(echo $read1 | sed'" 's/_R1.fastq.gz/_R2.fastq.gz/'); "'echo $read1 >> trimarrayread1; echo $read2 >> trimarrayread2; done' >> 1a.trimadaptors.sh
printf '\n' >> 1a.trimadaptors.sh
echo 'mapfile -t read1 < trimarrayread1 # assign read1 as elements to $read1 variable' >> 1a.trimadaptors.sh
echo 'mapfile -t read2 < trimarrayread2 # assign read2 as elements to $read2 variable' >> 1a.trimadaptors.sh
printf '\n' >> 1a.trimadaptors.sh
echo "srun trim_galore --output_dir $trimdir --paired"' --fastqc ${read1[${SLURM_ARRAY_TASK_ID}]} ${read2[${SLURM_ARRAY_TASK_ID}]}' >> 1a.trimadaptors.sh

echo '# -- 1a. Adaptor trimming started -- #'

# assign '1a.trimadaptors.sh' to variable: JOBID1 and run
JOBID1=$( sbatch -W --array=$trimarray 1a.trimadaptors.sh | awk '{print $4}' ) # Run the first job and then store the first job to variable JOBID1 (taken by awk once run); Do not exit until the submitted job terminates.

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
echo "grep $spID $libids | sed 's/R1.fastq.gz/R1_val_1.fq.gz/g' | sed 's/R2.fastq.gz/R2_val_2.fq.gz/g' > $libids1 # grep the relevant species files from the long list" >> 1b.renamefiles.sh
echo "sed 's/^/mv /g' $libids1 > $libids2" >> 1b.renamefiles.sh
echo "sed -i '1 i\#!/bin/sh' $libids2" >> 1b.renamefiles.sh
printf '\n' >> 1b.renamefiles.sh
echo "sh $libids2" >> 1b.renamefiles.sh

echo '# -- 1a. Adaptor trimming completed -- #'

echo '# -- 1b. Renaming files started -- #'

# assign '1b.renamefiles.sh' to variable: JOBID2 and run
JOBID2=$( sbatch -W --dependency=afterok:${JOBID1} 1b.renamefiles.sh | awk '{print $4}' ) # JOB2 depends on JOB1 completing successfully

################################################################################################################

### 2. Mapping to genome - HiSat2
# Note: Final, name sorted BAM files will be created at step 4 and will be in the $htseqdir/*.sorted.bam

mkdir -p $readalign
cd $readalign

# 2a. build genome index for HiSat2 alignment
# NOTE: Since the pipeline will be ran several times, it will create genome indexes several times over - if running for one species then run here. Otherwise, run ATAC_Bioinf_pipeline_v2aA_gDNAindexes.sh BEFORE this pipeline.

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
echo 'source hisat2-2.2.1_CBG' >> 2a.genomeindex.sh
echo 'hisat2-build -f '$gFA $spG >> 2a.genomeindex.sh

echo '# -- 1b. Renaming files completed -- #'

echo '# -- 2a.'$spID' genome index building started -- #'

JOBID3=$( sbatch -W --dependency=afterok:${JOBID2} 2a.genomeindex.sh | awk '{print $4}' ) # JOB3 depends on JOB2 completing successfully

# 2b. genome alignment with HiSat2

echo '#!/bin/bash -e' > 2b.readalign.sh
echo '#SBATCH -p ei-largemem # partition (queue)' >> 2b.readalign.sh
echo '#SBATCH -N 1 # number of nodes' >> 2b.readalign.sh
echo '#SBATCH -c '$hisat_thread '# number of cores' >> 2b.readalign.sh
echo '#SBATCH --mem 512GB' >> 2b.readalign.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 2b.readalign.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 2b.readalign.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 2b.readalign.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 2b.readalign.sh
printf '\n' >> 2b.readalign.sh
echo "awk -F' ' '{print \$2}' $libids1 | sed -e"' "s|^|'$trimdir'/|g" > '"$reads" >> 2b.readalign.sh
echo 'mapfile -t reads < '$reads' # ${reads[0]} calls read1 AND ${reads[1]} calls read2' >> 2b.readalign.sh
echo "awk -F' ' '{print \$2}' " $libids1 " | awk -F'_' '{print \$1\"_\"\$2\"_\"\$3}' > "$prefix "# create a prefix file to iterate" >> 2b.readalign.sh
echo "mapfile -t prefixmap < $prefix"' # assign prefixes to $prefixmap' >> 2b.readalign.sh
printf '\n' >> 2b.readalign.sh
echo 'source hisat2-2.2.1_CBG' >> 2b.readalign.sh
echo 'ml samtools/1.7' >> 2b.readalign.sh
echo 'hisat2 -p '$hisat_thread' -x '$spG' -1 ${reads[0]} -2 ${reads[1]} --summary-file '${spID}'_'${spG}'.summary --met-file '${spID}'_'${spG}'.metrics | samtools view -Su /dev/stdin | samtools sort -o ${prefixmap[0]}'$bam >> 2b.readalign.sh

echo '# -- 2a.'$spID' genome index building completed -- #'

echo '# -- 2b.'$spID' read alignment started -- #'

JOBID4=$( sbatch -W --dependency=afterok:${JOBID3} 2b.readalign.sh | awk '{print $4}' ) # JOB4 depends on JOB2 completing successfully

################################################################################################################

### 3. Alignment QC - Qualimap2

mkdir -p $qualimap
cd $qualimap

echo '#!/bin/bash -e' > 3.qualimap.sh
echo '#SBATCH -p tgac-medium # partition (queue)' >> 3.qualimap.sh
echo '#SBATCH -N 1 # number of nodes' >> 3.qualimap.sh
echo '#SBATCH -n 1 # number of tasks' >> 3.qualimap.sh
echo '#SBATCH --mem 48000' >> 3.qualimap.sh
echo '#SBATCH -t 0-23:59' >> 3.qualimap.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 3.qualimap.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 3.qualimap.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 3.qualimap.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 3.qualimap.sh
printf '\n' >> 3.qualimap.sh
echo 'ml java' >> 3.qualimap.sh
echo 'source qualimap-2.2.1' >> 3.qualimap.sh
echo 'unset DISPLAY # remove display to make qualimap run' >> 3.qualimap.sh
echo 'zcat '$transcr' > "$(basename "'$transcr'" .gz)"' >> 3.qualimap.sh
echo 'bamfile=$(ls -1 '$readalign'/*.bam)' >> 3.qualimap.sh
printf '\n' >> 3.qualimap.sh
echo 'qualimap rnaseq \' >> 3.qualimap.sh
echo '-outdir qualimap \' >> 3.qualimap.sh
echo '-a proportional \' >> 3.qualimap.sh
echo '-bam ${bamfile} \' >> 3.qualimap.sh
echo '-p strand-specific-reverse \' >> 3.qualimap.sh
echo '-gtf "$(basename "'${transcr}'" .gz)" \' >> 3.qualimap.sh
echo '--java-mem-size=8G' >> 3.qualimap.sh
printf '\n' >> 3.qualimap.sh
echo 'rm "$(basename "'${transcr}'" .gz)"' >> 3.qualimap.sh


echo '# -- 2b.'$spID' read alignment completed -- #'

echo '# -- 3.'$spID' qualimap mapping stats started -- #'

JOBID5=$( sbatch -W --dependency=afterok:${JOBID4} 3.qualimap.sh | awk '{print $4}' ) # JOB4 depends on JOB3 completing successfully

################################################################################################################

### 4. Counts - HTseq
# Note: Final, name sorted BAM files will be created at this step and will be in the $htseqdir - *.sorted.bam

mkdir -p $htseqdir
cd $htseqdir

echo '#!/bin/bash -e' > 4.htseqcounts.sh
echo '#SBATCH -p tgac-medium # partition (queue)' >> 4.htseqcounts.sh
echo '#SBATCH -N 1 # number of nodes' >> 4.htseqcounts.sh
echo '#SBATCH -n 1 # number of tasks' >> 4.htseqcounts.sh
echo '#SBATCH --mem 32000' >> 4.htseqcounts.sh
echo '#SBATCH -t 0-16:59' >> 4.htseqcounts.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 4.htseqcounts.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 4.htseqcounts.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 4.htseqcounts.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 4.htseqcounts.sh
printf '\n' >> 4.htseqcounts.sh
echo 'source samtools-1.7' >> 4.htseqcounts.sh
echo 'source HTSeq-0.6.1' >> 4.htseqcounts.sh
echo 'bamfile=$(ls -1 '$readalign'/*.bam)' >> 4.htseqcounts.sh
echo 'zcat '$transcr' > "$(basename "'$transcr'" .gz)"' >> 4.htseqcounts.sh
printf '\n' >> 4.htseqcounts.sh
echo 'samtools sort -n -T "$(basename "$bamfile" .bam).sorted" -o "$(basename "$bamfile" .bam).sorted.bam" ${bamfile}' >> 4.htseqcounts.sh
echo 'samtools view -q 3 "$(basename "$bamfile" .bam).sorted.bam" | htseq-count --mode=intersection-nonempty --stranded=yes --type=exon --idattr=gene_id --quiet - "$(basename "'${transcr}'" .gz)" > '${spID}'_htseq_counts.out' >> 4.htseqcounts.sh
echo 'rm "$(basename "'${transcr}'" .gz)"' >> 4.htseqcounts.sh

echo '# -- 3.'$spID' qualimap mapping stats completed -- #'

echo '# -- 4. '$spID' HTseq counts started -- #'

JOBID6=$( sbatch -W --dependency=afterok:${JOBID5} 4.htseqcounts.sh | awk '{print $4}' ) # JOB5 depends on JOB4 completing successfully
# JOBID6=$( sbatch -W 4.htseqcounts.sh | awk '{print $4}' ) # JOB5 depends on JOB4 completing successfully


################################################################################################################

### 5. Differential gene expression - DESeq2
# Outputs: normalized quantification tables, some important statistics for the whole gene or transcript list, and the list of significantly differentially expressed genes (with default threshold of FDR<0.05)
# The raw count is normalized based on median-of-ratios method (if DESeq2 is used) when the reads are mapped to a genome.

mkdir -p $deseqdir
cd $deseqdir

## i. create an experimental design

awk -F' ' '{print $2}' $libids | awk -F'_' '{print $1"_"$2"_"$3}' | sort -u > $prefixes # create a prefix file to iterate

printf 'sampleID\tspecies\tgroup\tsubject\n' > ${expdes} # create four columns: sampleID, species, condition (tissue), replicate (individual)
awk '{print $1}' $prefixes | awk -F'_' '{print $0,$1,$3,$2}' OFS='\t' >> ${expdes}

## ii. R-script to run DEseq2

# directory <- "/path/to/your/files/"
#
# # We specify which files to read in using list.files, and select those files which contain the string "treated" using grep. The sub function is used to chop up the sample filename to obtain the condition status, or you might alternatively read in a phenotypic table using read.table.
#
# sampleFiles <- grep("treated",list.files(directory),value=TRUE)
# sampleCondition <- sub("(.*treated).*","\\1",sampleFiles)
# sampleTable <- data.frame(sampleName = sampleFiles,
#                           fileName = sampleFiles,
#                           condition = sampleCondition)
# sampleTable$condition <- factor(sampleTable$condition)
#
# # Then we build the DESeqDataSet using the following function:
#
# library("DESeq2")
#
# dds <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
#                                        directory = directory,
#                                        design= ~ condition)
# # dds$condition <- factor(dds$condition, levels = c("mock.2hr", "infected.2hr", "mock.4hr", "infected.4hr", "mock.6hr", "infected.6hr", "mock.8hr", "infected.8hr"))
#
# # get rid of low count data
# keep <- rowSums(counts(dds)) >= 10 dds <- dds[keep,]
#
# #perform the dfe analysis based on negative binomial model and get the results dds <- DESeq(dds) res <- results(dds)
#
# resultsNames(dds)
#
# # Use the contrasts function in DESeq2; Contrasts enable the user to generate results for all possible differences: log2 fold change of B vs A, of C vs A, and of C vs B.
# # The contrast argument of results function is used to extract test results of log2 fold changes of interest, for example:
#
# results(dds, contrast=c("condition","mock.2hr","infected.2hr"))


echo '# -- 4. '$spID' HTseq counts completed -- #'

echo '# -- 5. DESeq2 analysis started -- #'

# ################################################################################################################
