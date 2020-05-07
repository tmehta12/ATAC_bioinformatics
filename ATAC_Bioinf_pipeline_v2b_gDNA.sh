#!/bin/sh

#!/bin/bash -e
#SBATCH -p tgac-long # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 8000 # memory pool for all cores
#SBATCH -t 3-15:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

################################################################################################################

# ATAC-seq pipeline - Part 2 (gDNA only - to be ran prior to ATAC pipeline 'ATAC_Bioinf_pipeline_v2b.sh')
# March 2020: Tarang K. Mehta, Earlham Institute, Norwich, UK

################################################################################################################

# Script usage: ./ATAC_Bioinf_pipeline_v2b_gDNA.sh -s "spID" -g "spG" -f "gFA"
# e.g. ./ATAC_Bioinf_pipeline_v2b.sh -s Mz1_L_gDNA -g M_zebra_UMD1 -f /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa
# Note: Script is adapted for SBATCH usage

## Place this script and the following files in $WD of each gDNA folder only e.g. Mz1_L_gDNA (created in previous script)
# 1. libids.txt in $scripts folder: this has been prepared for previous script './ATAC_Bioinf_pipeline_v2a.sh': a 2-column space-delimited table where col1='R1/R2 filename's col2='desired species renamed filename: species_tissue_experiment e.g. Mz_L_ATAC/gDNA'
# 2. Run as an sbatch script with 8Gb memory and >3 day runtime - will spawn off other jobs
# 3. Run this gDNA script first for each experiment and then the ATAC once all gDNA are successfuly completed

################################################################################################################

# ~ This pipeline should be ran after 'ATAC_Bioinf_pipeline_v2a.sh' that merges files sequenced over mutilple lanes, but prior to 'ATAC_Bioinf_pipeline_v2b.sh' that works on the ATAC data for final peak calling

################################################################################################################

# ~ This pipeline is ran as species-specific, and contains the following components:

# 1. Trim adaptors - trimgalore
# 2. Read alignment - bowtie2 > samtools sorted bam

################################################################################################################

# Setting parameters for command line input

helpFunction()
{
   echo ""
   echo "Usage: $0 -s spID -g spG -f gFA"
   echo -e "\t-s spID = Species ID, preferably two short letters, individual and tissue e.g. Metriaclima zebra Individual 1 Liver ATAC/gDNA = Mz1_L_ATAC/Mz1_L_gDNA Note: this naming convention needs to be the same as renaming in space delimited file"
   echo -e "\t-g spG = Species genome ID e.g. hg19 or M_zebra_UMD1"
   echo -e "\t-f gFA = Full path to genome assembly in FASTA format"
   exit 1 # Exit script after printing help
}

while getopts "s:g:f:" opt
do
   case "$opt" in
      s ) spID="$OPTARG" ;;
      g ) spG="$OPTARG" ;;
      f ) gFA="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$spID" ] || [ -z "$spG" ] || [ -z "$gFA" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
echo "$spID"
echo "$spG"
echo "$gFA"

################################################################################################################

# All variables are added (and can be amended) here

scripts=(/tgac/workarea/group-vh/Tarang/ATACseq/2.run2) # place all scripts in the topmost directory - create this separately
WD=(/tgac/workarea/group-vh/Tarang/ATACseq/2.run2/$spID) # insert the working directory
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
RAarray=0-0 # INSERT the number range of paired *fastq.merged.gz to align in zero base e.g. 10 pairs = 0-9
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
echo "srun trim_galore --output_dir $trimdir --paired"' --fastqc ${read1[${SLURM_ARRAY_TASK_ID}]} ${read2[${SLURM_ARRAY_TASK_ID}]}' >> 1a.trimadaptors.sh

echo '# -- 1a. Adaptor trimming started -- #'

# assign '1a.trimadaptors.sh' to variable: JOBID1 and run
JOBID1=$( sbatch -W --array=$trimarray 1a.trimadaptors.sh | awk '{print $4}' ) # Run the first job and then store the first job to variable JOBID1 (taken by awk once run); Do not exit until the submitted job terminates.

# rename the files according to species, tissue and experiment - provide this as a 2-column SPACE-delimited table that will be assigned as a variable above, placed in the trimdir
# only provide for the two paired files you are working on and not all files

# Example is (note, these are made up examples)
# echo 'PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R1.fastq.merged.gz Pn1_T_gDNA_TCGCCTGC-AACCGCCA_L001_R1.fastq.merged.gz' > libids.txt
# echo 'PRO1563_S1_lib_TCGCCTGC-AACCGCCA_L001_R2.fastq.merged.gz Pn1_T_gDNA_TCGCCTGC-AACCGCCA_L001_R2.fastq.merged.gz' >> libids.txt

# The trimmed files need calling to rename, not the original files!
PRO1563_S1_lib_CAGAATGC-TAACTCTA_L001_R1.fastq.merged.gz_trimmed.fq.gz
PRO1563_S1_lib_CAGAATGC-TAACTCTA_L001_R2.fastq.merged.gz_trimmed.fq.gz


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
echo "grep $spID $libids | sed 's/.fastq.merged.gz/.fastq.merged.gz_trimmed.fq.gz/g' > $libids1 # grep the relevant species files from the long list" >> 1b.renamefiles.sh
echo "sed 's/^/mv /g' $libids1 > $libids2" >> 1b.renamefiles.sh
echo "sed -i '1 i\\n' $libids2" >> 1b.renamefiles.sh
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
echo 'mapfile -t reads < '$reads' # ${reads[0]} calls read1 AND ${reads[1]} calls read2' >> 2b.readalign.sh
echo "awk -F' ' '{print \$2}' " $libids1 " | awk -F'_' '{print \$1\"_\"\$2\"_\"\$3}' > "$prefix "# create a prefix file to iterate" >> 2b.readalign.sh
echo 'mapfile -t prefixmap < '$prefix '# assign prefixes to $prefixmap' >> 2b.readalign.sh
echo '# run bowtie2 with multimapping and threading, then output sorted BAM file' >> 2b.readalign.sh
echo 'srun bowtie2 -k ' $multimapping ' -X2000 --mm --threads ' $bwt_thread ' -x ' $idx ' -1 ${reads[0]} -2 ${reads[1]} 2>'$prefixmap$log '| samtools view -Su /dev/stdin | samtools sort -o $prefixmap '$bam >> 2b.readalign.sh
printf '\n' >> 2b.readalign.sh
echo 'samtools flagstat ' $prefixmap$bam ' > ' $prefixmap$fgQC1 ' # output alignment stats' >> 2b.readalign.sh

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

# JOBID4=$( sbatch --dependency=afterok:${JOBID3} --array=$RAarray 2b.readalign.sh | awk '{print $4}' ) # JOB4 depends on JOB3 completing successfully

echo '# -- 2a.'$spID' genome index building completed -- #'

echo '# -- 2b.'$spID' read alignment started -- #'

JOBID4=$( sbatch -W --dependency=afterok:${JOBID3} 2b.readalign.sh | awk '{print $4}' ) # JOB4 depends on JOB3 completing successfully

################################################################################################################

### 3. Genomic DNA library alignments completion

echo '#!/bin/bash -e' > 3.gDNA_alignments_complete.sh
echo '#SBATCH -p tgac-short # partition (queue)' >> 3.gDNA_alignments_complete.sh
echo '#SBATCH -N 1 # number of nodes' >> 3.gDNA_alignments_complete.sh
echo '#SBATCH -n 1 # number of tasks' >> 3.gDNA_alignments_complete.sh
echo '#SBATCH --mem 2000' >> 3.gDNA_alignments_complete.sh
echo '#SBATCH -t 0-00:45' >> 3.gDNA_alignments_complete.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> 3.gDNA_alignments_complete.sh
echo "#SBATCH --mail-user=$email # send-to address" >> 3.gDNA_alignments_complete.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> 3.gDNA_alignments_complete.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> 3.gDNA_alignments_complete.sh
printf '\n' >> 3.gDNA_alignments_complete.sh
echo "echo 'Genomic DNA library alignment is complete' > gDNAcompleted.txt" >> 3.gDNA_alignments_complete.sh

echo '# -- 2b.'$spID' read alignment completed -- #'

JOBID5=$( sbatch -W --dependency=afterok:${JOBID4} 3.gDNA_alignments_complete.sh | awk '{print $4}' ) # JOB5 depends on JOB4 completing successfully

################################################################################################################

### Finish the script
exit 0
