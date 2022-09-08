#!/bin/sh

##############################################################################
### July 2021 - Running the RNAseq pipeline
##############################################################################

# ~ This pipeline will be ran as species-specific, and contains the following components:

# 1. Trimming using Trim Galore (https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) and subsequently, an additional QC report is generated.
# 2. Mapping to genome (can also map to transcriptome but opting for less bias this way) - HiSat2
# 3. Alignment QC - Qualimap2
# 4. Counts - featureCounts or HTseq
# 5. Differential gene expression - edgeR or DESeq2
  # Outputs: normalized quantification tables, some important statistics for the whole gene or transcript list, and the list of significantly differentially expressed genes or transcripts (with default threshold of FDR<0.05)
  # The raw count is normalized based on Trimmed Mean of M values (TMM) (if edgeR is used) or the median-of-ratios method (if DESeq2 is used) when the reads are mapped to a genome.

##############################################################################

# Variables

topdir=(/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch) # place all scripts in the topmost directory - create this separately
wd=(${topdir}/RNAseq) # insert the working directory
email=Tarang.Mehta@earlham.ac.uk # SBATCH out and err send to address

readsdir=(${wd}/0.reads)

libids=($wd/libids.txt) # 2-col space-delimited file where col1 is the *_{R1,R2}.fastq.gz and col2 is the species_tissueID_experiment_barcode_{R1,R2}.fastq.gz # NOTE:
libids2=($wd/libids2.txt)
libids3=($wd/libids3.txt) # 2-col space-delimited file where col1 is the *_{R1,R2}.fastq.merged.gz and col2 is the species_tissue_experiment only e.g. Pnm1_L_ATAC
prefix=($wd/prefix.txt)
prefix2=($wd/prefix2.txt)
prefix3=($wd/prefix3.txt)

RNAscr=(RNAseq_Bioinf_pipeline_v1a.sh)

# Define species prefix here according to filenames etc.
sp1=Mz
sp2=Pn
sp3=Ab
sp4=Nb
sp5=On

# Define species genome ID here
sp1ID=(M_zebra_UMD2a)
sp2ID=(PunNye1.0)
sp3ID=(AstBur1.0)
sp4ID=(NeoBri1.0)
sp5ID=(O_niloticus_UMD_NMBU)

user=mehtat # userID for cluster to download mitochondrial genomes on internet access node

## Genomes and Annotations
genomesdir=(/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/data/ATACseq/3.run2/genomes)
# Metriaclima zebra - M_zebra_UMD2a (PacBio)
mzeb=($genomesdir/Mzebra)
mzebgenomelg=$mzeb/dna/Maylandia_zebra.M_zebra_UMD2a.dna.primary_assembly.allLG.fa
mzebgenomenc=$mzeb/dna/Maylandia_zebra.M_zebra_UMD2a.dna.nonchromosomal.fa
mzebgenome=$mzeb/dna/Maylandia_zebra.M_zebra_UMD2a.dna.primary_assembly.allLG.and.nonchromosomal.fa
mzebgtf=$mzeb/current_gtf/maylandia_zebra/Maylandia_zebra.M_zebra_UMD2a.101.gtf.gz

# Pundamilia nyererei
pnye=($genomesdir/Pnyererei)
pnyegenome=$pnye/dna/Pundamilia_nyererei.PunNye1.0.dna.nonchromosomal.fa
pnyegtf=$pnye/current_gtf/pundamilia_nyererei/Pundamilia_nyererei.PunNye1.0.101.gtf.gz

# Astatotilapia burtoni
abur=($genomesdir/Aburtoni)
aburgenome=$abur/dna/Haplochromis_burtoni.AstBur1.0.dna.nonchromosomal.fa
aburgtf=$abur/current_gtf/haplochromis_burtoni/Haplochromis_burtoni.AstBur1.0.101.gtf.gz

# Neolamprologus brichardi
nbri=($genomesdir/Nbrichardi)
nbrigenome=$nbri/dna/Neolamprologus_brichardi.NeoBri1.0.dna.nonchromosomal.fa
nbrigtf=$nbri/current_gtf/neolamprologus_brichardi/Neolamprologus_brichardi.NeoBri1.0.101.gtf.gz

# Oreochromis niloticus - O_niloticus_UMD_NMBU (PacBio)
onil=($genomesdir/Oniloticus)
onilgenomelg=$onil/dna/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.primary_assembly.allLG.fa
onilgenomenc=$onil/dna/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.nonchromosomal.fa
onilgenome=$onil/dna/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.primary_assembly.allLG.and.nonchromosomal.fa
onilgtf=$onil/current_gtf/oreochromis_niloticus/Oreochromis_niloticus.O_niloticus_UMD_NMBU.101.gtf.gz

##############################################################################

## 0. Create working directory

mkdir -p $wd # DONE

##############################################################################

## 1. Create the appropriate directory structure and file paths for downstream analysis (and to run following scripts)

# 1a. Create appropriate directory structure - this uses the $libids file

# Example of $libids is (note, these are made up examples)
# R0635-S0001_9Pnm1B9_A56807_1_HF53CDSX2_CGCAACTA-GAATCCGA_L004_R1.fastq.gz	Pn_1_Brain_CGCAACTA-GAATCCGA_R1.fastq.gz
# R0635-S0001_9Pnm1B9_A56807_1_HF53CDSX2_CGCAACTA-GAATCCGA_L004_R2.fastq.gz	Pn_1_Brain_CGCAACTA-GAATCCGA_R2.fastq.gz

# A. read in the space delimited file and prepare working directories for all col2 entries
# B. within each species_tissue_experiment working directory, create a raw reads directory too

cd $wd

awk -F' ' '{print $2}' $libids | awk -F'_' '{print $1"_"$2"_"$3}' | sort -u > $prefix # create a prefix file to iterate

while IFS= read -r i; do
  # echo $i
  mkdir $wd/$i
  cd $i
  mkdir 0.rawreads
  cd ../
done < $prefix

# 1b. Create symbolic links to raw reads in per species_tissue - no file renaming required at this stage; will be done after trimming

cd $wd

awk -F' ' '{print $2}' $libids | awk -F'_' '{print $1"_"$2"_"$3}' > $prefix3
cut -f1 $libids > $libids2
paste -d' ' $libids2 $prefix3 > $libids3 # 2-col space-delimited file where col1 is the *_{R1,R2}.fastq.merged.gz and col2 is the species_tissue_experiment only e.g. Pnm1_L_ATAC

# while read, if statement to create symbolic link to raw reads for each species_tissue_experiment

while IFS= read -r i; do
  # echo $i
  spdir=$(echo $i | awk -F' ' '{print $2}') # pull out the species_tissue_experiment ID of each line
  # echo $spdir
  cd $wd/$spdir # change to each subdirectory
  dir=$(pwd) # create variable of the full path
  # echo $dir
  dir2=$(echo $dir | sed "s|$wd/||g") # pull out the species_tissue_experiment ID of the working path
  # echo $dir2
  read=$(echo $i | awk -F' ' '{print $1}') # pull out the corresponding read filename of each line
  # echo $read
  if [[ $spdir == $dir2 ]]; then # if the species_tissue_experiment ID in the list and the working directory match, then
    ln -s $readsdir/$read 0.rawreads/ # create symbolic link in raw reads dir to corresponding raw read
  fi
done < $libids3


################################################################################################################

## 2. Add relevant scripts to each directory

# Add each RNAseq_Bioinf_pipeline_v1a.sh script to each species folder

cd $wd

while IFS= read -r i; do
  cp $RNAscr $i
done < $prefix

##############################################################################

## 3. trimming, mapping, QC and generating counts of RNA-seq reads

# # Example
# cd $wd/Mz_1_Liver
# sbatch RNAseq_Bioinf_pipeline_v1a.sh -s Mz_1_Liver -g M_zebra_UMD1 -f $mzebgenome -u mehtat -a $mzebgtf

# A. Create a prefix file to iterate for running scripts
awk -F'_' '{print $1,$1"_"$2"_"$3}' OFS='\t' $prefix > $prefix2

# B. While loop with if statements to appropriately spawn off sbatch scripts
while read -r i1 i2; do
  if [[ "$i1" == "$sp1" ]]; then
    # echo "New if statement: "$i1 $i2 $sp1
    cd $wd/$i2
    # echo "sbatch $RNAscr -s $i2 -g $sp1ID -f $mzebgenome -m $sp1mitID -u $user -a $mzebgtf"
    sbatch $RNAscr -s $i2 -g $sp1ID -f $mzebgenome -u $user -t $mzebgtf
  fi
  if [[ "$i1" == "$sp2" ]]; then
    # echo "New if statement: "$i1 $i2 $sp2
    cd $wd/$i2
    # echo "sbatch $RNAscr -s $i2 -g $sp2ID -f $pnyegenome -m $sp2mitID -u $user -a $pnyegtf"
    sbatch $RNAscr -s $i2 -g $sp2ID -f $pnyegenome -u $user -t $pnyegtf
  fi
  if [[ "$i1" == "$sp3" ]]; then
    # echo "New if statement: "$i1 $i2 $sp3
    cd $wd/$i2
    # echo "sbatch $RNAscr -s $i2 -g $sp3ID -f $aburgenome -m $sp3mitID -u $user -a $aburgtf"
    sbatch $RNAscr -s $i2 -g $sp3ID -f $aburgenome -u $user -t $aburgtf
  fi
  if [[ "$i1" == "$sp4" ]]; then
    # echo "New if statement: "$i1 $i2 $sp4
    cd $wd/$i2
    # echo "sbatch $RNAscr -s $i2 -g $sp4ID -f $nbrigenome -m $sp4mitID -u $user -a $nbrigtf"
    sbatch $RNAscr -s $i2 -g $sp4ID -f $nbrigenome -u $user -t $nbrigtf
  fi
  if [[ "$i1" == "$sp5" ]]; then
    # echo "New if statement: "$i1 $i2 $sp5
    cd $wd/$i2
    # echo "sbatch $RNAscr -s $i2 -g $sp5ID -f $onilgenome -m $sp5mitID -u $user -a $onilgtf"
    sbatch $RNAscr -s $i2 -g $sp5ID -f $onilgenome -u $user -t $onilgtf
  fi
done < $prefix2 # this while loop will run RNA script in each RNA folder

## cancel failing jobs
squeue -u mehtat | grep 36975 | awk '{print $1}' | xargs -n 1 scancel

## checking jobs and the samples associated (change the grep accordingly)
cd $wd
jobids=$(squeue -u mehtat | grep '2a.genom' | awk -F' ' '{print $1}' | awk 'BEGIN { ORS = " " } { print }') # this will store the jobids into a single row variable to use in for loop

for i in $jobids; do
  com=$(scontrol show jobid -dd $i | grep 'Command=' )
  err=$(scontrol show jobid -dd $i | grep 'StdErr=')
  echo -e $com'\t'$err
done # this will output two columns - first col shows the command where error occured and second shows the stout error file associated. The folder ran before the command will invariably have the error with it.

################################################################################################################

## 4. Mapping to genome - checking alignment rates from native HiSat2 output:

for i in /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/RNAseq/*_*_{Brain,Eye,Liver,Gill,Testis}/2.mapping/*.summary; do
  aln=$(grep 'overall alignment rate' ${i})
  file="$(basename "${i}" .summary)"
  echo -e ${file}'\t'${aln}
done

# Ab_5_Brain_AstBur1.0	82.04% overall alignment rate
# Mz_1_Brain_M_zebra_UMD2a	96.64% overall alignment rate
# Mz_2_Brain_M_zebra_UMD2a	96.46% overall alignment rate
# Nb_4_Brain_NeoBri1.0	90.40% overall alignment rate
# Nb_5_Brain_NeoBri1.0	66.73% overall alignment rate
# On_1_Brain_O_niloticus_UMD_NMBU	90.50% overall alignment rate
# On_2_Brain_O_niloticus_UMD_NMBU	90.73% overall alignment rate
# On_3_Brain_O_niloticus_UMD_NMBU	92.59% overall alignment rate
# Pn_1_Brain_PunNye1.0	86.94% overall alignment rate
# Pn_2_Brain_PunNye1.0	85.46% overall alignment rate
# Pn_3_Brain_PunNye1.0	86.61% overall alignment rate
# Pn_4_Brain_PunNye1.0	86.88% overall alignment rate
# Mz_1_Eye_M_zebra_UMD2a	96.21% overall alignment rate
# Mz_2_Eye_M_zebra_UMD2a	97.16% overall alignment rate
# On_1_Eye_O_niloticus_UMD_NMBU	88.16% overall alignment rate
# On_2_Eye_O_niloticus_UMD_NMBU	90.57% overall alignment rate
# On_3_Eye_O_niloticus_UMD_NMBU	91.09% overall alignment rate
# Pn_2_Eye_PunNye1.0	83.69% overall alignment rate
# Pn_3_Eye_PunNye1.0	83.00% overall alignment rate
# Pn_4_Eye_PunNye1.0	79.65% overall alignment rate
# Ab_5_Liver_AstBur1.0	83.11% overall alignment rate
# Ab_6_Liver_AstBur1.0	82.79% overall alignment rate
# Mz_1_Liver_M_zebra_UMD2a	96.42% overall alignment rate
# Mz_2_Liver_M_zebra_UMD2a	96.15% overall alignment rate
# Nb_4_Liver_NeoBri1.0	73.01% overall alignment rate
# Nb_5_Liver_NeoBri1.0	71.49% overall alignment rate
# On_1_Liver_O_niloticus_UMD_NMBU	93.95% overall alignment rate
# On_2_Liver_O_niloticus_UMD_NMBU	93.07% overall alignment rate
# On_3_Liver_O_niloticus_UMD_NMBU	92.67% overall alignment rate
# Pn_1_Liver_PunNye1.0	82.10% overall alignment rate
# Pn_2_Liver_PunNye1.0	85.37% overall alignment rate
# Pn_3_Liver_PunNye1.0	86.33% overall alignment rate
# Pn_4_Liver_PunNye1.0	82.49% overall alignment rate
# On_1_Gill_O_niloticus_UMD_NMBU	92.87% overall alignment rate
# On_2_Gill_O_niloticus_UMD_NMBU	91.12% overall alignment rate
# Ab_5_Testis_AstBur1.0	86.56% overall alignment rate
# Ab_6_Testis_AstBur1.0	88.04% overall alignment rate
# Mz_1_Testis_M_zebra_UMD2a	96.82% overall alignment rate
# Mz_2_Testis_M_zebra_UMD2a	96.96% overall alignment rate
# On_2_Testis_O_niloticus_UMD_NMBU	95.15% overall alignment rate
# On_3_Testis_O_niloticus_UMD_NMBU	90.02% overall alignment rate
# Pn_1_Testis_PunNye1.0	88.10% overall alignment rate
# Pn_2_Testis_PunNye1.0	87.78% overall alignment rate
# Pn_3_Testis_PunNye1.0	89.57% overall alignment rate
# Pn_4_Testis_PunNye1.0	84.92% overall alignment rate

## check all mapping ran - we expect 45 outfiles
ls -1 /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/RNAseq/*_*_{Brain,Eye,Liver,Gill,Testis}/2.mapping/*.summary | wc -l # DONE - all good

################################################################################################################

## 5. Mapping to genome - check the number of missing mates encountered when running HTseq counts:

cd /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/RNAseq

for i in *_*_*/4.htseqcounts/slurm*.err; do
  mm=$(grep 'missing mate' ${i} | sed 's|Warning: ||g' | sed 's| encountered.||g')
  echo -e ${i}'\t'${mm}
done

# Ab_5_Brain/4.htseqcounts/slurm.t768n1.37039424.err	41151 reads with missing mate
# Ab_5_Liver/4.htseqcounts/slurm.t128n68.37039432.err	47013 reads with missing mate
# Ab_5_Testis/4.htseqcounts/slurm.t128n55.37039430.err	37229 reads with missing mate
# Ab_6_Liver/4.htseqcounts/slurm.t128n68.37039434.err	41040 reads with missing mate
# Ab_6_Testis/4.htseqcounts/slurm.t128n68.37039436.err	32122 reads with missing mate
# Mz_1_Brain/4.htseqcounts/slurm.t128n24.37039575.err	91290 reads with missing mate
# Mz_1_Eye/4.htseqcounts/slurm.t256n14.37039611.err	118564 reads with missing mate
# Mz_1_Liver/4.htseqcounts/slurm.t128n52.37039577.err	124664 reads with missing mate
# Mz_1_Testis/4.htseqcounts/slurm.t128n79.37039572.err	129067 reads with missing mate
# Mz_2_Brain/4.htseqcounts/slurm.t256n11.37039615.err	142000 reads with missing mate
# Mz_2_Eye/4.htseqcounts/slurm.t128n62.37039582.err	77983 reads with missing mate
# Mz_2_Liver/4.htseqcounts/slurm.t256n14.37039614.err	164838 reads with missing mate
# Mz_2_Testis/4.htseqcounts/slurm.t128n62.37039580.err	117146 reads with missing mate
# Nb_4_Brain/4.htseqcounts/slurm.t128n52.37039578.err	54152 reads with missing mate
# Nb_4_Liver/4.htseqcounts/slurm.t256n6.37039610.err	63025 reads with missing mate
# Nb_5_Brain/4.htseqcounts/slurm.t256n14.37039612.err	23463 reads with missing mate
# Nb_5_Liver/4.htseqcounts/slurm.t128n79.37039573.err	24592 reads with missing mate
# On_1_Brain/4.htseqcounts/slurm.t128n24.37039574.err	68406 reads with missing mate
# On_1_Eye/4.htseqcounts/slurm.t128n24.37039576.err	40481 reads with missing mate
# On_1_Gill/4.htseqcounts/slurm.t128n52.37039579.err	102522 reads with missing mate
# On_1_Liver/4.htseqcounts/slurm.t256n14.37039613.err	171688 reads with missing mate
# On_2_Brain/4.htseqcounts/slurm.t128n62.37039581.err	77113 reads with missing mate
# On_2_Eye/4.htseqcounts/slurm.t256n11.37039616.err	52426 reads with missing mate
# On_2_Gill/4.htseqcounts/slurm.t256n12.37040035.err	124987 reads with missing mate
# On_2_Liver/4.htseqcounts/slurm.t256n11.37040031.err	118782 reads with missing mate
# On_2_Testis/4.htseqcounts/slurm.t256n12.37040033.err	67810 reads with missing mate
# On_3_Brain/4.htseqcounts/slurm.t256n11.37040032.err	58927 reads with missing mate
# On_3_Eye/4.htseqcounts/slurm.t128n81.37040040.err	86436 reads with missing mate
# On_3_Liver/4.htseqcounts/slurm.t256n12.37040036.err	118161 reads with missing mate
# On_3_Testis/4.htseqcounts/slurm.t128n81.37040039.err	138001 reads with missing mate
# Pn_1_Brain/4.htseqcounts/slurm.t256n12.37040034.err	34805 reads with missing mate
# Pn_1_Liver/4.htseqcounts/slurm.t256n4.37040037.err	63960 reads with missing mate
# Pn_1_Testis/4.htseqcounts/slurm.t128n81.37040038.err	38297 reads with missing mate
# Pn_2_Brain/4.htseqcounts/slurm.t128n51.37040046.err	43169 reads with missing mate
# Pn_2_Eye/4.htseqcounts/slurm.t128n46.37040048.err	40457 reads with missing mate
# Pn_2_Liver/4.htseqcounts/slurm.t128n16.37040051.err	54486 reads with missing mate
# Pn_2_Testis/4.htseqcounts/slurm.t128n16.37040050.err	37956 reads with missing mate
# Pn_3_Brain/4.htseqcounts/slurm.t128n51.37040045.err	42886 reads with missing mate
# Pn_3_Eye/4.htseqcounts/slurm.t128n65.37040041.err	31224 reads with missing mate
# Pn_3_Liver/4.htseqcounts/slurm.t256n6.37052117.err	70800 reads with missing mate
# Pn_3_Testis/4.htseqcounts/slurm.t128n65.37040042.err	38482 reads with missing mate
# Pn_4_Brain/4.htseqcounts/slurm.t128n65.37040043.err	61445 reads with missing mate
# Pn_4_Eye/4.htseqcounts/slurm.t128n51.37040044.err	41334 reads with missing mate
# Pn_4_Liver/4.htseqcounts/slurm.t128n16.37040052.err	69158 reads with missing mate
# Pn_4_Testis/4.htseqcounts/slurm.t128n46.37040049.err	55529 reads with missing mate

################################################################################################################
