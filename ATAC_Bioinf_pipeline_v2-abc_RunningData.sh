#!/bin/sh

##############################################################################
### Aug 2020 - Running the ATAC/gDNA data using pipeline
## The pipeline is run with the following scripts in this order:
# 1. ATAC_Bioinf_pipeline_v2a.sh
# 2. ATAC_Bioinf_pipeline_v2b_gDNA.sh
# 3. ATAC_Bioinf_pipeline_v2b.sh (scripts below called within)
  # 3a. ATAC_Bioinf_pipeline_v2b_part3a.py
  # 3b. ATAC_Bioinf_pipeline_v2b_part3b.R
# 4. ATAC_Bioinf_pipeline_v2c.sh (script below called within)
  # 4a. ATAC_Bioinf_pipeline_v2c_part3a.R
##############################################################################

### Variables here:

WD=(/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/data/ATACseq/3.run2) # if you change path here then ensure to change in all other scripts
WD2=(/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2)
scripts=(/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/data/ATACseq/3.run2)
rawreadfolder3=(/ei/data/reads/PIP-2582/200717_A00478_0125_BHTF7MDMXX)
libids=($scripts/libids.txt) # 2-col space-delimited file where col1 is the *_{R1,R2}.fastq.merged.gz and col2 is the species_tissue_experiment_barcode_{R1,R2}.fastq.merged.gz e.g. Mz_L_ATAC/gDNA
mergefiles1=($scripts/mergefiles_1.txt) # create a tab-delimited file where col1 is the sampleID e.g. Pnm1_L_gDNA, and all proceeding columns of each line are all the files that you want to merge e.g. replicates sequenced on different lanes. This will be automatically merged below so column numbers for each row can be different.
mergefiles2=($scripts/mergefiles_2.txt) # this file removes the sample ID - each line is all the files that you want to merge e.g. replicates sequenced on different lanes. This will be automatically merged below so column numbers for each row can be different.

genomesdir=($WD/genomes)

ATACscript1=(ATAC_Bioinf_pipeline_v2a.sh)
ATACscript2=(ATAC_Bioinf_pipeline_v2b_gDNA.sh)
ATACscript3=(ATAC_Bioinf_pipeline_v2b.sh)
ATACscript4=(ATAC_Bioinf_pipeline_v2c.sh)

prefix=($scripts/prefix.txt)
prefixgDNA=($scripts/prefixgDNA.txt)
prefixgDNA2=($scripts/prefixgDNA2.txt)
prefixATAC=($scripts/prefixATAC.txt)
prefixATAC2=($scripts/prefixATAC2.txt)

# Define species prefix here according to filenames etc.
sp1=Mz
sp2=Pn
sp3=Ab
sp4=Nb
sp5=On
sp6=Ac

# Define species genome ID here
sp1ID=(M_zebra_UMD2a)
sp2ID=(PunNye1.0)
sp3ID=(AstBur1.0)
sp4ID=(NeoBri1.0)
sp5ID=(O_niloticus_UMD_NMBU)
sp6ID=(fAstCal1.2)

# Define each species NCBI mitochondrial accession numbers here
sp1mitID=(NC_027944.1)
sp2mitID=(NC_028011.1)
sp3mitID=(NC_027289.1)
sp4mitID=(NC_009062.1)
sp5mitID=(NC_013663.1)
sp6mitID=(NC_018560.1)

user=mehtat # userID for cluster to download mitochondrial genomes on internet access node

## Genomes and Annotations
# Metriaclima zebra - M_zebra_UMD2a (PacBio)
# mzeb=(/tgac/workarea/group-vh/Tarang/Reference_Genomes/ensembl/cichlids/Mzebra)
# mzebgenome=$mzeb/dna/Maylandia_zebra.M_zebra_UMD2a.dna.primary_assembly.allLG.fa
# mzebgtf=$mzeb/current_gtf/maylandia_zebra/Maylandia_zebra.M_zebra_UMD2a.100.gtf.gz
mzeb=($genomesdir/Mzebra)
mzebgenome=$mzeb/dna/Maylandia_zebra.M_zebra_UMD2a.dna.primary_assembly.allLG.fa
mzebgtf=$mzeb/current_gtf/maylandia_zebra/Maylandia_zebra.M_zebra_UMD2a.101.gtf.gz

# Pundamilia nyererei
# pnye=(/tgac/workarea/group-vh/Tarang/Reference_Genomes/ensembl/cichlids/Pnyererei)
# pnyegenome=$pnye/dna/Pundamilia_nyererei.PunNye1.0.dna.nonchromosomal.fa
# pnyegtf=$pnye/current_gtf/pundamilia_nyererei/Pundamilia_nyererei.PunNye1.0.100.gtf.gz
pnye=($genomesdir/Pnyererei)
pnyegenome=$pnye/dna/Pundamilia_nyererei.PunNye1.0.dna.nonchromosomal.fa
pnyegtf=$pnye/current_gtf/pundamilia_nyererei/Pundamilia_nyererei.PunNye1.0.101.gtf.gz

# Astatotilapia burtoni
# abur=(/tgac/workarea/group-vh/Tarang/Reference_Genomes/ensembl/cichlids/Aburtoni)
# aburgenome=$abur/dna/Haplochromis_burtoni.AstBur1.0.dna.nonchromosomal.fa
# aburgtf=$abur/current_gtf/haplochromis_burtoni/Haplochromis_burtoni.AstBur1.0.100.gtf.gz
abur=($genomesdir/Aburtoni)
aburgenome=$abur/dna/Haplochromis_burtoni.AstBur1.0.dna.nonchromosomal.fa
aburgtf=$abur/current_gtf/haplochromis_burtoni/Haplochromis_burtoni.AstBur1.0.101.gtf.gz

# Neolamprologus brichardi
# nbri=(/tgac/workarea/group-vh/Tarang/Reference_Genomes/ensembl/cichlids/Nbrichardi)
# nbrigenome=$nbri/dna/Neolamprologus_brichardi.NeoBri1.0.dna.nonchromosomal.fa
# nbrigtf=$nbri/current_gtf/neolamprologus_brichardi/Neolamprologus_brichardi.NeoBri1.0.100.gtf.gz
nbri=($genomesdir/Nbrichardi)
nbrigenome=$nbri/dna/Neolamprologus_brichardi.NeoBri1.0.dna.nonchromosomal.fa
nbrigtf=$nbri/current_gtf/neolamprologus_brichardi/Neolamprologus_brichardi.NeoBri1.0.101.gtf.gz

# Oreochromis niloticus - O_niloticus_UMD_NMBU (PacBio)
# onil=(/tgac/workarea/group-vh/Tarang/Reference_Genomes/ensembl/cichlids/Oniloticus)
# onilgenome=$onil/dna/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.primary_assembly.allLG.fa
# onilgtf=$onil/current_gtf/oreochromis_niloticus/Oreochromis_niloticus.O_niloticus_UMD_NMBU.100.gtf.gz
onil=($genomesdir/Oniloticus)
onilgenome=$onil/dna/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.primary_assembly.allLG.fa
onilgtf=$onil/current_gtf/oreochromis_niloticus/Oreochromis_niloticus.O_niloticus_UMD_NMBU.101.gtf.gz

# Astatotilapia calliptera
# acal=(/tgac/workarea/group-vh/Tarang/Reference_Genomes/ensembl/cichlids/Acalliptera)
# acalgenome=$acal/dna/Astatotilapia_calliptera.fAstCal1.2.dna.nonchromosomal.fa
# acalgtf=$acal/current_gtf/astatotilapia_calliptera/Astatotilapia_calliptera.fAstCal1.2.100.gtf.gz
acal=($genomesdir/Acalliptera)
acalgenome=$acal/dna/Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.fa
acalgtf=$acal/current_gtf/astatotilapia_calliptera/Astatotilapia_calliptera.fAstCal1.2.101.gtf.gz

##############################################################################

## 0. Create working directory
mkdir -p $WD
cd $WD

##############################################################################

## 1. Merging files and create appropriate data structure - ATAC_Bioinf_pipeline_v2a.sh

# a. Creating the $mergefiles file in command line/excel and then port over (change line endings!)
ls -1 $rawreadfolder3/*.fastq.gz | awk -F'_' '{print $5"_"$6}' | sed 's/A14791_//g' | sed 's/_.*//g' | sed 's/63ATACNb4L64ATACNb4T/63ATACNb4L/g' > rawreadfolder3_files.tmp1
ls -1 $rawreadfolder3/*.fastq.gz > rawreadfolder3_files.tmp2
paste rawreadfolder3_files.tmp1 rawreadfolder3_files.tmp2 > rawreadfolder3_files # then manually duplicated some index hopped files using nano
rm rawreadfolder3_files.tmp*
# b. For files that could be 'mixed' samples due to index hopping - either create duplicates of the file but with different sample prefixes e.g. 63ATACNb4L_On2T and On2T_63ATACNb4L (but this will not work with the pipeline)
# Instead, going to make the assumption that since no On2T indvidual reads were created, that 'R0267-S0018_63ATACNb4LOn2T_A14708_HTF7MDMXX_GTCAACCG-TCCTTCTT_L001_{R1,R2}.fastq.gz' is actually On2T so rename accordingly

# c. Create the $libids file using $mergefiles1 - this is specific for your files so do change
awk '{print $2}' $mergefiles1 | sed 's/\/tgac\/data\/reads\/TarangMehta_EI_EI_TM_ENQ-1771_A_03\/170214_D00507_0252_AHG5KGBCXY\///g' | sed 's/\/ei\/data\/reads\/PIP-2582\/200717_A00478_0125_BHTF7MDMXX\///g' | sed 's/fastq.gz/fastq.merged.gz/g' > $libids.tmp1
awk '{print $1"_"$2}' $mergefiles1 | sed 's/\/tgac\/data\/reads\/TarangMehta_EI_EI_TM_ENQ-1771_A_03\/170214_D00507_0252_AHG5KGBCXY\///g' | sed 's/\/ei\/data\/reads\/PIP-2582\/200717_A00478_0125_BHTF7MDMXX\///g' | sed 's/fastq.gz/fastq.merged.gz/g' > $libids.tmp2
paste $libids.tmp1 $libids.tmp2 > $libids
rm $libids.tmp*

sbatch $ATACscript1 # {DONE - 24/08/2020}

##############################################################################

## 2a. For building genome indexes

# # this just moves any updated gDNA scripts from $WD to sample dirs
# while IFS= read -r i; do
#   cp $ATACscript2 $i
# done < $prefixgDNA

# A. download required genomes (already stored in workarea but having permission problems!)

software
WD=(/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/data/ATACseq/3.run2) # if you change path here then ensure to change in all other scripts
genomesdir=($WD/genomes)
cd $WD

# Metriaclima zebra
mkdir -p $genomesdir/Mzebra
wget --recursive --no-host-directories --cut-dirs=3 ftp://ftp.ensembl.org/pub/current_fasta/maylandia_zebra -P $genomesdir/Mzebra
wget --recursive --no-host-directories --cut-dirs=1 ftp://ftp.ensembl.org/pub/{current_gff3,current_gtf}/maylandia_zebra -P $genomesdir/Mzebra

# Pundamilia nyererei
mkdir -p $genomesdir/Pnyererei
wget --recursive --no-host-directories --cut-dirs=3 ftp://ftp.ensembl.org/pub/current_fasta/pundamilia_nyererei -P $genomesdir/Pnyererei
wget --recursive --no-host-directories --cut-dirs=1 ftp://ftp.ensembl.org/pub/{current_gff3,current_gtf}/pundamilia_nyererei -P $genomesdir/Pnyererei

# Astatotilapia burtoni
mkdir -p $genomesdir/Aburtoni
wget --recursive --no-host-directories --cut-dirs=3 ftp://ftp.ensembl.org/pub/current_fasta/haplochromis_burtoni -P $genomesdir/Aburtoni
wget --recursive --no-host-directories --cut-dirs=1 ftp://ftp.ensembl.org/pub/{current_gff3,current_gtf}/haplochromis_burtoni -P $genomesdir/Aburtoni

# Neolamprologus brichardi
mkdir -p $genomesdir/Nbrichardi
wget --recursive --no-host-directories --cut-dirs=3 ftp://ftp.ensembl.org/pub/current_fasta/neolamprologus_brichardi -P $genomesdir/Nbrichardi
wget --recursive --no-host-directories --cut-dirs=1 ftp://ftp.ensembl.org/pub/{current_gff3,current_gtf}/neolamprologus_brichardi -P $genomesdir/Nbrichardi

# Astatotilapia calliptera
mkdir -p $genomesdir/Acalliptera
wget --recursive --no-host-directories --cut-dirs=3 ftp://ftp.ensembl.org/pub/current_fasta/astatotilapia_calliptera -P $genomesdir/Acalliptera
wget --recursive --no-host-directories --cut-dirs=1 ftp://ftp.ensembl.org/pub/{current_gff3,current_gtf}/astatotilapia_calliptera -P $genomesdir/Acalliptera

# Oreochromis niloticus
mkdir -p $genomesdir/Oniloticus
wget --recursive --no-host-directories --cut-dirs=3 ftp://ftp.ensembl.org/pub/current_fasta/oreochromis_niloticus -P $genomesdir/Oniloticus
wget --recursive --no-host-directories --cut-dirs=1 ftp://ftp.ensembl.org/pub/{current_gff3,current_gtf}/oreochromis_niloticus -P $genomesdir/Oniloticus

exit

cat $mzeb/dna/Maylandia_zebra.M_zebra_UMD2a.dna.primary_assembly.LG*.fa.gz > $mzeb/dna/Maylandia_zebra.M_zebra_UMD2a.dna.primary_assembly.allLG.fa.tmp1.gz
zcat $mzeb/dna/Maylandia_zebra.M_zebra_UMD2a.dna.primary_assembly.allLG.fa.tmp1.gz | awk '{print $1}' | gzip > $mzeb/dna/Maylandia_zebra.M_zebra_UMD2a.dna.primary_assembly.allLG.fa.gz
rm $mzeb/dna/Maylandia_zebra.M_zebra_UMD2a.dna.primary_assembly.allLG.fa.tmp1.gz
gunzip -c $mzeb/dna/Maylandia_zebra.M_zebra_UMD2a.dna.primary_assembly.allLG.fa.gz > $mzeb/dna/Maylandia_zebra.M_zebra_UMD2a.dna.primary_assembly.allLG.fa

zcat $pnyegenome.gz | awk '{print $1}' > $pnyegenome
zcat $aburgenome.gz | awk '{print $1}' > $aburgenome
zcat $nbrigenome.gz | awk '{print $1}' > $nbrigenome
zcat $acalgenome.gz | awk '{print $1}' > $acalgenome

cat $acal/dna/Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.*.fa.gz > $acal/dna/Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.fa.tmp1.gz
zcat $acal/dna/Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.fa.tmp1.gz | awk '{print $1}' | gzip > $acal/dna/Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.fa.gz
rm $acal/dna/Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.fa.tmp1.gz
gunzip -c $acal/dna/Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.fa.gz > $acal/dna/Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.fa

cat $onil/dna/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.primary_assembly.LG*.fa.gz > $onil/dna/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.primary_assembly.allLG.fa.tmp1.gz
zcat $onil/dna/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.primary_assembly.allLG.fa.tmp1.gz | awk '{print $1}' | gzip > $onil/dna/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.primary_assembly.allLG.fa.gz
rm $onil/dna/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.primary_assembly.allLG.fa.tmp1.gz
gunzip -c $onil/dna/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.primary_assembly.allLG.fa.gz > $onil/dna/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.primary_assembly.allLG.fa

##############################################################################

## 2b. Trimming and alignment of gDNA reads

# # Example
# cd $WD/Ab5_L_gDNA # DONE
# sbatch ATAC_Bioinf_pipeline_v2b_gDNA.sh -s Ab5_L_gDNA -g AstBur1.0 -f /tgac/workarea/group-vh/Tarang/Reference_Genomes/ensembl/cichlids/Aburtoni/dna/Haplochromis_burtoni.AstBur1.0.dna.nonchromosomal.fa

# A. Create a prefix file to iterate for running scripts
awk -F'_' '{print $1}' $prefixgDNA | sed 's/1a//g' | sed 's/1b//g' | sed 's/2a//g' | sed 's/2b//g' | sed 's/3a//g' | sed 's/3b//g' | sed 's/[0-9]//g' | sed 's/Pnm/Pn/g' > $prefixgDNA.tmp1 # strip all other characters to expose only species ID
awk '{print $1}' $prefixgDNA > $prefixgDNA.tmp2
paste $prefixgDNA.tmp1 $prefixgDNA.tmp2 > $prefixgDNA2 # create a new file that has species ID alongside the sample ID
rm $prefixgDNA.tmp1 $prefixgDNA.tmp2

# B. While loop with if statements to appropriately spawn off sbatch scripts
while read -r i1 i2; do
  if [[ "$i1" == "$sp1" ]]; then
    # echo "New if statement: "$i1 $i2 $sp1
    cd $WD/$i2
    # echo "sbatch $ATACscript2 -s $i2 -g $sp1ID -f $mzebgenome"
    sbatch $ATACscript2 -s $i2 -g $sp1ID -f $mzebgenome
  fi
  if [[ "$i1" == "$sp2" ]]; then
    # echo "New if statement: "$i1 $i2 $sp2
    cd $WD/$i2
    # echo "sbatch $ATACscript2 -s $i2 -g $sp2ID -f $pnyegenome"
    sbatch $ATACscript2 -s $i2 -g $sp2ID -f $pnyegenome
  fi
  if [[ "$i1" == "$sp3" ]]; then
    # echo "New if statement: "$i1 $i2 $sp3
    cd $WD/$i2
    # echo "sbatch $ATACscript2 -s $i2 -g $sp3ID -f $aburgenome"
    sbatch $ATACscript2 -s $i2 -g $sp3ID -f $aburgenome
  fi
  if [[ "$i1" == "$sp4" ]]; then
    # echo "New if statement: "$i1 $i2 $sp4
    cd $WD/$i2
    # echo "sbatch $ATACscript2 -s $i2 -g $sp4ID -f $nbrigenome"
    sbatch $ATACscript2 -s $i2 -g $sp4ID -f $nbrigenome
  fi
  if [[ "$i1" == "$sp5" ]]; then
    # echo "New if statement: "$i1 $i2 $sp5
    cd $WD/$i2
    # echo "sbatch $ATACscript2 -s $i2 -g $sp5ID -f $onilgenome"
    sbatch $ATACscript2 -s $i2 -g $sp5ID -f $onilgenome
  fi
  if [[ "$i1" == "$sp6" ]]; then
    # echo "New if statement: "$i1 $i2 $sp6
    cd $WD/$i2
    # echo "sbatch $ATACscript2 -s $i2 -g $sp6ID -f $acalgenome"
    sbatch $ATACscript2 -s $i2 -g $sp6ID -f $acalgenome
  fi
done < $prefixgDNA2 # this while loop will run gDNA script in each gDNA folder

# C. Check whether the output files have written their completion status ok:
cd $WD
slurmoutstotal=$(ls -1 *_*_gDNA/slurm*.out | wc -l)
completeslurm=$(ls -1 *_*_gDNA/2.read_alignment/gDNAcompleted.txt | wc -l)
if [[ $slurmoutstotal -eq $completeslurm ]]
then
  echo -e "All gDNA library alignments have completed successfully....\nTotal runs = $slurmoutstotal\nTotal completed runs = $completeslurm"
else
  echo -e "Not all gDNA library alignments have completed successfuly....\nTotal runs = $slurmoutstotal\nTotal completed runs = $completeslurm"
fi

# D. Check and print all gDNA mapping rates

for f in *_*_gDNA/2.read_alignment/*.align.log; do
  sample=$(echo $f | sed 's/_gDNA\/2.read_alignment.*/_gDNA/g')
  ar=$(grep 'overall alignment rate' $f | sed 's/ overall alignment rate//g')
  echo -e $sample'\t'$ar >> gDNA_alignment_stats_mainsummary.txt
done

for f in *_*_gDNA/2.read_alignment/*.align.log; do
  echo '==================================================' >> gDNA_alignment_stats_summary.txt
  echo $f >> gDNA_alignment_stats_summary.txt
  grep 'overall alignment rate' $f >> gDNA_alignment_stats_summary.txt
  echo '==================================================' >> gDNA_alignment_stats_summary.txt
done
for i in *_*_gDNA/2.read_alignment/*_flagstat_qc1.txt; do
  echo '==================================================' >> gDNA_alignment_stats_summary.txt
  echo $i >> gDNA_alignment_stats_summary.txt
  grep 'QC-passed reads' $i >> gDNA_alignment_stats_summary.txt
  grep 'mapped (' $i >> gDNA_alignment_stats_summary.txt
  echo '==================================================' >> gDNA_alignment_stats_summary.txt
done

##############################################################################

## 3. Trimming and alignment of ATAC reads, and then filtering, calling peaks, and annotation

# # Example
# cd $WD/Ab5_L_ATAC
# sbatch ATAC_Bioinf_pipeline_v2b.sh -s Ab5_L_ATAC -g AstBur1.0 -f /tgac/workarea/group-vh/Tarang/Reference_Genomes/ensembl/cichlids/Aburtoni/dna/Haplochromis_burtoni.AstBur1.0.dna.nonchromosomal.fa -m NC_027289.1 -u mehtat -a /tgac/workarea/group-vh/Tarang/Reference_Genomes/ensembl/cichlids/Aburtoni/current_gtf/haplochromis_burtoni/Haplochromis_burtoni.AstBur1.0.100.gtf.gz

# A. Create a prefix file to iterate for running scripts
awk -F'_' '{print $1}' $prefixATAC | sed 's/1a//g' | sed 's/1b//g' | sed 's/2a//g' | sed 's/2b//g' | sed 's/3a//g' | sed 's/3b//g' | sed 's/[0-9]//g' | sed 's/Pnm/Pn/g' > $prefixATAC.tmp1 # strip all other characters to expose only species ID
awk '{print $1}' $prefixATAC > $prefixATAC.tmp2
paste $prefixATAC.tmp1 $prefixATAC.tmp2 > $prefixATAC2 # create a new file that has species ID alongside the sample ID
rm $prefixATAC.tmp1 $prefixATAC.tmp2

# B. While loop with if statements to appropriately spawn off sbatch scripts
while read -r i1 i2; do
  if [[ "$i1" == "$sp1" ]]; then
    # echo "New if statement: "$i1 $i2 $sp1
    cd $WD/$i2
    # echo "sbatch $ATACscript3 -s $i2 -g $sp1ID -f $mzebgenome -m $sp1mitID -u $user -a $mzebgtf"
    sbatch $ATACscript3 -s $i2 -g $sp1ID -f $mzebgenome -m $sp1mitID -u $user -a $mzebgtf
  fi
  if [[ "$i1" == "$sp2" ]]; then
    # echo "New if statement: "$i1 $i2 $sp2
    cd $WD/$i2
    # echo "sbatch $ATACscript3 -s $i2 -g $sp2ID -f $pnyegenome -m $sp2mitID -u $user -a $pnyegtf"
    sbatch $ATACscript3 -s $i2 -g $sp2ID -f $pnyegenome -m $sp2mitID -u $user -a $pnyegtf
  fi
  if [[ "$i1" == "$sp3" ]]; then
    # echo "New if statement: "$i1 $i2 $sp3
    cd $WD/$i2
    # echo "sbatch $ATACscript3 -s $i2 -g $sp3ID -f $aburgenome -m $sp3mitID -u $user -a $aburgtf"
    sbatch $ATACscript3 -s $i2 -g $sp3ID -f $aburgenome -m $sp3mitID -u $user -a $aburgtf
  fi
  if [[ "$i1" == "$sp4" ]]; then
    # echo "New if statement: "$i1 $i2 $sp4
    cd $WD/$i2
    # echo "sbatch $ATACscript3 -s $i2 -g $sp4ID -f $nbrigenome -m $sp4mitID -u $user -a $nbrigtf"
    sbatch $ATACscript3 -s $i2 -g $sp4ID -f $nbrigenome -m $sp4mitID -u $user -a $nbrigtf
  fi
  if [[ "$i1" == "$sp5" ]]; then
    # echo "New if statement: "$i1 $i2 $sp5
    cd $WD/$i2
    # echo "sbatch $ATACscript3 -s $i2 -g $sp5ID -f $onilgenome -m $sp5mitID -u $user -a $onilgtf"
    sbatch $ATACscript3 -s $i2 -g $sp5ID -f $onilgenome -m $sp5mitID -u $user -a $onilgtf
  fi
  if [[ "$i1" == "$sp6" ]]; then
    # echo "New if statement: "$i1 $i2 $sp6
    cd $WD/$i2
    # echo "sbatch $ATACscript3 -s $i2 -g $sp6ID -f $acalgenome -m $sp6mitID -u $user -a $acalgtf"
    sbatch $ATACscript3 -s $i2 -g $sp6ID -f $acalgenome -m $sp6mitID -u $user -a $acalgtf
  fi
done < $prefixATAC2 # this while loop will run ATAC script in each ATAC folder

## checking jobs and the samples associated (change the grep accordingly)
cd $WD
jobids=$(squeue -u mehtat | grep '5.peakca' | awk -F' ' '{print $1}' | awk 'BEGIN { ORS = " " } { print }') # this will store the jobids into a single row variable to use in for loop

for i in $jobids; do
  com=$(scontrol show jobid -dd $i | grep 'Command=' )
  err=$(scontrol show jobid -dd $i | grep 'StdErr=')
  echo -e $com'\t'$err
done # this will output two columns - first col shows the command where error occured and second shows the stout error file associated. The folder ran before the command will invariably have the error with it.

############################################ ~~~~~~~ ############################################
## several jobs failed due to insufficent space > move completed jobs to scratch area (except corresponding gDNA required for ATAC runs)
cd $WD
mkdir -p $WD2
failedids=($WD/failed_IDs.txt)
failedidsATAC=($WD/failed_IDs_ATAC.txt)
prefixcomb=($WD/prefix_all.txt)
successids=($WD/successful_IDs.txt)
prefixfailed=($WD/prefix_failed.txt)
successids2=($WD/successful_IDsv2.txt)

for i in $jobids; do
  workdirATAC=$(scontrol show jobid -dd $i | grep 'WorkDir=' | sed 's/5.peak_calling//g' | sed 's/WorkDir=//g' | sed "s/\/ei\/projects\/9\/9e238063-c905-4076-a975-f7c7f85dbd56\/data\/ATACseq\/3.run2\///g" | sed 's/\///g')
  workdirgDNA=$(echo $workdirATAC | sed 's/_ATAC/_gDNA/g' )
  echo $workdirATAC >> $failedids
  echo $workdirgDNA >> $failedids
done # this will provide the sample IDs for those that failed plus corresponding gDNA folders

cat $prefixgDNA $prefixATAC > $prefixcomb

grep -vf $failedids $prefixcomb > $successids # create a file to iterate sample IDs that finished successfuly
# grep -v '_gDNA' $successids | grep -v '1aAc_3dpf_ATAC' | grep -v '2aAc_7dpf_ATAC' > $successids2 # this removes any already done

echo '#!/bin/bash -e' > mv_files.sh
echo '#SBATCH -p tgac-short # partition (queue)' >> mv_files.sh
echo '#SBATCH -N 1 # number of nodes' >> mv_files.sh
echo '#SBATCH -n 1 # number of tasks' >> mv_files.sh
echo '#SBATCH --array=0-101' >> mv_files.sh
echo '#SBATCH --mem-per-cpu 10000' >> mv_files.sh
echo '#SBATCH -t 0-00:45' >> mv_files.sh
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> mv_files.sh
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> mv_files.sh
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> mv_files.sh
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> mv_files.sh
printf '\n' >> mv_files.sh
echo 'source rsync-3.1.1' >> mv_files.sh
echo "mapfile -t successids < $successids" >> mv_files.sh
echo "#mapfile -t successids < $successids2" >> mv_files.sh
echo '#mv ${successids[${SLURM_ARRAY_TASK_ID}]}'" $WD2" >> mv_files.sh
echo '#mv ${successids[${SLURM_ARRAY_TASK_ID}]}'" $WD2" >> mv_files.sh
echo 'rsync -avhP --checksum ${successids[${SLURM_ARRAY_TASK_ID}]} '" $WD2" >> mv_files.sh

sbatch mv_files.sh # run the above

rm -r $WD2/On3_B_ATAC/{1.adaptor_trimming,2.read_alignment,3.Mtfilt_fragcnt,4.postalign_filt,5.peak_calling} # needs a full re-run
rm -r $WD2/On3_T_ATAC/{1.adaptor_trimming,2.read_alignment,3.Mtfilt_fragcnt,4.postalign_filt,5.peak_calling} # needs a re-run (peak calling only but re-run all)
mv $WD2/On3_B_ATAC $WD
mv $WD2/On3_T_ATAC $WD
sbatch $ATACscript3 -s On3_B_ATAC -g $sp5ID -f $onilgenome -m $sp5mitID -u $user -a $onilgtf
sbatch $ATACscript3 -s On3_T_ATAC -g $sp5ID -f $onilgenome -m $sp5mitID -u $user -a $onilgtf
rm -r Pnm4_E_ATAC/{1.adaptor_trimming,2.read_alignment,3.Mtfilt_fragcnt,4.postalign_filt,5.peak_calling}
sbatch $ATACscript3 -s Pnm4_E_ATAC -g $sp2ID -f $pnyegenome -m $sp2mitID -u $user -a $pnyegtf

while read -r h; do
  rm -r $h
done < $successids # once the above is done, remove the folders from $WD

# remove failed files folders - best to start again!
grep ATAC $failedids > $failedidsATAC

while read -r g; do
  rm -r $WD/${g}/{1.adaptor_trimming,2.read_alignment,3.Mtfilt_fragcnt,4.postalign_filt,5.peak_calling}
  rm $WD/${g}/slurm*
done < $failedidsATAC

# re-run
grep -f $failedids $prefixATAC2 > $prefixfailed # create a new file to iterate for sample IDs

# While loop with if statements to appropriately spawn off sbatch scripts
while read -r i1 i2; do
  if [[ "$i1" == "$sp1" ]]; then
    # echo "New if statement: "$i1 $i2 $sp1
    cd $WD/$i2
    # echo "sbatch $ATACscript3 -s $i2 -g $sp1ID -f $mzebgenome -m $sp1mitID -u $user -a $mzebgtf"
    sbatch $ATACscript3 -s $i2 -g $sp1ID -f $mzebgenome -m $sp1mitID -u $user -a $mzebgtf
  fi
  if [[ "$i1" == "$sp2" ]]; then
    # echo "New if statement: "$i1 $i2 $sp2
    cd $WD/$i2
    # echo "sbatch $ATACscript3 -s $i2 -g $sp2ID -f $pnyegenome -m $sp2mitID -u $user -a $pnyegtf"
    sbatch $ATACscript3 -s $i2 -g $sp2ID -f $pnyegenome -m $sp2mitID -u $user -a $pnyegtf
  fi
  if [[ "$i1" == "$sp3" ]]; then
    # echo "New if statement: "$i1 $i2 $sp3
    cd $WD/$i2
    # echo "sbatch $ATACscript3 -s $i2 -g $sp3ID -f $aburgenome -m $sp3mitID -u $user -a $aburgtf"
    sbatch $ATACscript3 -s $i2 -g $sp3ID -f $aburgenome -m $sp3mitID -u $user -a $aburgtf
  fi
  if [[ "$i1" == "$sp4" ]]; then
    # echo "New if statement: "$i1 $i2 $sp4
    cd $WD/$i2
    # echo "sbatch $ATACscript3 -s $i2 -g $sp4ID -f $nbrigenome -m $sp4mitID -u $user -a $nbrigtf"
    sbatch $ATACscript3 -s $i2 -g $sp4ID -f $nbrigenome -m $sp4mitID -u $user -a $nbrigtf
  fi
  if [[ "$i1" == "$sp5" ]]; then
    # echo "New if statement: "$i1 $i2 $sp5
    cd $WD/$i2
    # echo "sbatch $ATACscript3 -s $i2 -g $sp5ID -f $onilgenome -m $sp5mitID -u $user -a $onilgtf"
    sbatch $ATACscript3 -s $i2 -g $sp5ID -f $onilgenome -m $sp5mitID -u $user -a $onilgtf
  fi
  if [[ "$i1" == "$sp6" ]]; then
    # echo "New if statement: "$i1 $i2 $sp6
    cd $WD/$i2
    # echo "sbatch $ATACscript3 -s $i2 -g $sp6ID -f $acalgenome -m $sp6mitID -u $user -a $acalgtf"
    sbatch $ATACscript3 -s $i2 -g $sp6ID -f $acalgenome -m $sp6mitID -u $user -a $acalgtf
  fi
done < $prefixfailed # this while loop will run ATAC script in each ATAC folder

## the following are struggling to complete:
# Pnm4_T_ATAC
# # /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/data/ATACseq/3.run2/Pnm4_T_ATAC/4.postalign_filt/ #sambamba-sort: Unable to write to stream - GIVE IT MORE MEMORY
# nano /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/data/ATACseq/3.run2/Pnm4_T_ATAC/4.postalign_filt/4.postalign_filt.sh # change memory to 120G
# scancel 28283356
# scancel 28296911
# cd /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/data/ATACseq/3.run2/Pnm4_T_ATAC/4.postalign_filt/
# sbatch 4.postalign_filt.sh # if that runs successfuly then run peak calling and bed to bigbed conversion for narrowpeaks separately too

# these are also struggling at the same stage
scancel -n ATAC_Bioinf_pipeline_v2b.sh
scancel -n 5.peakcall.sh

# these all failed for 'sambamba-sort: /tmp/sambamba-pid19093-imkg: No space left on device'
# create tmp dir for each sample in a place with lots of space - /tgac/workarea/group-vh/Tarang/ATACseq/3.run2/${i}_tmp
# up the memory to 120000 (in sbatch header and sambamba sort)
  # sambamba sort -m 120G -t 1 --tmpdir /tgac/workarea/group-vh/Tarang/ATACseq/3.run2/${i}_tmp -o $bam_file_sorted -u $bam_file
  # sambamba markdup -l 0 -t 1 --tmpdir /tgac/workarea/group-vh/Tarang/ATACseq/3.run2/${i}_tmp $bam_file_sorted $bam_file_dup

# edited the memory and tmp dir and ran using sbatch > that ran successfuly > then run peak calling and bed to bigbed conversion for narrowpeaks separately too
for i in On3_T_ATAC On3_B_ATAC On2_G_ATAC Pnm3_B_ATAC Pnm4_T_ATAC Mz2_T_ATAC 1bAc_3dpf_ATAC On1_G_ATAC; do
  # ls -tlrh $i/4.postalign_filt/
  # rm $i/4.postalign_filt/slurm*
  # echo $i
  # head $i/4.postalign_filt/slurm*.err
  # mkdir /tgac/workarea/group-vh/Tarang/ATACseq/3.run2/${i}_tmp
  # sed -i 's/88000/120000/g' $i/4.postalign_filt/4.postalign_filt.sh
  # sed -i 's/88G/120G/g' $i/4.postalign_filt/4.postalign_filt.sh
  # nano $i/4.postalign_filt/4.postalign_filt.sh
  # less $i/4.postalign_filt/4.postalign_filt.sh
  # sbatch $i/4.postalign_filt/4.postalign_filt.sh
  ls -tlrh $i/5.peak_calling/
  # cd $WD/${i}/5.peak_calling/
  # sbatch 5.peakcall.sh
done

for i in Ab5_L_ATAC Ab6_L_ATAC; do
  ls -tlrh $i/4.postalign_filt/
  # rm $i/4.postalign_filt/slurm*
  # echo $i
  # head $i/4.postalign_filt/slurm*.err
  # mkdir /tgac/workarea/group-vh/Tarang/ATACseq/3.run2/${i}_tmp
  # sed -i 's/88000/120000/g' $i/4.postalign_filt/4.postalign_filt.sh
  # sed -i 's/88G/120G/g' $i/4.postalign_filt/4.postalign_filt.sh
  # nano $i/4.postalign_filt/4.postalign_filt.sh
  # less $i/4.postalign_filt/4.postalign_filt.sh
  # cd $WD/${i}/4.postalign_filt
  # sbatch 4.postalign_filt.sh
  # ls -tlrh $i/5.peak_calling/
  # cd $WD/${i}/5.peak_calling/
  # sbatch 5.peakcall.sh
done
# the above failed to tmp isssues on mark-dup - re-ran supressing the sort as that is completed

# > On1_G_ATAC # restarted peak calling due to time limit {{{{{COMPLETED}}}}}

# > Ab5_L_ATAC and Ab6_L_ATAC failed again at stage 4 for the same reason:
# sambamba-markdup: Cannot open file `/tgac/workarea/group-vh/Tarang/ATACseq/3.run2/Ab6_L_ATAC_tmp/sambamba-pid48145-markdup-jfzq/PairedEndsInfoxzwx16' in mode `w+' (Too many open files)
# added '--overflow-list-size 600000' to the sambamba markdup command
# re-ran step 4 - Ab5_L_ATAC completed {{{{{CHECK Ab6_L_ATAC step4}}}}}
# amended peak call time limit and re-ran peak calling for: {{{{{{{RUNNING - Ab5_L_ATAC and Ab6_L_ATAC}}}}}}}

# some of the other jobs (moved to scratch) didn't finish (genrich) due to time limit:
for i in *_*_ATAC/5.peak_calling; do
  ls -tlrh $i
  # sample=$(echo $i | sed 's/_ATAC\/5.peak_calling.*/_ATAC/g')
  # time=$(grep 'DUE TO TIME LIMIT' $i)
  # echo -e $sample
done

# these need re-running peakcalling with longer run time: Pnm1_L_ATAC Nb5_T_ATAC Nb5_L_ATAC Nb4_T_ATAC Nb4_L_ATAC Ab5_T_ATAC
for i in Pnm1_L_ATAC Nb5_T_ATAC Nb5_L_ATAC Nb4_T_ATAC Nb4_L_ATAC Ab5_T_ATAC; do
  # rm ${i}/5.peak_calling/slurm*
  # rm ${i}/5.peak_calling/Pnm*
  # rm ${i}/5.peak_calling/Nb*
  # rm ${i}/5.peak_calling/Ab*
  # rm ${i}/5.peak_calling/*.bed
  # rm ${i}/5.peak_calling/*.as
  # sed -i 's/0-04:59/0-07:59/g' ${i}/5.peak_calling/5.peakcall.sh # change time to 0-07:59
  # sed -i 's/data/scratch/g' ${i}/5.peak_calling/5.peakcall.sh
  # sed -i 's/scratch\/ATACseq\/3.run2\/genomes/data\/ATACseq\/3.run2\/genomes/g' ${i}/5.peak_calling/5.peakcall.sh
  # cd $WD2/${i}/5.peak_calling
  ls -tlrh $WD2/${i}/5.peak_calling
  # sbatch 5.peakcall.sh ###### {{{{{CHECK}}}}}
done

# the above ran ok except for Nb4_L_ATAC (due to gDNA folder incorrect location) and Nb4_T_ATAC (time limit)
cd $WD2/Nb4_L_ATAC/5.peak_calling; sbatch 5.peakcall.sh # re-ran
cd $WD2/Nb4_T_ATAC/5.peak_calling; sbatch 5.peakcall.sh # re-ran

## last few to check once complete
ls -tlrh $WD/Ab5_L_ATAC/5.peak_calling
ls -tlrh $WD/Ab6_L_ATAC/5.peak_calling
ls -tlrh $WD2/Nb4_L_ATAC/5.peak_calling
ls -tlrh $WD2/Nb4_T_ATAC/5.peak_calling


# bed to bigbed conversion
source ucsc_utils-v333

for i in On3_T_ATAC On3_B_ATAC On2_G_ATAC Pnm3_B_ATAC Pnm4_T_ATAC Mz2_T_ATAC 1bAc_3dpf_ATAC On1_G_ATAC Ab5_L_ATAC Ab6_L_ATAC; do
  cd $WD/${i}/5.peak_calling
  bigbed=(${i}/5.peak_calling/${i}_peaks.narrowPeak.bb)
  peak=(${i}/5.peak_calling/${i}_peaks.narrowPeak)
  peakgz=(${i}/5.peak_calling/${i}_peaks.narrowPeak.gz)
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
  zcat ${peakgz} | sort -k1,1 -k2,2n > ${bigbed}.tmp
  bedClip ${bigbed}.tmp ${scafflen2} ${bigbed}.tmp2
  bedToBigBed -type=bed6+4 -as=narrowPeak.as ${bigbed}.tmp2 ${scafflen2} ${bigbed}
  rm -f ${bigbed}.tmp ${bigbed}.tmp2
done

############################################ ~~~~~~~ ############################################



# C. Check whether the output files have written their completion status ok:
cd $WD
ATACslurmoutstotal=$(ls -1 *_*_ATAC/slurm*.out | wc -l)
ATACcompleteslurm=$(grep 'EXITING SCRIPT' *_*_ATAC/slurm*.out | wc -l)
if [[ $ATACslurmoutstotal -eq $ATACcompleteslurm ]]
then
  echo -e "All ATAC library alignments to peak calling and bigbed conversion have completed successfully....\nTotal runs = $ATACslurmoutstotal\nTotal completed runs = $ATACcompleteslurm"
else
  echo -e "Not all ATAC library alignments to peak calling and bigbed conversion have completed successfuly....\nTotal runs = $ATACslurmoutstotal\nTotal completed runs = $ATACcompleteslurm"
fi

# D. Check and print all ATAC mapping rates

for f in *_*_ATAC/2.read_alignment/*.align.log; do
  sample=$(echo $f | sed 's/_ATAC\/2.read_alignment.*/_ATAC/g')
  ar=$(grep 'overall alignment rate' $f | sed 's/ overall alignment rate//g')
  echo -e $sample'\t'$ar >> ATAC_alignment_stats_mainsummary.txt
done

for f in *_*_ATAC/2.read_alignment/*.align.log; do
  echo '==================================================' >> ATAC_alignment_stats_summary.txt
  echo $f >> ATAC_alignment_stats_summary.txt
  grep 'overall alignment rate' $f >> ATAC_alignment_stats_summary.txt
  echo '==================================================' >> ATAC_alignment_stats_summary.txt
done
for i in *_*_ATAC/2.read_alignment/*_flagstat_qc1.txt; do
  echo '==================================================' >> ATAC_alignment_stats_summary.txt
  echo $i >> ATAC_alignment_stats_summary.txt
  grep 'QC-passed reads' $i >> ATAC_alignment_stats_summary.txt
  grep 'mapped (' $i >> ATAC_alignment_stats_summary.txt
  echo '==================================================' >> ATAC_alignment_stats_summary.txt
done

# E. report number of narrow peaks - MACS2
for file in *_*_ATAC/5.peak_calling/*_*_ATAC_peaks.narrowPeak.gz; do
  sample=$(echo $file | sed 's/_ATAC\/5.peak_calling.*/_ATAC/g')
  peaks=$(zcat $file | wc -l)
  echo -e $sample'\t'$peaks >> ATAC_narrowpeaks_mainsummary.txt
done

for file in *_*_ATAC/5.peak_calling/*_*_ATAC_peaks.narrowPeak; do
  sample=$(echo $file | sed 's/_ATAC\/5.peak_calling.*/_ATAC/g')
  peaks=$(wc -l $file | awk -F' ' '{print $1}')
  echo -e $sample'\t'$peaks >> ATAC_narrowpeaks_mainsummary.txt
done

# E. report number of Genrich peaks
for file in *_*_ATAC/5.peak_calling/*_*_Genrich.peaks; do
  sample=$(echo $file | sed 's/_ATAC\/5.peak_calling.*/_ATAC/g')
  peaks=$(wc -l $file | awk -F' ' '{print $1}')
  echo -e $sample'\t'$peaks >> ATAC_Genrichpeaks_mainsummary.txt
done


##############################################################################

## 4. ATAC downstream analysis - IDR, footprinting, peak annotation and differential analysis

##############################################################################

# Before you start running this section, consider creating a tar ball of the files that you will not need. This includes:
# *_*_gDNA
# *_*_ATAC/{1.adaptor_trimming,2.read_alignment,3.Mtfilt_fragcnt,4.postalign_filt,5.peak_calling}
