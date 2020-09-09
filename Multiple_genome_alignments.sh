#!/bin/sh

##############################################################################
### Sept 2020 - Multiple genome alignments
# USAGE: ./alignment.local.sh [reference].fasta
# REQUIREMENTS: Run in its own directory with nothing but *.fasta files
# SOFTWARE IT USES:
#       Sources mugsy
#       Uses nucmer (better to use 64-bit compiled version for larger genomes)
#       delta2maf (from mugsy MUMmer-3.20 package)
#       maf_sort from tba/multiz
#       single_cov2 from tba/multiz
#       multiz from tba/multiz
##############################################################################

### Instructions:
# 1 Copy this directory (/ei/scratch/craine/Tarang_align/) into your own space.
# 2. Copy or soft link the genomes you want to align into that directory.  I prefer to soft link in case my script destroys your files.
    # Your genomes must end in .fasta
    # I recommend keeping your genome names simple, i.e. the species name.  No underscores or dots or anything.  Such as Oleucostictus.fasta
    # Decide which of your species you want to be the overall reference that the others are mapped to.
# 4. Simply run ./alignment.local.sh [reference].fasta

##############################################################################
WD=(/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/1.MultiGenAlign_MzPnAbNbAc_OnREF) # if you change path here then ensure to change in all other scripts
copyscripts=(/ei/scratch/craine/Tarang_align/)

# Genomes
genomesdir=(/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/data/ATACseq/3.run2/genomes)

mzeb=($genomesdir/Mzebra)
mzebgenome=$mzeb/dna/Maylandia_zebra.M_zebra_UMD2a.dna.primary_assembly.allLG.fa
mzebgenome2=$WD/MetzebUMD2a.fasta

pnye=($genomesdir/Pnyererei)
pnyegenome=$pnye/dna/Pundamilia_nyererei.PunNye1.0.dna.nonchromosomal.fa
pnyegenome2=$WD/PunNye1.fasta

abur=($genomesdir/Aburtoni)
aburgenome=$abur/dna/Haplochromis_burtoni.AstBur1.0.dna.nonchromosomal.fa
aburgenome2=$WD/AstBur1.fasta

nbri=($genomesdir/Nbrichardi)
nbrigenome=$nbri/dna/Neolamprologus_brichardi.NeoBri1.0.dna.nonchromosomal.fa
nbrigenome2=$WD/NeoBri1.fasta

acal=($genomesdir/Acalliptera)
acalgenome=$acal/dna/Astatotilapia_calliptera.fAstCal1.2.dna.primary_assembly.allLG.fa
acalgenome2=$WD/AstCal12.fasta

onil=($genomesdir/Oniloticus)
onilgenome=$onil/dna/Oreochromis_niloticus.O_niloticus_UMD_NMBU.dna.primary_assembly.allLG.fa
onilgenome2=$WD/OrenilUMDNMBU.fasta

runscript=(alignment.local.sh)
cores=16
mem=512GB
sbatchrunscript=(sbatch_alignment.local.sh)
ref=OrenilUMDNMBU.fasta # insert the reference fasta you want (not absolute path)

##############################################################################
# 1 Copy this directory (/ei/scratch/craine/Tarang_align/) into your own space.

mkdir -p $WD
cd $WD

cp -r $copyscripts/* .

##############################################################################
# 2. Copy or soft link the genomes you want to align into that directory.  I prefer to soft link in case my script destroys your files.
    # Your genomes must end in .fasta
    # I recommend keeping your genome names simple, i.e. the species name.  No underscores or dots or anything.  Such as Oleucostictus.fasta

cp $mzebgenome $mzebgenome2
cp $pnyegenome $pnyegenome2
cp $aburgenome $aburgenome2
cp $nbrigenome $nbrigenome2
cp $acalgenome $acalgenome2
cp $onilgenome $onilgenome2

##############################################################################
# 3. Simply run ./alignment.local.sh [reference].fasta
# For the ATAC-seq work we want O. niloticus as the reference - basal and best genome

echo '#!/bin/bash -e' > $sbatchrunscript
echo '#SBATCH -p ei-largemem # partition (queue)' >> $sbatchrunscript
echo '#SBATCH -N 1 # number of nodes' >> $sbatchrunscript
echo "#SBATCH -c $cores # number of cores" >> $sbatchrunscript
echo "#SBATCH --mem $mem # memory pool for all cores" >> $sbatchrunscript
echo '#SBATCH -o slurm.%N.%j.out # STDOUT' >> $sbatchrunscript
echo '#SBATCH -e slurm.%N.%j.err # STDERR' >> $sbatchrunscript
echo '#SBATCH --mail-type=ALL # notifications for job done & fail' >> $sbatchrunscript
echo '#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address' >> $sbatchrunscript
printf '\n' >> $sbatchrunscript
echo "./$runscript $ref" >> $sbatchrunscript

sbatch $sbatchrunscript # run the above

# If you need to re-run then rename the original files of the fasta (that got changed) before re-running
# for file in *original; do mv $file ${file%.original}; done

##############################################################################
