########################## June 2022 ####################################
#
# Enrichment of motifs in gene promoter peaks across tissues and species
#
#########################################################################

###### Summary of what you need and need to do:
## 1. Prepare a motif ID and TF name list with col1-motifID and col2-motif_name for ALL input motifs used for calling TF footprints (all found here: /hpc-home/mehtat/rgtdata/motifs/*.mtf) > equivalent to ~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/data/JASPAR_Vertebrates_2018_motifname.txt > $MOTIFNAME
## 2. Using the *_ATAC_mpbs.bed output from TF footprinting, filter all TF footprints for those falling in gene promoter regions and associate with gene (based on gene promoter)
## 3. Create a two column file with col1-ensemblID and col2-tissue for the above but for each species > equivalent to ~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/results/k10_arb_output/prediction/Puny_speciesspecnames_tissueassign.txt > $CAFILE > used to create the $EAINFILE in step 6 below
## 4. Using the file from step 2 above, extract the species-specific TG ensemblID > col1 and TF motifID > col2 (like those used in file created for step 1) > equivalent to ~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/results/refined_motifs_Ab/Pn_motifnames_regnet.txt > $NET
## 5. Using the files from step 4, run PrepEAfiles_local.sh to prepare the filtered/refined and renamed motif instances as inputs for the enrichment program
## 6. Run the enrichment analysis - test the enrichment of TF motifs in open-chromatin regions overlapping gene promoter regions against a background of open-chromatin peaks in the whole genome

# ============================================
# Step 1: Prepare a motif ID and TF name list
# ============================================

## 1. Prepare a motif ID and TF name list with col1-motifID and col2-motif_name for ALL input motifs used for calling TF footprints (all found here: /hpc-home/mehtat/rgtdata/motifs/*.mtf) > equivalent to ~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/data/JASPAR_Vertebrates_2018_motifname.txt > $MOTIFNAME

mkdir -p /Users/mehtat/Documents/TGAC/Technology_Development/ATAC-Seq/ATAC_Bioinformatics/TFfootprints
cd /Users/mehtat/Documents/TGAC/Technology_Development/ATAC-Seq/ATAC_Bioinformatics/TFfootprints

# The motif ID and TF name lists are here: /hpc-home/mehtat/rgtdata/motifs/*.mtf
# copied the above over

# footprint *.mbps outputs are (so you know what motifID format to follow)
# LG1	39234	39245	MA1628.1.Zic1::Zic2	14.26072286231245	-
# LG1	39238	39246	MAFG_HUMAN.H11MO.1.A	11.03991633366865	+
# LG1	39238	39246	MAFG_MOUSE.H11MO.1.A	11.03991633366865	+
# LG1	39318	39324	MA0151.1.Arid3a	10.676982226424872	-
# LG1	39311	39322	MA0468.1.DUX4	12.679446248088677	-


# An example of a motif ID and TF name list is in: /Users/mehtat/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/data/JASPAR_Vertebrates_2018_motifname.txt
# Format:
# MA0004.1	ARNT
# MA0006.1	AHR::ARNT
# MA0019.1	DDIT3::CEBPA

cat *.mtf | awk '{print $2,toupper($4)}' OFS='\t' | sort -u -k1,1 > allTFmotifs_motifname.txt # DONE - 4,408 motifs across cichlid and vertebrates (input set)

cut -f2 allTFmotifs_motifname.txt | awk '!visited[$0]++' | wc -l # 1442 - number of UNIQUE motifs


# ============================================
# Step 2: Filter for TF footprints in promoter regions and annotate
# ============================================

## 2. Using the *_ATAC_mpbs.bed output from TF footprinting, filter all TF footprints for those falling in gene promoter regions and associate with gene (based on gene promoter)

# rgt-hint footprinting output you want is a *_mpbs.bed file
# bed file output: col1 - chr, col2 - start, col3 - end, col4 - motifID, col5 - bit-score of the motif matching, col6 - strand
# footprinting output bed files here: /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/3.TFfprint_SignalTrack/{mz,pn,ab,nb,on}_fp/*_{B,E,L,T}_ATAC/*_ATAC_mpbs.bed
# process these on HPC and then copy over local

mkdir -p /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files
cd /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files


# copied to local folder here: /Users/mehtat/Documents/TGAC/Technology_Development/ATAC-Seq/ATAC_Bioinformatics/TFfootprints
# Haven't obvs included On Gill here ..

nano 2.TFfp_annot.sh

#!/bin/bash -e
#SBATCH -p ei-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 36000 # memory pool for all cores
#SBATCH -t 0-05:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

fpdir=/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/3.TFfprint_SignalTrack
annotdir=/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/4b.peak_phylop/final_annot
mzannotbed=${annotdir}/Metriaclima_zebra.M_zebra_UMD2a.annot.bed
pnannotbed=${annotdir}/Pundamilia_nyererei.PunNye1.0.annot.bed
abannotbed=${annotdir}/Astatotilapia_burtoni.AstBur1.0.annot.bed
nbannotbed=${annotdir}/Neolamprologus_brichardi.NeoBri1.0.annot.bed
onannotbed=${annotdir}/Oreochromis_niloticus.O_niloticus_UMD_NMBU.annot.bed

source bedtools-2.30.0

for i in ${fpdir}/{mz,pn,ab,nb,on}_fp/*_{B,E,L,T}_ATAC/*_ATAC_mpbs.bed; do
  id=$(echo $(basename "${i}" ) | awk -F'_' '{print $1}' | sed 's/[0-9]*//g' | sed 's|Pnm|Pn|g')
	file=$(echo $(basename "${i}") | sed 's|_mpbs.bed|_mpbs.notab.bed|g')
  if [[ ${id} =~ "Mz" ]]; then
		sed 's/\t$//' $i > $file # mpbs files have tabs at the end of each line so remove
    echo -e ${id}'\t'"True"
    bedtools intersect -a $mzannotbed -b ${file} -wb |
		awk '$4=="5kb_gene_promoter"' |
		awk '{print $11,$12,$13,$14,$15,$16,$7,$10,$4,$6}' OFS='\t' > $(basename "${file}" .notab.bed).geneprom.bed
		rm $file
  fi
	if [[ ${id} =~ "Pn" ]]; then
		sed 's/\t$//' $i > $file # mpbs files have tabs at the end of each line so remove
    echo -e ${id}'\t'"True"
    bedtools intersect -a $pnannotbed -b ${file} -wb |
		awk '$4=="5kb_gene_promoter"' |
		awk '{print $11,$12,$13,$14,$15,$16,$7,$10,$4,$6}' OFS='\t' > $(basename "${file}" .notab.bed).geneprom.bed
		rm $file
  fi
	if [[ ${id} =~ "Ab" ]]; then
		sed 's/\t$//' $i > $file # mpbs files have tabs at the end of each line so remove
    echo -e ${id}'\t'"True"
    bedtools intersect -a $abannotbed -b ${file} -wb |
		awk '$4=="5kb_gene_promoter"' |
		awk '{print $11,$12,$13,$14,$15,$16,$7,$10,$4,$6}' OFS='\t' > $(basename "${file}" .notab.bed).geneprom.bed
		rm $file
  fi
	if [[ ${id} =~ "Nb" ]]; then
		sed 's/\t$//' $i > $file # mpbs files have tabs at the end of each line so remove
    echo -e ${id}'\t'"True"
    bedtools intersect -a $nbannotbed -b ${file} -wb |
		awk '$4=="5kb_gene_promoter"' |
		awk '{print $11,$12,$13,$14,$15,$16,$7,$10,$4,$6}' OFS='\t' > $(basename "${file}" .notab.bed).geneprom.bed
		rm $file
  fi
	if [[ ${id} =~ "On" ]]; then
		sed 's/\t$//' $i > $file # mpbs files have tabs at the end of each line so remove
    echo -e ${id}'\t'"True"
    bedtools intersect -a $onannotbed -b ${file} -wb |
		awk '$4=="5kb_gene_promoter"' |
		awk '{print $11,$12,$13,$14,$15,$16,$7,$10,$4,$6}' OFS='\t' > $(basename "${file}" .notab.bed).geneprom.bed
		rm $file
  fi
done

# run the above
sbatch 2.TFfp_annot.sh


# ============================================
# Step 3: Species-specific TFfp annot with col1-ensemblID and col2-tissue
# ============================================

## 3. Create a two column file with col1-ensemblID and col2-tissue for the above but for each species > equivalent to ~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/results/k10_arb_output/prediction/Puny_speciesspecnames_tissueassign.txt > $CAFILE > used to create the $EAINFILE in step 6 below

cd /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files

for i in *_ATAC_mpbs.geneprom.bed; do
	sp=$(echo $i | awk -F'_' '{print $1}' | sed 's/[0-9]*//g' | sed 's|Pnm|Pn|g')
	tissue=$(echo ${i} | awk -F'_' '{print $2}')
	# echo -e $sp'\t'$tissue
	awk -v var=$tissue -v var2=$sp '{print $7,var,var2}' OFS='\t' ${i} >> speciesspecnames_tissueassign.txt
done

for i in Mz Pn Ab Nb On; do
	awk -v var=$i '$3==var' speciesspecnames_tissueassign.txt | cut -f1,2 > ${i}_speciesspecnames_tissueassign.txt
done # DONE

# ============================================
# Step 4: Species-specific TFfp annot with TG ensemblID > col1 and TF motifID > col2
# ============================================

## 4. Using the file from step 2 above, extract the species-specific TG ensemblID > col1 and TF motifID name > col2 (like those used in file created for step 1) > equivalent to ~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs/motif_discovery/promoters_JASPAR2018/cichlids_motif_enrichments/results/refined_motifs_Ab/Pn_motifnames_regnet.txt > $NET

cd /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files

for i in *_ATAC_mpbs.geneprom.bed; do
	sp=$(echo $i | awk -F'_' '{print $1}' | sed 's/[0-9]*//g' | sed 's|Pnm|Pn|g')
	# tissue=$(echo ${i} | awk -F'_' '{print $2}')
	# echo -e $sp'\t'$tissue
	awk -v var=$sp '{print $7,$4,var}' OFS='\t' ${i} >> all_motifnames_regnet.tmp
done

awk 'BEGIN{OFS="\t"}NR==FNR{a[$1]=$2;next}{if(a[$2]){print $0,a[$2];}else{print $0,"NULL";}}' allTFmotifs_motifname.txt all_motifnames_regnet.tmp | awk '{print $1,$4,$2,$3}' OFS='\t' > all_motifnames_regnet.txt

for i in Mz Pn Ab Nb On; do
	awk -v var=$i '$4==var' all_motifnames_regnet.txt | cut -f1,2 > ${i}_motifnames_regnet.txt
done # DONE

rm all_motifnames_regnet.tmp # DONE

# ============================================
# Step 5: Prepare the filtered/refined and renamed motif instances
# ============================================

cd /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files


wd=/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files

for SPC in Ab Pn On Mz Nb; do
	export NET=${wd}/${SPC}_motifnames_regnet.txt
	#head ${NET}
	cat ${NET} | grep -v NULL | awk 'NF==2{printf("%s\t%s\t1\n",$1,$2)}' > ${wd}/${SPC}_gotermap.txt
	# cat ${NET} | grep -v NULL  | cut -f2 | sort | uniq -c | head
	cat ${NET} | grep -v NULL  | cut -f2 | sort | uniq -c | awk '{printf("%s\t%s\t1\n",$2,$1)}' > ${wd}/${SPC}_genecnt.txt
done

# for local unix
# rename 's/Ab_/Asbu_/g' ${wd}/Ab_*
# rename 's/Mz_/Meze_/g' ${wd}/Mz_*
# rename 's/Pn_/Puny_/g' ${wd}/Pn_*
# rename 's/On_/Orni_/g' ${wd}/On_*
# rename 's/Nb_/Nebr_/g' ${wd}/Nb_*

# for the cluster
rename Ab_ Asbu_ ${wd}/Ab_*
rename Mz_ Meze_ ${wd}/Mz_*
rename Pn_ Puny_ ${wd}/Pn_*
rename On_ Orni_ ${wd}/On_*
rename Nb_ Nebr_ ${wd}/Nb_*

# Final motif annotation files are in /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files

# ============================================
# Step 6: Run the enrichment analysis
# ============================================

# The tool for applying the enrichment analysis is /tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer, and this is implemented with the steps in the EnrichmentAnalysis.sh script

# All files need to arranged accordingly as per below
mkdir -p /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files/enrichment/{input1,input2}

cd /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files/enrichment

mv /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files/*_{gotermap,genecnt}.txt /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files/enrichment/input1 # these files have to be in a separate folder
mv /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files/*_speciesspecnames_tissueassign.txt /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files/enrichment/input2 # these files have to be in a separate folder
mv /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files/{Asbu,Meze,Puny,Nebr,Orni}_motifnames_regnet.txt /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files/enrichment/input2


nano EnrichmentAnalysis.sh

#!/bin/bash -e
#SBATCH -p ei-medium # partition (queue)
#SBATCH -N 1 # number of nodes
#SBATCH -n 1 # number of tasks
#SBATCH --mem 48000 # memory pool for all cores
#SBATCH -t 0-05:59 # time (D-HH:MM)
#SBATCH -o slurm.%N.%j.out # STDOUT
#SBATCH -e slurm.%N.%j.err # STDERR
#SBATCH --mail-type=ALL # notifications for job done & fail
#SBATCH --mail-user=Tarang.Mehta@earlham.ac.uk # send-to address

wd=/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files

ml gcc
ml zlib

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/tgac/software/testing/enrichAnalyzer/0.1/x86_64/lib:/tgac/software/testing/enrichAnalyzer/0.1/x86_64/lib2/
export PYGOINPUT=/ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files/enrichment/ReworkClusterAssign.py
export GOEXE=/tgac/software/testing/enrichAnalyzer/0.1/src/enrichAnalyzer

#changed the following lines manually
export Species_List_All="Meze Puny Asbu Nebr Orni"
export OUTDIR=${wd}/enrichment/input1
export OUTDIR2=${wd}/enrichment/input2
export EAOUTDIR=${wd}/enrichment/output

echo -e "Starting GO analysis"
mkdir -p ${EAOUTDIR}

for SN in Meze Puny Asbu Nebr Orni; do
	CAFILE=${OUTDIR2}/${SN}_speciesspecnames_tissueassign.txt
  EAINFILE=${EAOUTDIR}/${SN}_MOTIFS_INPUT.txt
  python3 ${PYGOINPUT} ${CAFILE} > ${EAINFILE}
done

for SN in Meze Puny Asbu Nebr Orni
do
	export EAINFILE=${EAOUTDIR}/${SN}_MOTIFS_INPUT.txt
	export EAOUTPUT=${EAOUTDIR}/${SN}_MOTIFS_OUTPUT
	export CAFILE=${OUTDIR2}/${SN}_tissueassign.txt
	if [[ -e ${OUTDIR2}/${SN}_speciesspecnames_tissueassign.txt ]]
	then
		CAFILE=$(echo ${OUTDIR2}/${SN}_speciesspecnames_tissueassign.txt | cut -f1 | sort -u)
	fi
	EADAT=${OUTDIR}/${SN}_
        ${GOEXE} ${EAINFILE} ${CAFILE} ${EADAT} 1 ${EAOUTPUT} persg
done

# run the above
sbatch EnrichmentAnalysis.sh


## Column headers for the output are:

#1. Tissue
#2. Enriched TF motif for tissue col1 in open-chromatin gene promoter regions
#3. Hypergeometric test p-value
#4. Corrected p-value (Benjamini-Hochberg procedure, FDR)
#5. Total number of genes with any term that intersect the total numbers associated with the tissue
#6. Total number of genes with the term in column 2
#7. Total number of target gene promoter regions for enrichment
#8. Total number of target gene promoter regions that have the term in col2
#9. Fold enrichment (must be >1, otherwise things wouldn't be significant). But it is essentially the ratio of observed over expected fraction of motif hits.
#10. Target gene promoters (ensemblID) that have enrichment for the term in col 2


# ============================================
# Step 7: Prepare files for enrichment plots
# ============================================

# TF enrichment outputs here:
# /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files/enrichment/output/*_MOTIFS_OUTPUT_details.txt

# 1. remove 10th col, extend naming for tissue column e.g. ClusterB > Brain, and add species column
# 2. plot
	# adj. p-val (FDR) < 0.05 (x-axis);
	# enriched TF motif (y-axis);
	# point size as number of open-chromatin gene promoters for the enriched term (col9)
	# point colour as fold enrichment
	# species as main facet and tissues as sub-facets

cd /ei/projects/9/9e238063-c905-4076-a975-f7c7f85dbd56/scratch/ATACseq/3.run2/XX.manuscriptfiles/2.fig2files/enrichment/output

for i in *_MOTIFS_OUTPUT_details.txt; do
	sp=$(echo ${i} | awk -F'_' '{print $1}' | sed 's|Meze|M. zebra|g' | sed 's|Puny|P. nyererei|g' | sed 's|Asbu|A. burtoni|g' | sed 's|Nebr|N. brichardi|g' | sed 's|Orni|O. niloticus|g')
	awk -v var="$sp" '{print var,$1,$2,$3,$4,$5,$6,$7,$8,$9}' OFS='\t' ${i} |
	sed 's|ClusterB|Brain|g' | sed 's|ClusterE|Eye|g' | sed 's|ClusterL|Liver|g' | sed 's|ClusterT|Testis|g' > $(basename "${i}" .txt).simp.txt
done

# above added local here: ~/Documents/TGAC/Technology_Development/ATAC-Seq/ATAC_Bioinformatics/TFfootprints/enrichment/*_MOTIFS_OUTPUT_details.simp.txt


# how many TFs (of the 1,442 input) are enriched in each species
for i in *_MOTIFS_OUTPUT_details.simp.txt; do
  count=$(cut -f3 ${i} | awk '!visited[$0]++' | wc -l)
  echo -e ${i}'\t'$count
done
# Asbu_MOTIFS_OUTPUT_details.simp.txt	824
# Meze_MOTIFS_OUTPUT_details.simp.txt	824
# Nebr_MOTIFS_OUTPUT_details.simp.txt	824
# Orni_MOTIFS_OUTPUT_details.simp.txt	824
# Puny_MOTIFS_OUTPUT_details.simp.txt	824

# 3. The above has too many motifs for a main figure - insert the motif enrichment results as a table in supp tables > Supp Table S10 # DONE
# filter the motif enrichment for top 5 enriched motifs (based on FE) in each species tissue

for i in Brain Eye Liver Testis; do
  for j in *_MOTIFS_OUTPUT_details.simp.txt; do
    grep ${i} ${j} | sort -t $'\t' -k10,10rn | head -5 > $(basename "$j" .txt).top5_${i}.txt
  done
done

for i in Asbu Meze Puny Nebr Orni; do
  cat ${i}_MOTIFS_OUTPUT_details.simp.top5_* > ${i}_MOTIFS_OUTPUT_details.simp.top5.txt
done

cat *_MOTIFS_OUTPUT_details.simp.top5.txt > All_MOTIFS_OUTPUT_details.simp.top5.txt # use this file for plotting top 5 TF enrichment for each species and tissue

# 4. The above does not highlight particularly interesting examples of motif enrichment
# Instead, show the TF motif enrichment of TFs that have tissue-specific functions: can look at compendium and other papers

# Brain
# Cichlid brains evolve diversity via subtle modification of conserved gene regulatory networks (ref: https://www.nature.com/articles/ncomms2753; https://www.pnas.org/doi/full/10.1073/pnas.1000395107)
# foxg1 - differentially expressed in rock- vs. sand-dwelling Malawi cichlids - https://www.nature.com/articles/ncomms2753
# ap2a - ectodermal migration during neural tube closure and cell fate specification - https://bmcdevbiol.biomedcentral.com/articles/10.1186/s12861-017-0146-0
# egr1 - GnRH neurons in newly dominant A. burtoni males show increased expression of a marker of neural activation (the transcription factor EGR1) within 20 min of these males expressing dominance-related behaviors. This ascent to dominance is followed by rapid increases in GnRH neuron size and testes size and activity. https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0030363

# Eye:
# CRX - https://bmcecolevol.biomedcentral.com/articles/10.1186/1471-2148-11-120
# NR2C2 - your genome bio paper
# Rx1 - https://academic.oup.com/mbe/article/31/9/2297/2925686
# RARα/β/γ and RXRα/β/γ - Browman HI, Hawryshyn CW. Retinoic acid modulates retinal development in the juveniles of a teleost fish. J Exp Biol. 1994;193:191–207.

# Liver:
# hnf4a - https://www.nature.com/articles/s41467-021-26166-2
# foxk1 - https://www.nature.com/articles/s41467-021-26166-2; https://www.nature.com/articles/nature03649
# foxa1 - https://www.nature.com/articles/nature03649
# foxa2 - https://www.nature.com/articles/nature03649

# Testis:
# egr1 - GnRH neurons in newly dominant A. burtoni males show increased expression of a marker of neural activation (the transcription factor EGR1) within 20 min of these males expressing dominance-related behaviors. This ascent to dominance is followed by rapid increases in GnRH neuron size and testes size and activity. https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0030363
# arα/β - control of testis mass; https://pubmed.ncbi.nlm.nih.gov/34643776/
# dmrt1 - male sex differentiation factor; https://pubmed.ncbi.nlm.nih.gov/30415011/

cd ~/Documents/TGAC/Technology_Development/ATAC-Seq/ATAC_Bioinformatics/TFfootprints/enrichment

for i in foxg1 ap2a egr1 crx nr2c2 rx1 hnf4a foxk1 foxa1 ar dmrt1; do
  for j in *_MOTIFS_OUTPUT_details.simp.txt; do
    grep -wiF ${i} ${j} >> All_MOTIFS_OUTPUT_details.simp.cand.txt
  done
done
