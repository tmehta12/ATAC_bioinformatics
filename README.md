# ATAC_bioinformatics
Multi tissue and multi species ATAC bioinformatics analysis pipeline

### A. Merging reads sequenced across multiple lanes
1. Create a topmost directory separately, placing this script in there and running from that directory e.g. /tgac/workarea/group-vh/Tarang/ATACseq/2.run2
2. sbatch ATAC_Bioinf_pipeline_v2a.sh

### B. Trimming and alignment of gDNA reads
Script usage: ./ATAC_Bioinf_pipeline_v2b.sh -s "spID" -g "spG" -f "gFA" -m "mtID" -u "Usr" -a "annot"
e.g. ./ATAC_Bioinf_pipeline_v2b.sh -s Mz1_L_ATAC -g M_zebra_UMD1 -f /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -m KT221043 -u mehtat
Note: Script is adapted for SBATCH usage

This pipeline is ran as species-specific, and contains the following components:
1. Trim adaptors - trimgalore
2. Read alignment - bowtie2 > samtools sorted bam

### C. Trimming and alignment of ATAC reads, and then filtering, calling peaks, and annotation
Script usage: ./ATAC_Bioinf_pipeline_v2b.sh -s "spID" -g "spG" -f "gFA" -m "mtID" -u "Usr" -a "annot"
e.g. ./ATAC_Bioinf_pipeline_v2b.sh -s Mz1_L_ATAC -g M_zebra_UMD1 -f /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -m KT221043 -u mehtat
Note: Script is adapted for SBATCH usage

Place the following files in $scripts - to be created prior to running
- As used in 'ATAC_Bioinf_pipeline_v2a.sh': a 2-column space-delimited table where col1='R1/R2 filename's col2='desired species renamed filename: species_tissue_experiment e.g. Mz_L_ATAC/gDNA'
- Scripts:
  - ATAC_Bioinf_pipeline_v2b_part3a.py
  - ATAC_Bioinf_pipeline_v2b_part3b.R
  - ATAC_Bioinf_pipeline_v2b_part5bD.py
- Run as an sbatch script with 8Gb memory and >15 day runtime - will spawn off other jobs

~ This pipeline is ran as species-specific, and contains the following components:

1. Trim adaptors - trimgalore
2. Read alignment - bowtie2 > samtools sorted bam
3. Remove mitochondrial mapped - samtools > new bam
  - 3a. Fragment/insert size distribution
4. Post alignment filtering
  - 4a. Filter reads (Sort, Map and remove duplicates, Remove reads unmapped, not primary alignment, reads failing platform, duplicates) - samtools > final bam
5. ATAC peak-calling
  - 5a. Convert PE Bam to tagalign (start/end positions of each read) - bedtools
  - 5b. TSS enrichment - plot
  - 5c. TN5 shifting of tagaligns - shift reads +4 bp for the +strand and -5 bp for the -strand
  - 5d. peak calling - macs2 NOTE: consider running another peak-calling program and take the intersection. Also, for analysis only consider open-chromatin so filter based on that?
6. IDR on all pairs of replicates, self-pseudoreplicates and pooled pseudoreplicates - IDR is optional. The IDR peaks are a subset of the naive overlap peaks that pass a specific IDR threshold of 10%.
  - 6a. IDR of true replicates
  - 6b. Compute Fraction of Reads in Peaks (FRiP)
7. Create signal tracks - bedtools
8. Annotation:
  - 8a. TSS enrichment
  - 8b. Fraction of Reads in annotated regions
