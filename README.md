# ATAC_bioinformatics
Multi tissue and multi species ATAC bioinformatics analysis pipeline

### A. Merging reads sequenced across multiple lanes
1. Create a topmost directory separately, placing this script in there and running from that directory e.g. /tgac/workarea/group-vh/Tarang/ATACseq/2.run2
2. Run the script on slurm
      
        sbatch ATAC_Bioinf_pipeline_v2a.sh

### B. Trimming and alignment of gDNA reads
1. Place the script 'ATAC_Bioinf_pipeline_v2b_gDNA.sh' in each gDNA ($spID) folder only e.g. Mz_L_gDNA
2. Script usage: 

        ./ATAC_Bioinf_pipeline_v2b_gDNA.sh -s "spID" -g "spG" -f "gFA"

Example input:

    ./ATAC_Bioinf_pipeline_v2b_gDNA.sh -s Mz1_L_gDNA -g M_zebra_UMD1 -f /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa
    
Note: Script is adapted for SBATCH usage

This pipeline is ran as species-specific, from each folder, and contains the following components:
1. Trim adaptors - trimgalore
2. Read alignment - bowtie2 > samtools sorted bam

### C. Trimming and alignment of ATAC reads, and then filtering, calling peaks, and annotation
1. Place the script 'ATAC_Bioinf_pipeline_v2b.sh' in each ATAC ($spID) folder only e.g. Mz_L_ATAC
2. Place this script and the following files in $WD (created in first script)

        - As used in './ATAC_Bioinf_pipeline_v2a.sh': a 2-column space-delimited table where col1='R1/R2 filename's col2='desired species renamed filename: species_tissue_experiment e.g. Mz_L_ATAC/gDNA'
        - Scripts:
          ATAC_Bioinf_pipeline_v2b_part3a.py
          ATAC_Bioinf_pipeline_v2b_part3b.R
        - Run as an sbatch script with 8Gb memory and ~6 days runtime - will spawn off other jobs

2. Script usage: 

        ./ATAC_Bioinf_pipeline_v2b.sh -s "spID" -g "spG" -f "gFA" -m "mtID" -u "Usr" -a "annot"
        
Example input:

     ./ATAC_Bioinf_pipeline_v2b.sh -s Mz1_L_ATAC -g M_zebra_UMD1 -f /tgac/workarea/group-vh/Tarang/Reference_Genomes/cichlids/Assemblies_12092016/Maylandia_zebra/mze_ref_M_zebra_UMD1_chrUn.fa -m KT221043 -u mehtat

Note: Script is adapted for SBATCH usage

This pipeline is ran as species-specific, and contains the following components:

    1. Trim adaptors - trimgalore
    2. Read alignment - bowtie2 > samtools sorted bam
    3. Remove mitochondrial mapped - samtools > new bam
      3a. Fragment/insert size distribution
    4. Post alignment filtering
      4a. Filter reads (Sort, Map and remove duplicates, Remove reads unmapped, not primary alignment, reads failing platform, duplicates) - samtools > final bam
    5. ATAC peak-calling
      5a. Convert PE Bam to tagalign (start/end positions of each read) - bedtools
      5b. TN5 shifting of tagaligns - shift reads +4 bp for the +strand and -5 bp for the -strand
      5c. count-based peak calling using Poisson distribution - macs2
      5d. count-based peak calling using a program that considers properly paired, unpaired and secondary alignments (unlike MACS2) - Genrich
      5e. markov model based peak calling specific for ATAC-seq data - HMMRATAC (this is currently surpressed to run)
    6. bed to bigbed conversion for narrowpeaks - bedClip and bedToBigbed from ucsc_tools
