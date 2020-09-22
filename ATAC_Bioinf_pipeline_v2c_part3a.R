# ----------------------------------------------------------------------------
# title: Running ATACseqQC
# author: Tarang K. Mehta (EI)
# date: July 2020
# Description: Running ATACseqQC for diagnostic plots of TSS enrichment
# How to run: Rscript ATAC_Bioinf_pipeline_v2c_part3a.R -i ${i} -o "$(echo ${i} | sed 's/.meme/.tmp.pwm/g')
# ----------------------------------------------------------------------------

## Quick start - install ATACseqQC and other packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install(version = "3.11")
# library(BiocManager)
# BiocManager::install(c("RMySQL","rtracklayer","GenomicFeatures","GLAD","gsl","ensembldb","GenomicRanges","MotIV","motifStack","ATACseqQC","ChIPpeakAnno", "MotifDb", "GenomicAlignments","Rsamtools","BSgenome","Biostrings","ggplot2","optparse"))

## load the libraries
# library("RMySQL")
library("rtracklayer")
library("GenomicFeatures")
library("GLAD")
library("gsl")
library("ensembldb")
library("GenomicRanges")
library("MotIV")
library("motifStack")
library("ATACseqQC")
library("ChIPpeakAnno")
library("MotifDb")
library("GenomicAlignments")
library("Rsamtools")
library("BSgenome")
library("Biostrings")
library("ggplot2")
library("optparse")

## set the options list for command line input

option_list = list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character',
              help="input *.bam file (nochrM-nodup-filtered-sorted; non-shifted!!)"),
  make_option(c("-g", "--gtf"), action="store", default=NA, type='character',
              help="input *.gtf file"),
  make_option(c("-s", "--shiftbam"), action="store", default=NA, type='character',
              help="path to output shifted and split BAMs - make and name path folder according to sample and tissue"),
  make_option(c("-p", "--ptp"), action="store", default=NA, type='character',
              help="output *.tiff filename for Promoter-Transcript score plot"),
  make_option(c("-n", "--nfrp"), action="store", default=NA, type='character',
              help="output *.tiff filename for NFR score plot"),
  make_option(c("-b", "--bsgenome"), action="store", default=NA, type='character',
              help="BSgenome package dir as built in shell script e.g. BSgenome.Abur.Ensembl.AstBur1.0"),
  make_option(c("-t", "--tssscore"), action="store", default=NA, type='character',
              help="output *.txt filename for TSS enrichment score summary"),
  make_option(c("-c", "--cpp"), action="store", default=NA, type='character',
              help="output *.tiff filename for cumulative percentage plot"),
  make_option(c("-h", "--hmp"), action="store", default=NA, type='character',
              help="output *.tiff filename for log-transformed signal around TSSs"),
  make_option(c("-r", "--rsp"), action="store", default=NA, type='character',
              help="output *.tiff filename for rescaled signal around TSSs")
)
opt = parse_args(OptionParser(option_list=option_list))

# setwd("/Users/mehtat/github/ATAC_bioinformatics/test_data/")

## input the bamFile from command line (1)
bamfile <- opt$i
# bamfile <- ("Ab5_L_ATAC.nochrM.nodup.filt.sorted.JH425323.1.bam")
bamfile.labels <- gsub(".bam", "", basename(bamfile))

## input genome/gtf paths
gtf_gr <- import(opt$g)
# gtf_gr <- import("Haplochromis_burtoni.AstBur1.0.100.gtf")
gtf_txdb <- makeTxDbFromGRanges(gtf_gr)

## Estimate the library complexity
estimateLibComplexity(readsDupFreq(bamfile))

## Nucleosome positioning

### Adjust the read start sites

# Tn5 transposase has been shown to bind as a dimer and inserts two adaptors into accessible DNA locations separated by 9 bp.
# Therefore, for downstream analysis, such as footprinting, all reads in input bamfile need to be shifted.
# The function `shiftGAlignmentsList` can be used to shift the reads.
# By default, all reads aligning to the positive strand are offset by +4bp,
# and all reads aligning to the negative strand are offset by -5bp.

# The adjusted reads will be written into a new bamfile for TSS enrichment and footprinting.

## bamfile tags to be read in
possibleTag <- combn(LETTERS, 2)
possibleTag <- c(paste0(possibleTag[1, ], possibleTag[2, ]),
                 paste0(possibleTag[2, ], possibleTag[1, ]))

bamTop100 <- scanBam(BamFile(bamfile, yieldSize = 100),
                     param = ScanBamParam(tag=possibleTag))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)==100]
# tags

## files will be output into outPath
# outPath <- ("/Users/mehtat/github/ATAC_bioinformatics/test_data/shiftedBAMout")
outPath <- opt$s
dir.create(outPath)
## shift the coordinates of 5'ends of alignments in the bam file
gal <- readBamFile(bamfile, tag=tags, asMates=TRUE, bigFile=TRUE)
shiftedbamfile.labels <- gsub(".bam", ".shifted.bam", basename(bamfile))
shiftedBamfile <- file.path(outPath, shiftedbamfile.labels)
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)


### Promoter/Transcript body (PT) score
# PT score is calculated as the coverage of promoter divided by the coverage of its transcript body.
# PT score will show if the signal is enriched in promoters.

txs <- transcripts(gtf_txdb)

pt <- PTscore(gal1, txs)

tiff(opt$p, units="in", width=5, height=5, res=100)
# tiff("xx.tiff", units="in", width=5, height=5, res=100)
plot(pt$log2meanCoverage, pt$PT_score,
     xlab="log2 mean coverage",
     ylab="Promoter vs Transcript",
     main="Promoter/Transcript body (PT) score")
dev.off()

### Nucleosome Free Regions (NFR) score

# NFR score is a ratio between cut signal adjacent to TSS and that flanking the corresponding TSS.
# Each TSS window of 400 bp is first divided into 3 sub-regions: the most upstream 150 bp (n1), the most downstream of 150 bp (n2), and the middle 100 bp (nf).
# Then the number of fragments with 5' ends overlapping each region are calculated for each TSS.
# The NFR score for each TSS is calculated as NFR-score = log2(nf) - log2((n1+n2)/2).
# A plot can be generated with the NFR scores as Y-axis and the average signals of 400 bp window as X-axis, very like a MA plot for gene expression data.


nfr <- NFRscore(gal1, txs)

tiff(opt$n, units="in", width=5, height=5, res=100)
# tiff("yy.tiff", units="in", width=5, height=5, res=100)
plot(nfr$log2meanCoverage, nfr$NFR_score,
     xlab="log2 mean coverage",
     ylab="Nucleosome Free Regions score",
     main="NFRscore for 200bp flanking TSSs",
     xlim=c(-10, 0), ylim=c(-5, 5))
dev.off()


### Transcription Start Site (TSS) Enrichment Score

# TSS enrichment score is a ratio between aggregate distribution of reads centered on TSSs and that flanking
# the corresponding TSSs. TSS score = the depth of TSS (1000 bp each side) / the depth of end flanks (100bp each end).
# TSS enrichment score is calculated according to the definition at [https://www.encodeproject.org/data-standards/terms/#enrichment](https://www.encodeproject.org/data-standards/terms/#enrichment).
# Transcription start site (TSS) enrichment values are dependent on the reference files used; cutoff values for high quality data are listed below and in the following table from [https://www.encodeproject.org/atac-seq/](https://www.encodeproject.org/atac-seq/).
# # GRCh38 Refseq TSS annotation
# #     below 5: Concerning
# #     5 - 7: Acceptable
# #     Above 7: Ideal

tsse <- TSSEscore(gal1, txs)
tsse$TSS.enrichment.score[sapply(tsse$TSS.enrichment.score, is.infinite)] <- NA
tssenrichsum <- summary(tsse$TSS.enrichment.score)
tssenrichsum2 <- data.frame(lapply(tssenrichsum, function(x) t(data.frame(x))))
sampleID <- gsub("ATAC.*", "ATAC", basename(bamfile))
rownames(tssenrichsum2)[rownames(tssenrichsum2) == "x"] <- sampleID
tssenrichsum2$cutoff <- cut(tssenrichsum2$Mean, c(-Inf,5,7,Inf), c("concerning", "acceptable", "ideal"))
tssenrichsum2 <- cbind(rownames(tssenrichsum2), data.frame(tssenrichsum2, row.names=NULL))
colnames(tssenrichsum2)[1] <- "sample"
# write.table(tssenrichsum2, file = "Ab5_L_ATAC_JH425323.1_TSSenrichmentscore.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
write.table(tssenrichsum2, file = opt$t, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

### Split reads

# The shifted reads will be split into different bins, namely
# nucleosome free, mononucleosome, dinucleosome, and trinucleosome.
# Shifted reads that do not fit into any of the above bins will
# be discarded. Splitting reads is a time-consuming step
# because we are using random forest to classify the fragments
# based on fragment length, GC content and conservation scores
# [@chen2013danpos].

# By default, we assign the top 10% of short reads (reads below 100_bp)
# as nucleosome-free regions and the top 10% of intermediate length reads
# as (reads between 180 and 247 bp) mononucleosome.
# This serves as the training set to classify the rest of the fragments
# using random forest. The number of the tree will be set to 2 times
# of square root of the length of the training set.

## Since we need to use custom genomes that are different to those available, forge a BSgenome package using bare sequences
# The seedfile, fasta sequences and BSgenome package build for this have been prepped in ATAC_Bioinf_pipeline_v2c.sh
# library(BSgenome.Abur.Ensembl.AstBur1.0)
library(opt$b)
genome <- Abur

## split the reads into NucleosomeFree, mononucleosome, dinucleosome and trinucleosome.
objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = outPath)

## save the binned alignments into bam files.
null <- writeListOfGAlignments(objs, outPath)
# dir(outPath) # list the files generated


### Heatmap and coverage curve for nucleosome positions

# By averaging the signal across all active TSSs, we should observe that
# nucleosome-free fragments are enriched at the TSSs,
# whereas the nucleosome-bound fragments should be enriched both upstream
# and downstream of the active TSSs and display characteristic phasing of upstream and
# downstream nucleosomes. Because ATAC-seq reads are concentrated at regions of
# open chromatin, should expect to see a strong nucleosome signal at the +1
# nucleosome, but the signal decreases at the +2, +3 and +4 nucleosomes.

bamfiles <- file.path(outPath,
                     c("NucleosomeFree.bam",
                     "mononucleosome.bam",
                     "dinucleosome.bam",
                     "trinucleosome.bam"))

## Plot the cumulative percentage of tag allocation in nucleosome-free and mononucleosome bam files.
tiff(opt$c, units="in", width=5, height=5, res=100)
# tiff("zz.tiff", units="in", width=5, height=5, res=100)
cumulativePercentage(bamfiles[1:2], as(seqinfo(Abur), "GRanges"))
dev.off()

## define the promoter and TSS regions
TSS <- promoters(txs, upstream=0, downstream=1)
TSS <- unique(TSS)

## estimate the library size for normalization
(librarySize <- estLibSize(bamfiles))

## calculate the signals around TSSs.
NTILE <- 101
dws <- ups <- 1010
# seqlev <- "JH425323.1" # only used for subsampling for a quick run; remove if running on whole BAM
sigs <- enrichedFragments(gal=objs[c("NucleosomeFree",
                                     "mononucleosome",
                                     "dinucleosome",
                                     "trinucleosome")],
                          TSS=TSS,
                          librarySize=librarySize,
                          # seqlev=seqlev, # this is only used for subsampling
                          TSS.filter=0.5,
                          n.tile = NTILE,
                          upstream = ups,
                          downstream = dws)

## log2 transformed signals
sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))

#plot heatmap
tiff(opt$h, units="in", width=5, height=7, res=100)
# tiff("aa.tiff", units="in", width=5, height=7, res=100)
featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                      zeroAt=.5, n.tile=NTILE)
dev.off()

## get signals normalized for nucleosome-free and nucleosome-bound regions.
out <- featureAlignedDistribution(sigs,
                                  reCenterPeaks(TSS, width=ups+dws),
                                  zeroAt=.5, n.tile=NTILE, type="l",
                                  ylab="Averaged coverage")

## rescale the nucleosome-free and nucleosome signals to 0~1
range01 <- function(x){(x-min(x))/(max(x)-min(x))}
out <- apply(out, 2, range01)

tiff(opt$r, units="in", width=8, height=6, res=100)
# tiff("bb.tiff", units="in", width=8, height=6, res=100)
matplot(out, type="l", xaxt="n",
        xlab="Position (bp)",
        ylab="Fraction of signal")
axis(1, at=seq(0, 100, by=10)+1,
     labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
dev.off()

# ## plot Footprints
#
# # ATAC-seq footprints infer factor occupancy genome-wide. The `factorFootprints`
# # function uses `matchPWM` to predict the binding sites using the input position
# # weight matrix (PWM).
# # Then it calculates and plots the accumulated coverage for those binding sites
# # to show the status of the occupancy genome-wide.
# # Unlike CENTIPEDE[@pique2011accurate], the footprints generated here
# # do not take the conservation (PhyloP) into consideration.
# # `factorFootprints` function could also accept the
# # binding sites as a GRanges object.
#
#
# ## footprints e.g. CTCF
# CTCF <- query(MotifDb, c("CTCF"))
# CTCF <- as.list(CTCF)
# # print(CTCF[[1]], digits=2)
# sigs <- factorFootprints(shiftedBamfile, pfm=CTCF[[1]],
#                          genome=genome,
#                          min.score="90%", seqlev=seqlev,
#                          upstream=100, downstream=100)
#
# featureAlignedHeatmap(sigs$signal,
#                       feature.gr=reCenterPeaks(sigs$bindingSites,
#                                                width=200+width(sigs$bindingSites[1])),
#                       annoMcols="score",
#                       sortBy="score",
#                       n.tile=ncol(sigs$signal[[1]]))
# sigs$spearman.correlation
# sigs$Profile.segmentation
#
#
# ### V-plot
# # V-plot is a plot to visualize fragment midpoint vs length for a given transcription factors.
# vp <- vPlot(shiftedBamfile, pfm=CTCF[[1]],
#             genome=genome, min.score="90%", seqlev=seqlev,
#             upstream=200, downstream=200,
#             ylim=c(30, 250), bandwidth=c(2, 1))
# distanceDyad(vp, pch=20, cex=.5)
#
# # Plot correlations for multiple samples
# path <- system.file("extdata", package="ATACseqQC", mustWork=TRUE)
# bamfiles <- dir(path, "*.bam$", full.name=TRUE)
# gals <- lapply(bamfiles, function(bamfile){
#                readBamFile(bamFile=bamfile, tag=character(0),
#                           which=GRanges("chr1", IRanges(1, 1e6)),
#                           asMates=FALSE)
#          })
# library(TxDb.Hsapiens.UCSC.hg19.knownGene)
# txs <- transcripts(TxDb.Hsapiens.UCSC.hg19.knownGene)
# library(GenomicAlignments)
# plotCorrelation(GAlignmentsList(gals), txs, seqlev="chr1")
