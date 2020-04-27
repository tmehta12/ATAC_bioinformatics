#!/usr/bin/env Rscript

## ---------------------------
##
## Script name: ATAC_Bioinf_pipeline_v2b_part3b.R
##
## Purpose of script: Plot fragment distribution from ATAC-seq data
##
## Author: Dr. Tarang K. Mehta
##
## Date Created: 11-03-2020
##
##
## ---------------------------
##
## Notes:
## This is ran from 'ATAC_Bioinf_pipeline_v2b.sh' at part 3
## Usage:
    ## source R-3.5.2
    ## R CMD BATCH --no-save --no-restore '--args input_BAM' ATAC_Bioinf_pipeline_v2b_part3b.R ATAC_Bioinf_pipeline_v2b_part3b.Rout
##
## ---------------------------
## Install relevant packages
# install.packages(knitr)
# install.packages(rmdformats)
# install.packages(dplyr)
# install.packages(DT)
# install.packages(tidyr)
# install.packages(ggplot2)
# install.packages(magrittr)
# install.packages(devtools)
## ---------------------------
## Load relevant packages
library("Rsamtools")
library("GenomicRanges")
library("ggplot2")
## ---------------------------
## ATAC-seq should represent a mix of fragment lengths corresponding to nucleosome free, mononucleosome and poly-nucleosome fractions.
## We can use insert lengths from across the genome to plot the distribution of all fragment lengths.

#' @title fragment size distribution
#' @description estimate the fragment size of bams
#' @param bamFiles A vector of characters indicates the file names of bams.
#' @param index The names of the index file of the 'BAM' file being processed;
#'        This is given without the '.bai' extension.
#' @param bamFiles.labels labels of the bam files, used for pdf file naming.
#' @param ylim numeric(2). ylim of the histogram.
#' @param logYlim numeric(2). ylim of log-transformed histogram for the insert.
#' @return Invisible fragment length distribution list.
#' @importFrom Rsamtools ScanBamParam scanBamFlag scanBam idxstatsBam
#' @importFrom graphics axis par plot
#' @import GenomicRanges
#' @export
#' @author Tarang K. Mehta
#' @examples
#' bamFiles <- dir(system.file("extdata", package="ATACseqQC"), "GL.*.bam$", full.names=TRUE)
#' bamFiles.labels <- sub(".bam", "", basename(bamFiles))
#' fragSizeDist(bamFiles, bamFiles.labels)

# setwd("/tgac/workarea/group-vh/Tarang/ATACseq/1.run1jan2017_twolanes/3a.post_alignment_filtering/Ab")

## allow for command line parameters
args = commandArgs(trailingOnly=TRUE)

# ## test if there is at least one argument: if not, return an error
# if (length(args)==0) {
#   stop("At least one argument must be supplied (input file).n", call.=FALSE)
# } else if (length(args)==1) {
#   # default output file
#   args[2] = "out.txt"
#   print(paste0("Input BAM is: ", args[1]))
# }

## Function for plotting the fragment size distribution
fragSizeDist <- function(bamFiles, bamFiles.labels, index=bamFiles, ylim=NULL,
                         logYlim=NULL){
  opar <- par(c("fig", "mar"))
  on.exit(par(opar))
  pe <- mapply(testPairedEndBam, bamFiles, index)
  if(any(!pe)){
    stop(paste(bamFiles[!pe], collapse = ", "),
         "is not Paired-End file.")
  }
  summaryFunction <- function(seqname, seqlength, bamFile, ind, ...) {
    param <-
      ScanBamParam(what=c('isize'),
                   which=GRanges(seqname, IRanges(1, seqlength)),
                   flag=scanBamFlag(isSecondaryAlignment = FALSE,
                                    isUnmappedQuery=FALSE,
                                    isNotPassingQualityControls = FALSE))
    table(abs(unlist(sapply(scanBam(bamFile, index=ind, ..., param=param),
                            `[[`, "isize"), use.names = FALSE)))
  }

  idxstats <- unique(do.call(rbind, mapply(function(.ele, .ind)
    idxstatsBam(.ele, index = .ind)[, c("seqnames", "seqlength")], bamFiles, index, SIMPLIFY=FALSE)))
  seqnames <- as.character(idxstats[, "seqnames"])
  seqlen <- as.numeric(idxstats[, "seqlength"])
  fragment.len <- mapply(function(bamFile, ind) summaryFunction(seqname=seqnames, seqlength=seqlen, bamFile, ind),
                         bamFiles, index, SIMPLIFY=FALSE)

  names(fragment.len) <- bamFiles.labels

  minor.ticks.axis <- function(ax,n=9,t.ratio=0.5,mn,mx,...){

    lims <- par("usr")
    lims <- if(ax %in% c(1,3)) lims[1:2] else lims[3:4]

    major.ticks <- pretty(lims,n=5)
    if(missing(mn)) mn <- min(major.ticks)
    if(missing(mx)) mx <- max(major.ticks)

    major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]

    labels <- sapply(major.ticks,function(i)
      as.expression(bquote(10^ .(i)))
    )
    axis(ax,at=major.ticks,labels=labels,
         las=ifelse(ax %in% c(2, 4), 2, 1), ...)

    n <- n+2
    minors <- log10(pretty(10^major.ticks[1:2],n))-major.ticks[1]
    minors <- minors[-c(1,n)]

    minor.ticks = c(outer(minors,major.ticks,`+`))
    minor.ticks <- minor.ticks[minor.ticks > mn & minor.ticks < mx]


    axis(ax,at=minor.ticks,tcl=par("tcl")*t.ratio,labels=FALSE)
  }

  null <- mapply(function(frag.len, frag.name){
    x <- 1:1010
    frag.len <- frag.len[match(x, names(frag.len))]
    frag.len[is.na(frag.len)] <- 0
    y <- frag.len / sum(frag.len)
    y <- as.numeric(y)
    names(y) <- x
    par(mar=c(5, 5, 4, 2) +.1)
    plot(x, y*10^3, main=paste(frag.name, "fragment sizes"),
         xlim=c(0, 1010), ylim=ylim,
         xlab="Fragment length (bp)",
         ylab=expression(Normalized ~ read ~ density ~ x ~ 10^-3),
         type="l")
    par(fig=c(.4, .95, .4, .95), new=TRUE)
    plot(x, log10(y), xlim=c(0, 1010), ylim=logYlim,
         xlab="Fragment length (bp)", ylab="Norm. read density",
         type="l", yaxt="n")
    minor.ticks.axis(2)
    par(opar)
  }, fragment.len, names(fragment.len))

  return(invisible(fragment.len))
}

# plotting
bamfile <- (args[1])
bamfile.labels <- gsub(".bam", "", basename(args[1]))

# Ab5L_bamfile <- ("/tgac/workarea/group-vh/Tarang/ATACseq/1.run1jan2017_twolanes/3a.post_alignment_filtering/Ab/PRO1563_S1_lib_CAGAATGC-GAACTGAG_L001_.nodup.filt.bam")
# Ab5L_bamfile.labels <- gsub(".bam", "", basename(Ab5L_bamfile))

temp <- fragSizeDist(bamfile, bamfile.labels)
ggsave(plot = temp, height = 7, width = 7) # this will output the default Rplots.pdf image which can be renamed in command line

# temp <- fragSizeDist(Ab5L_bamfile, Ab5L_bamfile.labels)
# ggsave(plot = temp, height = 7, width = 7)

# temp <- fragSizeDist(Ab5L_bamfile, Ab5L_bamfile.labels)
# ggsave(filename = "/tgac/workarea/group-vh/Tarang/ATACseq/1.run1jan2017_twolanes/3a.post_alignment_filtering/Ab/2Ba_Ab5L_FragmentSizeDist.pdf", plot = temp, height = 7, width = 7)

# # below is the traditional way of generating an image but this will not work (in this way) on the HPC as X11 is not available
# tiff('/tgac/workarea/group-vh/Tarang/ATACseq/1.run1jan2017_twolanes/3a.post_alignment_filtering/Ab/2Ba_Ab5L_FragmentSizeDist.tiff', units="in", width=10, height=10, res=300)
# temp
# dev.off()

#### An alternative method but has not been amended for the input format of fragment counts:
# ## ---------------------------
# ## Install relevant packages
# # install.packages(knitr)
# # install.packages(rmdformats)
# # install.packages(dplyr)
# # install.packages(DT)
# # install.packages(tidyr)
# # install.packages(ggplot2)
# # install.packages(magrittr)
# # install.packages(devtools)
# ## ---------------------------
# ## Load relevant packages
# library(magrittr)
# library(dplyr)
# library(ggplot2)
# ## ---------------------------
# ## ATAC-seq should represent a mix of fragment lengths corresponding to nucleosome free, mononucleosome and poly-nucleosome fractions.
# ## We can use insert lengths from across the genome to plot the distribution of all fragment lengths.
#
# fragLenPlot <- table(insertSizes) %>% data.frame %>% rename(InsertSize = insertSizes,
#     Count = Freq) %>% mutate(InsertSize = as.numeric(as.vector(InsertSize)),
#     Count = as.numeric(as.vector(Count))) %>% ggplot(aes(x = InsertSize, y = Count)) +
#     geom_line()
#
# # native scale plot
# A <- fragLenPlot + theme_bw()
#
# # log-scale plot
# B <- fragLenPlot + scale_y_continuous(trans = "log2") + theme_bw()
#
# ## ---------------------------
# ## Plotting insert sizes with open, mono- and di-nucleosome profiles
# ## We can now annotate our nucleosome free (< 100bp), mono-nucleosome (180bp-247bp) and di-nucleosome (315-437)
#
# # native scale
# C <- fragLenPlot + geom_vline(xintercept = c(180, 247), colour = "red") + geom_vline(xintercept = c(315,
#     437), colour = "darkblue") + geom_vline(xintercept = c(100), colour = "darkgreen") +
#     theme_bw()
#
# tiff('XX.tiff', units="in", width=10, height=10, res=300)
# C
# dev.off()
#
# # log-scale
# D <- fragLenPlot + scale_y_continuous(trans = "log2") + geom_vline(xintercept = c(180,
#     247), colour = "red") + geom_vline(xintercept = c(315, 437), colour = "darkblue") +
#     geom_vline(xintercept = c(100), colour = "darkgreen") + theme_bw()
#
# tiff('XX.tiff', units="in", width=10, height=10, res=300)
# D
# dev.off()
