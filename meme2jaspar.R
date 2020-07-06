#!/usr/bin/env Rscript

##########################################################################################################
# Author: T K Mehta
# Date: July 2020
# Decsription: converts meme format PWM to JASPAR 2016 PFM format
##########################################################################################################

# setwd("~/Documents/TGAC/Projects/Cichlid_GRNs/Arboretum_GT_v3/1.TFBSs_v2/FINAL_cichlidPWM_motifs/MouseDerived/CS/ab")

# if (!requireNamespace("BiocManager", quietly=TRUE))
#   install.packages("BiocManager")

# BiocManager::install("optparse")
library(optparse)
# BiocManager::install("universalmotif")
library(universalmotif)

option_list = list(
  make_option(c("-i", "--input"), action="store", default=NA, type='character',
              help="input *.meme file"),
  make_option(c("-o", "--output"), action="store", default=NA, type='character',
              help="output filename as *.pwm")
)
opt = parse_args(OptionParser(option_list=option_list))

memepwm = read_meme(opt$i)
write_jaspar(memepwm, opt$o)
# test1 <- read_meme("RG-cich_GTRDdata_mouse_sites_TF_ig_ab.gene.s0.84.ig_1_CS_ab.meme")
# write_jaspar(test1, 'test.pwm')

