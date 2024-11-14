#!/usr/bin/env Rscript
# Define input arguments --------------------------------------------#####
args <- commandArgs(trailingOnly = TRUE)
inputvcf = args[1]
nloci = as.numeric(args[2])
outputpdf = args[3]

# This script plots a histogram of the DP values in a VCF file
#
# The purpose of this script is to see the distribution DP in the dataset
# given the GQ filters. The main purpose is to see what "good" loci 
# is supposed to look like when GQ or QUAL is high. 
# This script aims to identify cutoff points for DP filters
#
# This script compares the results for what GQ, and especially DP should
# look like in a dataset with "good" loci, that is with minimum GQ scores
# of 20 and was called in all individuals.
# For the BCFtools dataset, I also filter for only QUAL > 30.

# Load libraries ----------------------------------------------------#####
library(vcfR) 
library(ggplot2) 
# Read input VCF data -----------------------------------------------#####
inputfile <- inputvcf
vcf <- read.vcfR(inputfile, verbose = FALSE) 
# Parse DP, GQ and GT from the GT region of the VCF 
dp <- extract.gt(vcf, element="DP", as.numeric = TRUE) 

dp_vec <- c(dp[!is.na(dp)])
dp_vec_quant <- quantile(dp_vec, c(0,0.01, 0.025, 0.05, 0.25, 0.5, 0.75, 0.95, 0.975, 0.99))
idx <- sample(seq(1,length(dp_vec)), nloci)
dt_dp <- data.frame(DP = dp_vec[idx])
p1 <- ggplot(dt_dp, aes(DP)) +
  geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) +
  xlim(c(0,250)) +
  geom_vline(aes(xintercept=dp_vec_quant['50%']),
             linetype=1, size=1) + 
  geom_vline(aes(xintercept=dp_vec_quant['2.5%']),
             linetype=2, size=1) +
  geom_vline(aes(xintercept=dp_vec_quant['97.5%']),
             linetype=4, size=1) +
  annotate("text", x = round(dp_vec_quant['2.5%']), y = 0.01, 
           label = paste("2.5% quantile =", round(dp_vec_quant['2.5%'])), angle = 90, vjust = -0.5) +
  annotate("text", x = round(dp_vec_quant['50%']), y = 0.01, 
           label = paste("Median =", round(dp_vec_quant['50%'])), angle = 90, vjust = -0.5) +
  annotate("text", x = round(dp_vec_quant['97.5%']), y = 0.01, 
           label = paste("97.5% quantile =", round(dp_vec_quant['97.5%'])), angle = 90, vjust = -0.5) +
  theme_bw() 

pdf(outputpdf, width = 8, height = 5)
p1
dev.off()
