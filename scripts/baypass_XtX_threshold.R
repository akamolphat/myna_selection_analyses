#!/usr/bin/env Rscript
# Define input arguments --------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
in_pi_xtx = args[1]
threshold = as.numeric(args[2])
# wdir = args[3]
# print(c(in_pi_xtx, threshold, wdir))

# setwd(wdir)
# Libraries ---------------------------------------------------------
library("ape")
library("corrplot")

pod.xtx <- read.table(in_pi_xtx, header = T)

pod.thresh <- quantile(pod.xtx$M_XtX ,probs = threshold)

cat(unname(pod.thresh))