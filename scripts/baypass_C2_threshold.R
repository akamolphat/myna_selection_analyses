#!/usr/bin/env Rscript
#
# This script calculate outlier thresholds based on particular column
#
# Define input arguments --------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
infile = args[1]
colname = args[2]
contrast_id = as.integer(args[3])
threshold = as.numeric(args[4])
# Library -----------------------------------------------------------
library(data.table)

# Read table --------------------------------------------------------
dt <- fread(infile, header = T) 
dt <- dt[dt$CONTRAST == contrast_id,]

# Get threshold -----------------------------------------------------
pod.thresh <- quantile(dt[[colname]], probs = threshold)

# Output to command line --------------------------------------------
cat(unname(pod.thresh))