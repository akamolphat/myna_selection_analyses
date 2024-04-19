#!/usr/bin/env Rscript
#
# This script calculate outlier thresholds based on particular column
#
# Define input arguments --------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
infile = args[1]
bp_dist = as.integer(args[2])
bp_region_dist = as.integer(args[3])
outclusterfile = args[4]
outsummaryfile = args[5]
outbedfile= args[6]

print(c(infile, bp_dist, bp_region_dist, outclusterfile, outsummaryfile, outbedfile))
# infile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/combined/myna_baypass_mac10_combined_xtx_outliers_0.99999.txt"
# bp_dist <- 100000  # distance between outliers to be grouped together
# bp_region_dist <- 250000  # distance around outlier regions to extract genes from
# outclusterfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/combined/myna_baypass_mac10_combined_xtx_outliers_0.99999_cluster.txt"
# outsummaryfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/combined/myna_baypass_mac10_combined_xtx_outliers_0.99999_cluster_summary.txt"
# outbedfile <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/combined/myna_baypass_mac10_combined_xtx_outliers_0.99999_outlier_region.bed"

# Library -----------------------------------------------------------------
library(data.table)
library(tidyverse)

# Define functions --------------------------------------------------------
clustPOS <- function(dt, thresh = 100000, clustname = "cluster1"){
  #
  # This function cluster points based on distances between outlier SNPs
  #
  ctr <- 1
  # thresh <- 100000
  # clustname <- "cluster1"
  for (i in 1:(nrow(dt) - 1)){
    if (dt$chr[i] != dt$chr[i + 1]){
      dt[[clustname]][i] <- ctr
      ctr <- ctr + 1
      dt[[clustname]][i + 1] <- ctr
    } else {
      if (dt$POS[i + 1] - dt$POS[i] > thresh){
        dt[[clustname]][i] <- ctr
        ctr <- ctr + 1
        dt[[clustname]][i + 1] <- ctr
      } else {
        dt[[clustname]][i] <- ctr
      }
    }
  }
  return(dt)
}




# Read file ---------------------------------------------------------------
dt <- fread(infile, header = T)

# Cluster points ----------------------------------------------------------
dt_outlier_cluster <- clustPOS(dt, thresh = bp_dist, clustname = "cluster1")

dt_outlier_region <- dt_outlier_cluster |>
  dplyr::group_by(chr, cluster1) |>
  dplyr::summarise(start_pos = min(POS),
            end_pos = max(POS),
            n = n())

dt_outlier_bed <- dt_outlier_region |>
  dplyr::mutate(start_pos = start_pos - bp_region_dist,
                end_pos = end_pos + bp_region_dist) 
dt_outlier_bed$start_pos[dt_outlier_bed$start_pos < 0] <- 0

fwrite(x = dt_outlier_cluster, file = outclusterfile, sep = "\t")
fwrite(x = dt_outlier_region, file = outsummaryfile, sep = "\t")
fwrite(x = dt_outlier_bed[,c("chr", "start_pos", "end_pos")], file = outbedfile, col.names = F, sep = "\t") 
# write.table(x = dt_outlier_bed, file = outbedfile, col.names = F, sep = "\t", row.names = F, quote = F)  # fwrite gives some issues, likely to do with the missing column names
