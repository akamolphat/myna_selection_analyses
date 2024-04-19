#!/usr/bin/env Rscript
# Define input arguments --------------------------------------------#####
hapfile <- "data/processed/BCFtools/vcftools_filtered/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.SVsnpsremoved.phased_shapeit5.snps.vcf.gz"
scanout <- "data/processed/EHHS/chr8_SVsnpsremoved_pos_shifted_scanhh.txt"


# hapfile <- "data/processed/BCFtools/vcftools_filtered/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.phased_shapeit5.snps.vcf.gz"
# IHSout <- "data/processed/EHHS/chr8_IHS.txt"

# Load libraries ----------------------------------------------------
library(rehh)
library(data.table)
library(tidyverse)

hh <- data2haplohh(hap_file = hapfile,
                   polarize_vcf = FALSE,
                   vcf_reader = "data.table")

hh@positions[hh@positions > 20678727] <- hh@positions[hh@positions > 20678727] - 8798

res.scan <- scan_hh(hh)
fwrite(res.scan, file = scanout, sep = "\t", quote = F, row.names = T)


