#!/usr/bin/env Rscript
# Define input arguments --------------------------------------------#####
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2){
  stop("Number of arguments are not correct. Seven arguments required: \n
       metafile \n
       cfilepath")
}

metafile = args[1]   # "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/metadata/myna_WGS_meta_pop_env.csv"
cfilepath = args[2]   # "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/contrastfile.txt"

# metafile = "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/metadata/myna_WGS_meta_pop_env.csv"
# cfilepath = "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/contrastfile.txt"


# This script generate contrast files
#
# Library ---------------------------------------------------------------------
library(tidyverse)

# Read input files ------------------------------------------------------------
dt_ind <- read_csv(metafile)

# Create contrast column ------------------------------------------------------
dt_popenvsub$pop_inv_nat1 <- -1  # Define all that is not native as -1
## Define populations ---------------------------------------------------------
nat_pop <- c("Madhya_Pradesh", "Maharashtra", "Tamil_Nadu", "Maharashtra_subpopulation_A")
indother <- c("Madhya_Pradesh", "Maharashtra", "Tamil_Nadu")
subpopAinv <- c("Maharashtra_subpopulation_A", "Melbourne", "Cairns", "Fiji", "Leigh", "Great_Barrier_Island", "Napier")
melbinv <- c("Melbourne", "Cairns", "Leigh", "Great_Barrier_Island", "Napier")


## Native vs Invasive ---------------------------------------------------------
dt_ind$pop_inv_nat1[dt_ind$pop3_abb %in% natpop] <- 1
## IND: Other vs Maharashtra subpop. A branch of invasion ---------------------
dt_ind$pop_inv_nat2 <- 0
dt_ind$pop_inv_nat2[dt_ind$pop3_abb %in% indother] <- 1
dt_ind$pop_inv_nat2[dt_ind$pop3_abb %in% subpopAinv] <- -1
## Native vs Melbourne branch of invasion -------------------------------------
dt_ind$pop_inv_nat3 <- 0
dt_ind$pop_inv_nat3[dt_ind$pop3_abb %in% nat_pop] <- 1
dt_ind$pop_inv_nat3[dt_ind$pop3_abb %in% melbinv] <- -1
## Native + Melb vs NZ --------------------------------------------------------
# NZ invasion specific outliers 
dt_ind$pop_inv_nat4 <- 0
dt_ind$pop_inv_nat4[dt_ind$pop3_abb %in% c(nat_pop, "Melbourne")] <- 1
dt_ind$pop_inv_nat4[dt_ind$pop3_abb %in% c("Leigh", "Great_Barrier_Island", "Napier")] <- -1
## Native vs Cairns + Fiji ----------------------------------------------------
# Tropical invasion specific outliers
dt_ind$pop_inv_nat5 <- 0
dt_ind$pop_inv_nat5[dt_ind$pop3_abb %in% nat_pop] <- 1
dt_ind$pop_inv_nat5[dt_ind$pop3_abb %in% c("Cairns", "Fiji")] <- -1
## Native vs Melb + NZ --------------------------------------------------------
# Temperate invasive specific outliers
dt_ind$pop_inv_nat6 <- 0
dt_ind$pop_inv_nat6[dt_ind$pop3_abb %in% nat_pop] <- 1
dt_ind$pop_inv_nat6[dt_ind$pop3_abb %in% c("Melbourne", "Napier", "Leigh", "Great_Barrier_Island")] <- -1

## Native vs Fiji + SA + Melb
dt_ind$pop_inv_nat7 <- 0
dt_ind$pop_inv_nat7[dt_ind$pop3_abb %in% nat_pop] <- 1
dt_ind$pop_inv_nat7[dt_ind$pop3_abb %in% c("Melbourne", "Fiji", "South_Africa")] <- -1

popcontrast <- c("pop_inv_nat1", "pop_inv_nat2", "pop_inv_nat3", "pop_inv_nat4", "pop_inv_nat5", "pop_inv_nat6", "pop_inv_nat7")

t_dt_popcontrast <- data.frame(t(dt_popenvsub[, popcontrast]))

write.table(t_dt_popcontrast, cfilepath, sep = "\t", row.names = F, col.names = F)


