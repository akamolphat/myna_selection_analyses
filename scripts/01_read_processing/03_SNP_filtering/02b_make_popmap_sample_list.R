# Libraries ---------------------------------------------------------
library(tidyverse)
# This script makes the population map file 
#
# On first glance, there are no samples with outlier heterozygosity or 
# missingness but closer look found one sample with outlier 
# heterozygosity. This was sample M02
# 
# However, there are one sample without coordinates
# These were labelled as "CAI" in the WGS dataset, but the
# ID number suggests that these were "Australian (unknown)"
# Kyle did not provide location for ID number 11626,
# but 11625 was from Cairns. There are a few of these samples
# such as sample 11737 which appeared together with a set of 
# samples. Perhaps, Annabel had some additional information on 
# where these samples are from (e.g. based on original sample ID?)
# 
# Based on heterozygosity
het <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/vcftools_sum_stats/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.het"
ind_het <- read_delim(het, delim = "\t",
                      col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)
# Make a list of individuals to remove, if any
ls_het_outlier <- ind_het$ind[ind_het$f < -0.2] # 10519, 10550, 12718, 12719, M0208, M0271 from Arthur's Creek, Arthur's Creek, Sydney, Sydney, Fiji, Odisha
#
#
# Read in individual missingness file 
imiss <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/vcftools_sum_stats/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.imiss"
ind_miss  <- read_delim(imiss, delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)
# Remove sample with outlier Fis (J800)
ind_miss <- ind_miss[!ind_miss$ind %in% ls_het_outlier,]
# Read in metadata file 
dt_meta <- read_csv("/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/metadata/myna_WGS_meta_pop_env.csv")
dt_meta <- merge(dt_meta, ind_miss[,"ind"], by.x = "Name", by.y = "ind")


# Extract in population map file, based on dt_meta
dt_popmap <- dt_meta[,c("Name", "pop3_abb")] 
colnames(dt_popmap) <- c("ind", "pop3_abb")

popmap <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/pop_map_WGS_ALL.txt"
write_delim(dt_popmap, popmap, delim = "\t", col_names = F)

dt_ind2keep <- dt_popmap[!dt_popmap$pop3_abb %in% c("Australia_unknown"),]
dt_ind2keep <- as_tibble(dt_ind2keep)
# Write a list of samples to keep -----------------------------------
write_delim(dt_ind2keep[,c("ind")], file = "data/processed/BCFtools/vcftools_filtered/filtsamp_2keep.txt", delim = "\t", col_names = F)

dir.create("data/processed/BCFtools/vcftools_filtered/pop_popmaps", showWarnings = F, recursive = T)

for (i in unique(dt_ind2keep$pop3_abb)){
  print(i)
  dt_sub <- dt_ind2keep[dt_ind2keep$pop3_abb == i,]
  outputfile <- paste("data/processed/BCFtools/vcftools_filtered/pop_popmaps/", i, "_samples.txt", sep = "")
  write_delim(dt_sub[,c("ind")], outputfile, delim = "\t", col_names = F)
}
