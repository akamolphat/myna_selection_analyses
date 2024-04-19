# Library ---------------------------------------------------------------------
library(tidyverse)
# Read in .idxstats file to get cnotig lengths
#
# Only subset contigs which are part of the superscaffold
# 
# Remove chrZ and chrW
#
# Read in .idxstats file ------------------------------------------------------
i <- "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/DArTseq/data/processed/align_to_ref/idx_stats/10201_rep_a.idxstats"
dt <- read_delim(i, 
                  col_names = c("contig", "contig_length", 
                                "mapped_reads", "unmapped_reads") )

# Read in lmiss file ----------------------------------------------------------
# This is to see if the scaffold list made from .idxstats are present in the 
# WGS data as WGS data was called by Kat
dtx <- read_delim(file = "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/vcftools_sum_stats/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10/variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.lmiss")

# dt1 <- dt[((dt$contig_length > 100000 | grepl(x = dt$contig, pattern = "chr")) & !(dt$contig %in% c("Superscaffold_chrZ", "Superscaffold_chrW")) ), "contig"] 
# dt1$contig[!dt1$contig %in% unique(dtx$CHR)]
# write.table(dt1, file = "data/processed/BCFtools/vcftools_filtered/Contig_1kbp_chr.list",
#             quote = F, row.names = F, col.names = F)

dt2 <- dt[grepl(x = dt$contig, pattern = "Superscaffold") & !(dt$contig %in% c("Superscaffold_chrZ", "Superscaffold_chrW")), "contig"] 

ls <- dt2$contig[!dt2$contig %in% unique(dtx$CHR)]
if (length(ls) == 0) {
  write.table(dt2, file = "/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/Superscaffold.list",
              quote = F, row.names = F, col.names = F)
} else {
  print(paste("The following contigs are missing:", paste(ls, collapse = ", ")))
}




