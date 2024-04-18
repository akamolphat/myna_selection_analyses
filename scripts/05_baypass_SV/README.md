# 05_baypass_SV

This folder performs baypass on the genotyped SV

* `01_SV_GT_to_baypass.R` - This script converts the manually genotyped SV (from TableS1.2.csv) into the BayPass format. Essentially, the BayPass input file contains the allele count of SV insertions vs SV deletions in each population.
* `02_Baypass_core_C2_SV.sl` - This script performs BayPass, calculating XtX and C2.
