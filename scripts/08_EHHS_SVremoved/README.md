# 08_EHHS_SVremoved

This folder contain scripts which performed the following:

* `01_Shapeit5_chr8_repeat482137removed.sl` - performs statistical phasing using SHAPEIT5 on chromosome 8 with SNPs from the SV region remove (AKA repeat 482137 in RepeatMasker).
* `02a_get_IHS_chr8_repeat482137removed.sl` - calculates EHH, iHH, and iHS based on the outputs from `01*`
* `02b_get_IHH_chr8_repeat482137removed_pos_shifted.R` - shifts the SNPs after the SV region by the length of the SV, then calculates EHH and iHH.
* `03_combine_IHS_WGS_chr8_SVremoved.R` - calculates iHS (normalised with the rest of the genome) using iHH outputted from both `02a*` and `02b*`. This result in two different files. Only the results from `02b*` were reported in the manuscript.
