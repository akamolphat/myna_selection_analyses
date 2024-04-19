# 03_EHHS_WGS

This folder performs the following:

* `01_Shapeit5.sl` - performs statistical phasing based on SHAPEIT5 (each chromosome separately).
* `02_get_IHS.sl` - calculates EHH, iHH, and iHS for each chromosome separately and outputs to a file.
* `03_combine_IHS_WGS.R` - combine iHH outputs from `02_*` and calculates IHS based on whole-genome values instead of each chromosome. This effects the normalisation of iHS.
* `04_pval_hist.R` - plots p-value histogram to see if FDR correction is appropriate.
