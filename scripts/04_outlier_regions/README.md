# 04_outlier_regions

This folder define outlier regions and was used to help determine the "key" outlier regions.

* `01*.sl` - extract outliers based on the POD 99.999% threshold for XtX and C2, and -log<sub>10</sub> <i>p</i>-value > 6. This is performed on both MAC10 and MAC3 dataset.
* `02*.sl` - calls `02_cluster_outliers.R` to cluster outliers within 50 kb of one another. This is performed on both MAC10 and MAC3 dataset.
* `03*.sl` - makes a list of genes within 50 kb of the outliers, or the SNPs analysed, to be used in GO and protein class enrichment analysis.
* `04a_summarise_outlier_regions_mac10.R` - defines cluster of outliers defined in `02*` containing more than 5 outliers, defined as "outlier region". This was then used to see if there are outlier regions which were identified by XtX, and one of the C2-contrast statistics. Only CON001 and CON004 were used in the study. CON001 refers to C2-contrast statistics between native and invasive populations. CON004 refers to C2-contrast statistics between native and invasive populations founded by one-step introduction (i.e., Melbourne, Fiji, and South Africa) and is only mentioned in the discussion of the paper.
* `04b_extract_genelist_summarised_outlier_regions_mac10.sl` - extract list of genes within 50 kb of outlier regions defined in `04a*`.
* `04c_summarise_outlier_regions_genenames_mac10.R` - attach list of genes from `04a` to corresponding outlier regions.
