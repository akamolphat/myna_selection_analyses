# 02_baypass_WGS

The scripts in this folder are used to perform BayPass core model with the `-contrastfile` flag.

*`01*_vcf2plink_*.sl` - convert VCF file to plink format
*`02*_makebaypass_*.sl` - convert plink format to BayPass format
*`03*_subset_baypass_*.sl` - create sub-datasets, 1 in every 100 SNPs for MAC10, and 1 in every 170 SNPs for MAC3 dataset.
*`04*_subset_SNP_meta_*.sl` - create sub-metadata in the same way as `03*`, to make sure that SNP metadata can be re-assembled.
*`05_create_contrastfile.*` - create contrast file for the `-contrastfile` flag.
*`06a*_Baypass_core_contrast_*.sl` - perform BayPass core model with the contrastfile on each sub-dataset.
*`06b*_simulate_POD_core_contrast_*.sl` - simulate POD data based on omega of the first sub-dataset calculated in `06a*`, and performs the same BayPass core model with the contrastfile on the POD as in `06a*`.
*`07a*_combine_XtX_*.sl` - reassemble the XtX outputs from all the sub-datasets.
*`07b*_combine_C2_*.sl` - assemble the C2-contrast outputs from all the sub-datasets. Results from each C2-contrast comparisons were separated into different files.

