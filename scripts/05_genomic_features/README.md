# 05_genomic_features

This folder contain the scripts that do the following:

* `01_VEP_extract_missense_chr8.sl` - performs Variant Effect Predictor (VEP) on all SNPs in the MAC10 dataset, then extracts the missense mutations within the key outlier region on chromosome 8.
* `02_extract_SV_chr8.sl` - extracts previously identified SVs in the key outlier region on chromosome 8 from a VCF containing SVs.
* `03_extract_TE_chr8.sl` - extracts previously identified TEs in the key outlier region on chromosome 8 from a VCF containing TEs.
* `04_extract_RepeatMasker_chr8.sl` - extracts repeat regions within the proximity to the key outlier region, previously identified in Stuart et al. (2024).
* `05_get_gene_sequence_from_ref.sl` - extract nucleotide sequences of genes from within the key outlier region.
* `06_get_sequence_SV_to_20700000_chr8.sl` - extract nucleotide sequence of the peak outlier region after the SV (after the SV to position 20.7 Mb).
* `07_extract_maf.sl` - extract allele frequency of each population for all SNPs on chromosome 8. This is only to help with plotting of allele frequency heatmap.

