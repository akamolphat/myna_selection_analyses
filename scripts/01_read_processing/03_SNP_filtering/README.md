# 03_SNP_filtering

The scripts in this folder are arranged as follows:

## 01* scripts
The scripts with this prefix are used for basic filtering and determine read depth threshold for subsequent filters. 

* `01a_filter_bial_SNPs_QUAL30.sl` - filters for biallelic SNPs with QUAL > 30 (`--minQ 30`)
* `01ai_plot_bial_SNPs_QUAL30.sl` - **not necessary** - plots some basic stats after filtering. Require much higher RAM for shorter amount of time.
* `01aii_filter_minDP5.sl` - subsequently recodes genotype with DP < 5 to NA (`--minDP --recode --recode-INFO-all`)
* `01aiii_plot_minDP5.sl` - **not necessary** - plots some basic stats after filtering. Require much higher RAM for shorter amount of time.
* `01aiv_filter_nomissingness.sl` - subsequently filters for SNPs with no missingness. This is to see get a DP distribution of what "good" SNPs are supposed to look like.
* `01av_plot_nomissingness.sl` - plots read depth distribution and calculates 95% quantiles based on a subsample of 50000 SNPs to speed up the process.

## 02* scripts
The scripts with this prefix are used for the rest of the filtering, except for the MAC filter where two datasets were made.

* `02a_filter_DP5-35_nosingletonsdoubletons.sl` - recodes outputs from `01a_filter_bial_SNPs_QUAL30.sl` so that genotypes with DP < 5 and > 35 are recoded to NA. Singletons and doubletons were also removed.
* `02ai_plot_DP5-35_nosingletonsdoubletons.sl` - **not necessary** - plots some basic stats after filtering. Require much higher RAM for shorter amount of time.
* `02b_make_popmap_sample_list.R` - makes population map file (sample ID in first column, and population in second column)
* `02c_filter_Cairns_uncertain.sl` - removes sample J772 which were supposedly from Cairns, but not present in previous DArTseq analysis. It also does not cluster with other Cairns samples.
* `02d_split_VCF_by_pop.sl` - splits the VCF by population as defined in the population map file.
* `02di_plot_VCF_by_pop.sl` - **not necessary** - plots some basic stats for each population.
* `02ei_filter_missingness5n_per_pop.sl` - filters for SNPs genotyped in at least 5 individuals per population. This is done on the VCF file containing data from each population. 
* `02eii_plot_missingness5n_per_pop.sl` - **not necessary** - plots some basic stats for each population.
* `02fi_filter_fullVCFSNPlist5n.sl` - creates a list of SNPs that are present in all the filtered VCF of each population (after `02ei_filter_missingness5n_per_pop.sl`). The list of SNPs are then used to filter the VCF file (after `02c_filter_Cairns_uncertain.sl`) to retain only SNPs present in at least 5 individuals per population.
* `02fii_plot_fullVCFSNPlist5n.sl` - **not necessary** - plots some basic stats for each population.
* `02g_make_contig_list.*` - makes a list of contigs to keep (contigs at autosomal level)
* `02h_filter_contig_list.sl` - filter for SNPs only on autosomal contigs.

## 03* scripts
MAC filters to create two datasets used for downstream analyses.

* `03ai_filter_mac10.sl` - filter for MAC >= 10, creating the MAC10 dataset
* `03aii_filter_mac3.sl` - filter for MAC >= 3, creating the MAC3 dataset
