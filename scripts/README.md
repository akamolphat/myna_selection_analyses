# scripts

This folder contain the scripts that were used to process the reads, call SNPs, and perform analyses on the data.

Folders are named sequentially with a two digit prefix and the scripts in each folder are executed in order. For example, codes in the `01_*` folder are executed before codes in the `02_*` folder. Within each folder, the codes are executed sequentially. Some folders contain more scripts than others and contain their own README.md file.

## Structure

```
├── 01_read_processing
│   ├── 01_trim_reads
│   ├── 02_align_reads_variant_calling
│   └── 02_SNP_filtering
|
├── 02_baypass_WGS
│   ├── 01_xxx.sl
│   ├── ......
│   └── README.md
|
├── 03_EHHS_WGS
│   ├── 01_xxx.sl
│   └── README.md
|
├── 04_outlier_regions
│   ├── 01_xxx.sl
│   ├── ......
│   └── README.md
|
├── 05_genomic_features
│   ├── 01_xxx.sl
│   ├── ......
│   └── README.md
|
├── 06_samplot
│   ├── 00_10X_to_ref
│   ├── 01_xxx.sl
│   ├── ......
│   └── README.md
|
├── 07_baypass_SV
│   ├── 01_xxx.sl
│   ├── ......
│   └── README.md
|
├── 08_phylogenetic_tree_and_gene_tree
│   ├── 01_xxx.sl
│   ├── ......
│   └── 03_genetrees.md
|
├── 09_make_pots
│   ├── 01_xxx.sl
│   ├── ......
│   └── README.md
|
├── 01_xxx.sl
└── README.md
```
## General overview

* `01_read_processing` - Processes the WGS reads, align them to the reference genome, mark duplicate reads, and call SNPs. 
* `02_baypass_WGS` - Converts SNP data into BayPass format, and performs BayPass core model with the `--contrastfile` argument
* `03_EHHS` - Performs EHHS.
* `04_outlier_regions` - Define outlier regions.
* `05_genomic_features` - Perform VEP, extract SV, TE, RepeatMasker, and sequences of interest for BLASTing on NCBI.
* `06_samplot` - Create samplot to identify and investigate SV region in all samples, including reference individual.
* `07_baypass_SV` - Converts manually genotyped SV to BayPass format and performs BayPass on it.
* `08_phylogenetic_tree_and_gene_tree` - Creates phylogenetic tree based on thinned WGS data and SNPs from key outlier region, and create gene trees based on different copies of AMY2A in multiple reference genomes.
* `09_make_plots` - Make figures used in the manuscript and supplementary materials

