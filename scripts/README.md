# scripts

This folder contain the scripts that were used to process the reads, call SNPs, and perform analyses on the data.

Folders are named sequentially with a two digit prefix and the scripts in each folder are executed in order. For example, codes in the `01_*` folder are executed before codes in the `02_*` folder. Within each folder, the codes are executed sequentially. Some folders contain more scripts than others and contain their own README.md file.

## Structure

```
├── 01_read_processing
│   ├── 01_xxx.sl
│   └── 02_xxx.sl
|
├── 02_baypass_WGS
│   ├── 01_xxx.sl
│   └── 02_xxx.sl
|
├── 03_EHHS
│   ├── 01_xxx.sl
│   └── 02_xxx.sl
|
├── 04_samplot
│   ├── 00_10X_to_ref
│   └── README.md
|
├── 05_baypass_SV
│   └── README.md
|
├── 06_genomic_features
│   └── README.md
|
├── 01_xxx.sl
└── README.md
```
## General overview

* `01_read_processing` - Processes the WGS reads, align them to the reference genome, mark duplicate reads, and call SNPs. 
* `02_baypass_WGS` - Converts SNP data into BayPass format, and performs BayPass core model with the `--contrastfile` argument
* `03_EHHS` - Performs IHS, and EHHS.
* `04_samplot` - Create samplot to identify and investigate SV region in all samples, including reference individual.
* `05_baypass_SV` - Converts manually genotyped SV to BayPass format and performs BayPass on it.
* `06_genomic_features` - Perform VEP, extract SV, TE, RepeatMasker, and sequences of interest for BLASTing on NCBI.

