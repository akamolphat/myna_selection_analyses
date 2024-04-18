# scripts

This folder contain the scripts that were used to process the reads, call SNPs, and perform analyses on the data.

The scripts in each folder are run sequentially. 

This folder is organised as follows:

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
│   ├── 10X_to_ref
│   │   ├── 01a_map_10X_to_ref.sl
│   │   └── 01b_markduplicates.sl
│   ├── 00_markdup_rename_refindiv_10X.sl
│   ├── 01_rename_Superscafoldchr8.sl
│   ├── 02a_samplot_iter1.sl
│   ├── 02b_samplot_example_plot.sl
│   └── README.md
|
├── 05_baypass_SV
│   ├── 01_xxx.sl
│   └── 02_xxx.sl
|
├── 06_VEP_TE_SV_RepeatMasker
│   ├── 01_xxx.sl
│   └── 02_xxx.sl
|
├── 01_xxx.sl
└── README.md
```
