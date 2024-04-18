# 04_samplot

The scripts in this folder are arranged as follows:

```
├── 00_10X_to_ref
│   ├── 01a_map_10X_to_ref.sl
│   └── 01b_markduplicates.sl
├── 00_markdup_rename_refindiv_10X.sl
├── 01_rename_Superscafoldchr8.sl
├── 02a_samplot_iter1.sl
├── 02b_samplot_example_plot.sl
└── README.md
```

The scripts in this folder were executed in the following order:

* Scripts within the `00_10X_to_ref/` folder - the scripts in this folder aligns the Chromium 10X reads to the reference genome and mark the duplicates.
* `00_markdup_rename_refindiv_10X.sl` - This script renames the contig name in the sorted.dup.bam file of the reference individual. This is somehow necessary in the version of samplot that was used. When this is not performed, the samplot command does not work as it appears to assume that the contig names start with "chr".
* `01_rename_Superscaffoldchr8.sl` - This script renames the contig name in the sorted.dup.bam file of all the WGS samples, for the same reason as for `00_*.sl `.
* `02a_samplot_iter1.sl` - This script plots samplot for each sample within each population, creating one file per individual, and one file for the reference individual.
* `02b_samplot_example_plot.sl` - This script plots samplot for the samples used to present in the Supplementary.
