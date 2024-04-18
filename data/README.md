# data

This folder houses the data used in the analysis. The folders are currently empty, but the raw data were downloaded from the following:
* [Reference genome and annotations](https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_037013685.1/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GFF&include_annotation_type=RNA_FASTA&include_annotation_type=CDS_FASTA&include_annotation_type=PROT_FASTA&include_annotation_type=SEQUENCE_REPORT&hydrated=FULLY_HYDRATED) = NCBI accession GCA_037013685.1
* [Illumina WGS reads previously submitted with reference genome paper](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1054049) = NCBI Bioproject accession PRJNA1054049
* Illumina WGS reads
* Chromium 10X reads
* VGP reference genome
* TableS1.2.csv = metadata file

The folder is arranged as follows:

```
├── raw_data
│   └── README.md
│   
├── processed
│   └── README.md
|
├── TableS1.2.csv
|
└── README.md
```

TableS1.2.csv is part of supple
