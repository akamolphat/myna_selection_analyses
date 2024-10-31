# 01_trim_reads

Barcodes and adapters were trimmed using TRIM_GALORE version 0.6.7 (Krueger et al., 2021). For each paired-end sample, the following command was run:

```
file1=fileR1.fasta.gz # file containing R1
file2=fileR2.fasta.gz # file containing R2
trim_galore --2colour 20 --paired --fastqc ${file1} ${file2} 
```

The output files have “val_1” and “val_2” appended to the file name by default by trimgalore (e.g. fileR1_val_1.fq.gz).

