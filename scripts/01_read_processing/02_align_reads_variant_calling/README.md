# 02_align_reads_variant_calling

## Read alignments
Trimmed reads were aligned to the common myna reference genome, available on NCBI (accession GCA_037013685.1) (Stuart et al., 2024) using BWA version 0.7.17 mem (Li, 2013), before being processed by SAMTOOLS into sorted BAM files for each paired-end samples as follows:
```
GENOME=path/to/ref_genome.fasta # Path to reference genome
TRIM_DATA_R1=path/to/trimmed_readsR1.fq.gz # e.g. fileR1_val_1.fq.gz
TRIM_DATA_R2=path/to/trimmed_readsR2.fq.gz # e.g. fileR2_val_2.fq.gz
OUT_DIR=path/to/outputfolder/   # Path to output folder
SAMPLE=sampleID # e.g. J772
# Align. NOTE that -t 16 command is just the no. of threads

bwa mem -t 16 -R “@RG\tID:${SAMPLE}\tLB:${SAMPLE}_WGS\tPL:ILLUMINA\tSM:${SAMPLE}” -M ${GENOME} ${TRIM_DATA_R1} ${TRIM_DATA_R2} | samtools sort | samtools view -O BAM -o ${OUT_DIR}/${SAMPLE}.sorted.bam 
```
The following command can also be run to check some of the outputs
```
samtools flagstat ${OUT_DIR}/${SAMPLE}.sorted.bam 
```

## Variant calling
Duplicate reads were marked using PICARD version 2.26.10 (Picard Toolkit, 2019) MarkDuplicates using the following command:

```
SAMPLE=sampleID # e.g. J772
OUT_DIR=path/to/outputfolder/   # Path to output folder

# Mark Duplicates

picard MarkDuplicates INPUT=${OUT_DIR}/${SAMPLE}.sorted.bam OUTPUT=${OUT_DIR}/${SAMPLE}.sorted.dup.bam METRICS_FILE=${OUT_DIR}/${SAMPLE}.metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000;

# Index sorted.dup.bam files

samtools index -@ 2 ${OUT_DIR}/${SAMPLE}.sorted.dup.bam
```

Variants were jointly called across samples using BCFtools v1.13 (Danecek et al., 2021) mpileup (-C 50 -q 20 -Q 25), call and view functions. The following commands were executed:

```
GENOME=path/to/ref_genome.fasta # Path to reference genome
OUT_DIR=path/to/outputfolder/   # Path to folder to store VCF
# Path to textfile storing list of path to bamfiles (one filepath per line), e.g. ${OUT_DIR}/sample_bamfiles_list.txt  
BAM_LIST=path/to/bamfilelist.txt 
OUT_VCF=${OUT_DIR}myna_82inds.vcf.gz #Path to output VCF

# Calling variants

bcftools mpileup -C 50 -q 20 -Q 25 -a 'DP,AD,ADF,ADR,SP' -Ou -f ${GENOME} -b ${BAM_LIST} | \
bcftools call -c | \
bcftools view --exclude-types indels | \
bcftools sort --temp-dir ${OUT_DIR}/temp -Oz -o ${OUT_VCF}
```
