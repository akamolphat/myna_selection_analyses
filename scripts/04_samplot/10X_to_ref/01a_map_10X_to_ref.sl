#!/bin/bash -e
#SBATCH --job-name=01a
#SBATCH --account=uoa02613
#SBATCH --time=01-00:00:00
#SBATCH --mem=50GB
#SBATCH --output=01_errors/01a_%j.out
#SBATCH --error=01_errors/01a_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --profile task
#SBATCH --partition=milan

# load modules
ml BWA/0.7.17-GCC-11.3.0
ml SAMtools/1.16.1-GCC-11.3.0

DIR=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/
Xfolder=align2ref/10X_to_ref/
newreffolder=align2ref/ref/

echo "Copy over ref fasta file"

GENOME=${DIR}VEP/ref/AcTris_vAus2.1.

cp ${GENOME}* ${DIR}${newreffolder}

cd ${DIR}${Xfolder}

mkdir -p 10X_mapped

cd 10X_processed/

SAMPLE=13099

# set paths
REF=${DIR}${newreffolder}AcTris_vAus2.1.fasta
TRIM_DATA_R1=${DIR}${Xfolder}10X_processed/Myna_10x_processed_R1_val_1.fq.gz
TRIM_DATA_R2=${DIR}${Xfolder}10X_processed/Myna_10x_processed_R2_val_2.fq.gz
OUT_DIR=${DIR}${Xfolder}10X_mapped/
bamout=${OUT_DIR}/${SAMPLE}.10X.sorted.bam
bamcov=${OUT_DIR}/${SAMPLE}.10X.chr8.coverage
cd ${OUT_DIR}

echo "Map the reads"
bwa mem -t ${SLURM_CPUS_PER_TASK} \
-R "@RG\tID:${SAMPLE}\tLB:${SAMPLE}_WGS\tPL:ILLUMINA\tSM:${SAMPLE}" \
-M ${REF} ${TRIM_DATA_R1} ${TRIM_DATA_R2} | \
samtools sort | samtools view -O BAM -o ${bamout}

samtools index ${bamout}

echo "Calculate depths from bam file"
samtools depth -r Superscaffold_chr8 ${bamout} > ${bamcov}

awk '$2 > 20650000 && $2 < 20700000' 13099.10X.chr8.coverage > 13099.10X.chr8.20650000_20700000.coverage

# Check output
# samtools flagstat ${OUT_DIR}/${SAMPLE}.10X.sorted.bam 
