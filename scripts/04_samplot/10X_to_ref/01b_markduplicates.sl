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
ml SAMtools/1.16.1-GCC-11.3.0
ml picard/2.26.10-Java-11.0.4

DIR=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/
Xfolder=align2ref/10X_to_ref/
newreffolder=align2ref/ref/

echo "Copy over ref fasta file"

SAMPLE=13099

# set paths
REF=${DIR}${newreffolder}AcTris_vAus2.1.fasta
OUT_DIR=${DIR}${Xfolder}10X_mapped/
bamout=${OUT_DIR}/${SAMPLE}.10X.sorted.bam
bamdupout=${OUT_DIR}/${SAMPLE}.10X.sorted.dup.bam

cd ${OUT_DIR}


picard MarkDuplicates INPUT=${bamout} OUTPUT=${bamdupout}
samtools index ${bamdupout}
