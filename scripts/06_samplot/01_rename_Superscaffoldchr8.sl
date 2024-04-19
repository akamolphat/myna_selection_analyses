#!/bin/bash -e
#SBATCH --job-name=01
#SBATCH --output=01_errors/01_%j.out
#SBATCH --error=01_errors/01_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4
#SBATCH --array=1-82

module purge

ml SAMtools/1.16.1-GCC-11.3.0

folder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/mapped_reads/
cd ${folder}
infile=$(find J*[0-9].sorted.dup.bam -printf "%f\n" | sed "${SLURM_ARRAY_TASK_ID}q;d")
outfile=${infile//.sorted.dup.bam/.sorted.dup.reheader.bam}

samtools view -H ${infile} | sed -e 's/SN:Superscaffold_chr8/SN:chr8/' | samtools reheader - ${infile} > ${outfile}
samtools index ${outfile}

