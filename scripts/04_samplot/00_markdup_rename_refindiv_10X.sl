#!/bin/bash -e
#SBATCH --job-name=00
#SBATCH --output=00_errors/00_%j.out
#SBATCH --error=00_errors/00_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=1-00:00:00
#SBATCH --mem=100G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=2

module purge

ml SAMtools/1.16.1-GCC-11.3.0
ml picard/2.26.10-Java-11.0.4

folder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/align2ref/10X_to_ref/10X_mapped/
cd ${folder}

infile=13099.10X.sorted.bam
outfiledup=${infile//.sorted.bam/.sorted.dup.bam}
outfilereheader=${outfiledup//.sorted.dup.bam/.sorted.dup.reheader.bam}

picard MarkDuplicates INPUT=${folder}/${infile} OUTPUT=${folder}/${outfiledup} METRICS_FILE=${folder}/13099.10X.metrics.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000;
samtools index -@ 2 ${folder}/${outfiledup}

samtools view -H ${outfiledup} | sed -e 's/SN:Superscaffold_chr8/SN:chr8/' | samtools reheader - ${infile} > ${outfilereheader}
samtools index ${outfilereheader}

