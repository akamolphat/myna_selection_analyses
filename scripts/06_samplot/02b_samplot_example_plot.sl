#!/bin/bash -e
#SBATCH --job-name=02c
#SBATCH --output=02_errors/02c_%j.out
#SBATCH --error=02_errors/02c_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4

module purge
module load Miniconda3/22.11.1-1

source activate /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/programs/miniconda/envs/mamba 

chr=chr8
wdir=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/Potential_repeats/
mapdir1=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/mapped_reads/
cd ${wdir}
mapdir2=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/align2ref/10X_to_ref/10X_mapped/
outputfolder=repeat_ID_482137/example_genotyping
samplot plot -o ${outputfolder} -c ${chr} -s 20678727 -e 20687524 -t DEL -d 1000 -n Reference "Homozygous insertion" "Heterozygous" "Homozygous deletion" -b ${mapdir2}13099.10X.sorted.dup.reheader.bam ${mapdir1}J719.sorted.dup.reheader.bam ${mapdir1}J750.sorted.dup.reheader.bam ${mapdir1}J757.sorted.dup.reheader.bam
