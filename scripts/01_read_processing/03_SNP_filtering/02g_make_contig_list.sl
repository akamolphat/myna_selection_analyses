#!/bin/bash -e
#SBATCH --job-name=02g
#SBATCH --output=02_errors/02g_%j.out
#SBATCH --error=02_errors/02g_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task

module purge

reffai=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/align2ref/ref/AcTris_vAus2.1.fasta.fai
outcontiglist=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/Superscaffold.list

grep "Superscaffold_chr" ${reffai} | cut -f 1 | grep -v "Superscaffold_chrZ" | grep -v "Superscaffold_chrW" > ${outcontiglist}
