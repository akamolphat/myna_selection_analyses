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
ml R/4.2.1-gimkl-2022a

cd /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/scripts/

Rscript --vanilla 02g_make_contig_list.R
