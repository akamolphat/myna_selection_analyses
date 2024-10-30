#!/bin/bash -e
#SBATCH --job-name=11
#SBATCH --output=11_errors/11_%j.out
#SBATCH --error=11_errors/11_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=4:00:00
#SBATCH --mem=40G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task

module purge
ml R/4.2.1-gimkl-2022a

cd /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/scripts/

Rscript --vanilla 02_Fig2_manhattan_plots_XtX_CON001.R
