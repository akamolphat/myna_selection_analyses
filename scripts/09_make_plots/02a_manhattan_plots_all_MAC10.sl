#!/bin/bash -e
#SBATCH --job-name=11bi
#SBATCH --output=11_errors/11bi_%j.out
#SBATCH --error=11_errors/11bi_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=4:00:00
#SBATCH --mem=60G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --partition=milan

module purge
ml R/4.2.1-gimkl-2022a

cd /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/scripts/

Rscript --vanilla 11_manhattan_plots_all.R
