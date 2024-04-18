#!/bin/bash -e
#SBATCH --job-name=05
#SBATCH --output=05_errors/05_%j.out
#SBATCH --error=05_errors/05_%j.err
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

DIR=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/

for macval in mac3 mac10
do

    metafile=${DIR}metadata/myna_WGS_meta_pop_env.csv
    efolderpath=${DIR}processed/baypass/${macval}/

    Rscript --vanilla 05a_create_contrastfile.R ${metafile} ${cfilepath}

done
