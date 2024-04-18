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
    popindfile=${DIR}processed/BCFtools/vcftools_filtered/plink_QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.${macval}/pop_ind.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.${macval}.txt
    popenvfile=${DIR}metadata/myna_WGS_env_per_pop.csv  # This is kept the same as different mac filter files have the same individuals anyways
    pairspanelpng=${DIR}processed/baypass/${macval}/e_file_pairs_panel.png
    efilepath=${DIR}processed/baypass/${macval}/e_file.txt
    cfilepath=${DIR}processed/baypass/${macval}/contrastfile.txt
    efolderpath=${DIR}processed/baypass/${macval}/

    Rscript --vanilla 05a_create_contrastfile.R ${metafile} ${popindfile} ${popenvfile} ${pairspanelpng} ${efilepath} ${cfilepath} ${efolderpath}

done
