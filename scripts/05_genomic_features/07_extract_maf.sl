#!/bin/bash -e
#SBATCH --job-name=09a
#SBATCH --output=09_errors/09a_%j.out
#SBATCH --error=09_errors/09a_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=00:10:00
#SBATCH --mem=2G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4

module purge

frqfolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/plink_QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10/
outsnpfolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/mac10/combined/
frqfile=${frqfolder}myna.plink.frq.strat

cd ${outsnpfolder}

#grep "_chr8_" ${frqfile} > ${outsnpfolder}MAF_per_pop_chr8.out
#grep "_chr3_" ${frqfile} > ${outsnpfolder}MAF_per_pop_chr3.out
#grep "_chr22_" ${frqfile} > ${outsnpfolder}MAF_per_pop_chr22.out
grep "_chr2_" ${frqfile} > ${outsnpfolder}MAF_per_pop_chr2.out
