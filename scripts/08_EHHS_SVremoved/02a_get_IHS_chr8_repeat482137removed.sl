#!/bin/bash -e
#SBATCH --job-name=04b
#SBATCH --output=04_errors/04b_%j.out
#SBATCH --error=04_errors/04b_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=1:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4

module purge

ml R/4.2.1-gimkl-2022a

vcffolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/
IHSfolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/
IHSresfold=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/EHHS/

invcf=${vcffolder}variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.SVsnpsremoved.snps.vcf.gz
scanfile=${IHSfolder}chr8_SVsnpsremoved_scanhh.txt
outIHS=${IHSfolder}chr8_SVsnpsremoved_IHS.txt
outpng=${IHSresfold}chr8_SVsnpsremoved_IHS.png

cd /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/

Rscript --vanilla get_IHS.R ${invcf} ${scanfile} ${outIHS} ${outpng}

