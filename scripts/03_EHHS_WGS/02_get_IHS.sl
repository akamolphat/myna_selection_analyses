#!/bin/bash -e
#SBATCH --job-name=02
#SBATCH --output=02_errors/02_%j.out
#SBATCH --error=02_errors/02_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=6:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4
#SBATCH --array=1-35

module purge

ml BCFtools/1.16-GCC-11.3.0
ml R/4.2.1-gimkl-2022a

source ~/.bash_profile

chrno=$((SLURM_ARRAY_TASK_ID-1))
chrls=("chr1" "chr1A" "chr2" "chr3" "chr4" "chr4A" "chr5a" "chr5b" "chr5c" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chr23" "chr24" "chr25" "chr26" "chr27" "chr29" "chr30" "chr31" "chr32")
chr=${chrls[chrno]}
echo ${chr}

vcffolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/
IHSfolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/
IHSresfold=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/EHHS/

invcf=${vcffolder}variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.${chr}.phased_shapeit5.snps.vcf.gz
scanfile=${IHSfolder}${chr}_scanhh.txt
outIHS=${IHSfolder}${chr}_IHS.txt
outpng=${IHSresfold}${chr}_IHS.png

cd /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/

Rscript --vanilla get_IHS.R ${invcf} ${scanfile} ${outIHS} ${outpng}
