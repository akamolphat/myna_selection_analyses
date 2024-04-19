#!/bin/bash -e
#SBATCH --job-name=02c
#SBATCH --output=02_errors/02c_%j.out
#SBATCH --error=02_errors/02c_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=4:00:00
#SBATCH --mem=2G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613

ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
ml BCFtools/1.16-GCC-11.3.0
ml tabix/0.2.6-GCCcore-9.2.0
ml R/4.2.1-gimkl-2022a

cd /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/

Rscript --vanilla scripts/02b_make_popmap_sample_list.R 

cd /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/

inputa=variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.snps.vcf.gz
outputa=variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.snps.vcf.gz
keepsampletxt=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/filtsamp_2keep.txt
foldera=vcftools_filtered/
folderb=vcftools_sum_stats/
folderc="$(basename ${outputa} .snps.vcf.gz)"/

mkdir -p ${foldera}${folderb}${folderc}

cd ${foldera}

# Keep samples that are on the list
bcftools view -S ${keepsampletxt} --force-samples ${inputa} | bcftools filter -e 'AC==0 || AC=AN' -Ou | bcftools view -c 1 -Oz -o ${outputa}

bcftools index ${outputa}

# Only two sample is removed so I am not replotting the results
