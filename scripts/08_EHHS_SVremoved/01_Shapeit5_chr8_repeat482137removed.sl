#!/bin/bash -e
#SBATCH --job-name=04a
#SBATCH --output=04_errors/04a_%j.out
#SBATCH --error=04_errors/04a_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=6:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4

module purge

ml BCFtools/1.16-GCC-11.3.0

source ~/.bash_profile

vcffolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/

cd ${vcffolder}

bcftools view -r Superscaffold_chr8:0-20678726,Superscaffold_chr8:20687525-32000000 -Oz -o variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.SVsnpsremoved.snps.vcf.gz variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.snps.vcf.gz
bcftools index variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.SVsnpsremoved.snps.vcf.gz

invcf=variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.SVsnpsremoved.snps.vcf.gz
outbcf=${invcf//snps.vcf.gz/phased_shapeit5.snps.bcf}

phase_common --input ${invcf} --region "Superscaffold_chr8" --output ${outbcf}

outvcf=${outbcf//.snps.bcf/.snps.vcf.gz}
bcftools view -Oz -o ${outvcf} ${outbcf}
bcftools index ${outvcf} 
