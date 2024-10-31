#!/bin/bash -e
#SBATCH --job-name=WGS_5kbthin
#SBATCH --output=WGS_5kbthin_%j.out
#SBATCH --error=WGS_5kbthin_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=END
#SBATCH --time=6:00:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613

module purge
ml BCFtools/1.16-GCC-11.3.0
ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1

cd /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/

vcfa=variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10.snps.vcf.gz
vcfthin=variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10.5kbthinned.snps.vcf.gz

vcftools --gzvcf ${vcfa} --thin 5000 --recode --recode-INFO-all --stdout | bcftools view -Oz -o ${vcfthin}
bcftools index ${vcfthin}

# Extract key outlier region SNPs only
invcf=variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.snps.vcf.gz
outvcf=variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.mac10.chr8.20625-20699kb.snps.vcf.gz
bcftools --gzvcf ${invcf} -r Superscaffold_chr8:20625000-20699000 -Oz -o ${outvcf}
bcftools index ${outvcf}
# ${outvcf} refers to the path to the output VCF file containing only SNPs from the specified region.

