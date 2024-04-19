#!/bin/bash -e
#SBATCH --job-name=01
#SBATCH --output=01_errors/01_%j.out
#SBATCH --error=01_errors/01_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=6:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4
#SBATCH --array=1-35

module purge

ml BCFtools/1.16-GCC-11.3.0

source ~/.bash_profile

chrno=$((SLURM_ARRAY_TASK_ID-1))
chrls=("chr1" "chr1A" "chr2" "chr3" "chr4" "chr4A" "chr5a" "chr5b" "chr5c" "chr6" "chr7" "chr8" "chr9" "chr10" "chr11" "chr12" "chr13" "chr14" "chr15" "chr16" "chr17" "chr18" "chr19" "chr20" "chr21" "chr22" "chr23" "chr24" "chr25" "chr26" "chr27" "chr29" "chr30" "chr31" "chr32")
chr=${chrls[chrno]}
echo ${chr}
vcffolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/

cd ${vcffolder}

invcf=variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10.snps.vcf.gz
invcfchr8=${invcf//contigfilt1.mac10.snps.vcf.gz/mac10.${chr}.snps.vcf.gz}
outbcf=${invcfchr8//snps.vcf.gz/phased_shapeit5.snps.bcf}

chrfull=Superscaffold_${chr}

echo ${chrfull}

bcftools view -r ${chrfull} -Oz -o ${invcfchr8} ${invcf}
bcftools index ${invcfchr8}

phase_common --input ${invcfchr8} --region ${chrfull} --output ${outbcf}

outvcf=${invcfchr8//snps.vcf.gz/phased_shapeit5.snps.vcf.gz}
bcftools view -Oz -o ${outvcf} ${outbcf}
bcftools index ${outvcf} 

#ml SHAPEIT4/4.2.0-GCC-9.2.0
#shapeit4 --input ${invcfchr8} --region ${chrfull} --output ${outvcf}
