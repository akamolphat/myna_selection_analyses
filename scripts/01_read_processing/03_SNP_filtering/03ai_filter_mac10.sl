#!/bin/bash -e
#SBATCH --job-name=03ai
#SBATCH --output=03_errors/03ai_%j.out
#SBATCH --error=03_errors/03ai_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=12:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task

module purge
ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
ml BCFtools/1.16-GCC-11.3.0
ml tabix/0.2.6-GCCcore-9.2.0
ml R/4.2.1-gimkl-2022a
ml Stacks/2.61-gimkl-2022a


cd /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/
inputa=variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.snps.vcf.gz
outputa=${inputa//.snps.vcf.gz/.mac10.snps.vcf.gz}

folderb=vcftools_sum_stats/
folderc="$(basename ${outputa} .snps.vcf.gz)"/

mkdir -p ${folderb}${folderc}

bcftools view -c 10:minor ${inputa} -Oz -o ${outputa}
bcftools index ${outputa}

VCF=${outputa}
outfile="$(basename ${outputa} .snps.vcf.gz)"
outpath=${folderb}${folderc}
OUT=${outpath}${outfile}

# Allele frequency
vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2
# Mean depth per individual
vcftools --gzvcf $VCF --depth --out $OUT
# Mean depth per site
vcftools --gzvcf $VCF --site-mean-depth --out $OUT
# Mean site quality
vcftools --gzvcf $VCF --site-quality --out $OUT
# Proportion of missing data per individual
vcftools --gzvcf $VCF --missing-indv --out $OUT
# Proportion of missing data per site
vcftools --gzvcf $VCF --missing-site --out $OUT
# Heterozygousity and inbreeding coefficient per individual
vcftools --gzvcf $VCF --het --out $OUT
# Relatedness2
#vcftools --gzvcf $VCF --relatedness2 --out $OUT

# Plot the results

cd $outpath
outputpdf=${outfile}.all.stats.pdf

#Rscript --vanilla /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/01_data_exploration_site_indiv_stats.R $outfile $outputpdf
