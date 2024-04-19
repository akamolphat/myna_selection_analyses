#!/bin/bash -e
#SBATCH --job-name=02ai
#SBATCH --output=02_errors/02ai_%j.out
#SBATCH --error=02_errors/02ai_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=2:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613

module purge
ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
ml BCFtools/1.16-GCC-11.3.0
ml tabix/0.2.6-GCCcore-9.2.0
ml R/4.2.1-gimkl-2022a

cd /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/

inputa=variant.WGS.bial.QUAL30.minDP5.snps.vcf.gz
outputa=variant.WGS.bial.QUAL30_DP5-35.snps.vcf.gz
outputb=variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.snps.vcf.gz
folderb=vcftools_sum_stats/
folderc="$(basename ${outputb} .snps.vcf.gz)"/

mkdir -p ${folderb}${folderc}

# Filter DP5-35
minDP=5
maxDP=35

# vcftools --gzvcf ${inputa} --remove-indels --minDP ${minDP} --maxDP ${maxDP} --recode --recode-INFO-all --stdout | bcftools filter -e 'AC==0 || AC=AN' -Ou | bcftools view -c 1 -Oz -o ${outputa}

# Index using bcftools index
# bcftools index ${outputa}

# Filter out singletons and doubletons
singleton="$(basename ${outputa} .snps.vcf.gz)"
# 1. Make list of singletons and doubletons
# vcftools --gzvcf ${outputa} --singletons --out ${singleton}
# 2. Filter out singletons and doubletons. ID is set here as it was originally empty
# vcftools --gzvcf ${outputa} --exclude-positions ${singleton}.singletons --recode --recode-INFO-all --stdout | bcftools annotate --set-id "%CHROM\_%POS\_%REF\_%ALT" -O z -o ${outputb}

# Index using bcftools index
# bcftools index ${outputb}

VCF=${outputb}
outfile="$(basename ${outputb} .snps.vcf.gz)"
outpath=${folderb}${folderc}
OUT=${outpath}${outfile}

# Allele frequency
# vcftools --gzvcf $VCF --freq2 --out $OUT --max-alleles 2
# Mean depth per individual
# vcftools --gzvcf $VCF --depth --out $OUT
# Mean depth per site
# vcftools --gzvcf $VCF --site-mean-depth --out $OUT
# Mean site quality
# vcftools --gzvcf $VCF --site-quality --out $OUT
# Proportion of missing data per individual
# vcftools --gzvcf $VCF --missing-indv --out $OUT
# Proportion of missing data per site
# vcftools --gzvcf $VCF --missing-site --out $OUT
# Heterozygousity and inbreeding coefficient per individual
# vcftools --gzvcf $VCF --het --out $OUT
# Relatedness2
# vcftools --gzvcf $VCF --relatedness2 --out $OUT

# Plot the results

cd $outpath
outputpdf=${outfile}.all.stats.pdf

Rscript --vanilla /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/01_data_exploration_site_indiv_stats.R $outfile $outputpdf

# cd ../../
#
# Rscript --vanilla /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/plot_DP_hist.R ${outputa} 50000 ${outpath}${outfile}.DP_hist.pdf

