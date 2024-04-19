#!/bin/bash -e
#SBATCH --job-name=01aiii_plot_minDP5
#SBATCH --output=01_errors/01aiii_%j.out
#SBATCH --error=01_errors/01aiii_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=3:00:00
#SBATCH --mem=64G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613

ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
ml BCFtools/1.16-GCC-11.3.0
ml tabix/0.2.6-GCCcore-9.2.0
ml R/4.2.1-gimkl-2022a

cd /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/

inputa=variant.WGS.bial.QUAL30.snps.vcf.gz
outputa=variant.WGS.bial.QUAL30.minDP5.snps.vcf.gz
foldera=vcftools_filtered/
folderb=vcftools_sum_stats/
folderc="$(basename ${outputa} .snps.vcf.gz)"/

#mkdir -p ${foldera}${folderb}${folderc}

cd ${foldera}

#vcftools --gzvcf ${inputa} --remove-indels --minQ 30 --minDP 5 --recode --recode-INFO-all --stdout | bcftools filter -e 'AC==0 || AC=AN' -Ou | bcftools view -c 1 -Oz -o ${outputa} 

#bcftools index ${outputa}

VCF=${outputa}
outfile="$(basename ${outputa} .snps.vcf.gz)"
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

# Rscript --vanilla /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/01_data_exploration_site_indiv_stats.R $outfile $outputpdf

cd ../../

Rscript --vanilla /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/plot_DP_hist.R ${outputa} 50000 ${outpath}${outfile}.DP_hist.pdf
