#!/bin/bash -e
#SBATCH --job-name=09di
#SBATCH --output=09_errors/09di_%j.out
#SBATCH --error=09_errors/09di_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4

module purge
ml BEDTools/2.30.0-GCC-11.3.0

macval=mac10
bedfolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/${macval}/combined/outliers/
gfffolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/VEP/gff/
gffin=${gfffolder}Superscaffold_cleaned_VEP.gff.gz

echo "Create bed file based on all the SNPs input into BayPass and add/subtract bp_region_dist from end_pos/start_pos"
refsnpfolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/plink_QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.${macval}/

# Need to create bp_region_dist*2 as the end_pos column is based on the start_pos column which were subtracted by bp_region_dist at first
bp_region_dist=250000
bp_region_dist2=$((2*bp_region_dist))
infilerefsnp=${refsnpfolder}variant.baypass_SNPs.txt
outfilerefsnp=${bedfolder}myna_baypass_${macval}_all_snps_region_merged.bed
awk -v OFS="\t" -v var1="$bp_region_dist" -v var2="$bp_region_dist2" '{print $1,$2-=var1,$2+=var2}' ${infilerefsnp} | awk -v OFS='\t' '$2<0{$2=0}1' | less -S | bedtools merge > ${outfilerefsnp}

echo "Extract genes within the reference SNP region"
genelist=${bedfolder}myna_baypass_${macval}_all_snps_region_genes_list.txt
bedtools intersect -wb -a ${outfilerefsnp} -b ${gffin} | sort | uniq | grep -oP '(?<=Preferred_name=).*?(?=;)' | sort > ${genelist}
keggpathwaylist=${bedfolder}myna_baypass_${macval}_all_snps_region_KEGG_Pathway_list.txt
bedtools intersect -wb -a ${outfilerefsnp} -b ${gffin} | sort | uniq | grep -oP '(?<=em_KEGG_Pathway=).*?(?=;)' | tr ',' '\n' | sort > ${keggpathwaylist}
