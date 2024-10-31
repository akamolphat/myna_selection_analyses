#!/bin/bash -e
#SBATCH --job-name=09ci
#SBATCH --output=09_errors/09ci_%j.out
#SBATCH --error=09_errors/09ci_%j.err
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

echo "Get genes from XtX outlier regions"

bedin=${bedfolder}myna_baypass_${macval}_combined_xtx_outliers_0.99999_cluster1_outlier_region_merged.bed
genelist=${bedfolder}myna_baypass_${macval}_xtx_outliers_0.99999_cluster1_outlier_region_genes_list.txt
bedtools intersect -wb -a ${bedin} -b ${gffin} | sort | uniq | grep -oP '(?<=Preferred_name=).*?(?=;)' | sort > ${genelist}
keggpathwaylist=${bedfolder}myna_baypass_${macval}_xtx_outliers_0.99999_cluster1_outlier_region_KEGG_Pathway_list.txt
bedtools intersect -wb -a ${bedin} -b ${gffin} | sort | uniq | grep -oP '(?<=em_KEGG_Pathway=).*?(?=;)' | tr ',' '\n' | sort > ${keggpathwaylist}

for conid in {1..7}
do

    echo "Get genes from C2 outlier regions, contrast = "${conid}

    bedin=${bedfolder}myna_baypass_${macval}_combined_IS_C2_GEA_summary_contrast_CON_00${conid}_outliers_0.99999_cluster1_outlier_region_merged.bed
    genelist=${bedfolder}myna_baypass_${macval}_CON_00${conid}_outliers_0.99999_cluster1_outlier_region_genes_list.txt
    bedtools intersect -wb -a ${bedin} -b ${gffin} | sort | uniq | grep -oP '(?<=Preferred_name=).*?(?=;)' | sort > ${genelist}
    keggpathwaylist=${bedfolder}myna_baypass_${macval}_CON_00${conid}_outliers_0.99999_cluster1_outlier_region_KEGG_Pathway_list.txt
    bedtools intersect -wb -a ${bedin} -b ${gffin} | sort | uniq | grep -oP '(?<=em_KEGG_Pathway=).*?(?=;)' | tr ',' '\n' | sort > ${keggpathwaylist}
done
