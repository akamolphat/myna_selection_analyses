#!/bin/bash -e
#SBATCH --job-name=09bii
#SBATCH --output=09_errors/09bii_%j.out
#SBATCH --error=09_errors/09bii_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4

module purge
ml R/4.2.1-gimkl-2022a
ml BEDTools/2.30.0-GCC-11.3.0

echo "Get XtX outlier regions"

path2script=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/scripts/09b_define_outlier_regions.R
macval=mac3
wdir=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/${macval}/combined/outliers/
cd ${wdir}
infile=myna_baypass_${macval}_combined_xtx_outliers_0.99999.txt
bp_dist=100000  
bp_region_dist=250000 
outclusterfile=${infile//.txt/_cluster1.txt}
outsummaryfile=${infile//.txt/_cluster1_summary.txt}
outbedfile=${infile//.txt/_cluster1_outlier_region.bed}
outbedmerged=${outbedfile//.bed/_merged.bed}

#echo ${outclusterfile} ${outsummaryfile} ${outbedfile}

Rscript --vanilla ${path2script} ${infile} ${bp_dist} ${bp_region_dist} ${outclusterfile} ${outsummaryfile} ${outbedfile}

bedtools sort -i ${outbedfile} | bedtools merge -i stdin > ${outbedmerged}

rm ${outbedfile}

for conid in {1..7}
do

    echo "Get C2 outlier regions, contrast = "${conid}

    infile=myna_baypass_${macval}_combined_IS_C2_GEA_summary_contrast_CON_00${conid}_outliers_0.99999.txt
    outclusterfile=${infile//.txt/_cluster.txt}
    outsummaryfile=${infile//.txt/_cluster_summary.txt}
    outbedfile=${infile//.txt/_cluster1_outlier_region.bed}
    outbedmerged=${outbedfile//.bed/_merged.bed}
    Rscript --vanilla ${path2script} ${infile} ${bp_dist} ${bp_region_dist} ${outclusterfile} ${outsummaryfile} ${outbedfile}
    bedtools sort -i ${outbedfile} | bedtools merge -i stdin > ${outbedmerged}
    rm ${outbedfile}
done
