#!/bin/bash -e
#SBATCH --job-name=09ai
#SBATCH --output=09_errors/09ai_%j.out
#SBATCH --error=09_errors/09ai_%j.err
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

echo "Get M_XtX threshold"

path2script=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/baypass_XtX_threshold.R
macval=mac10
PODxtx=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/${macval}/subset_001/G.btapods_summary_pi_xtx.out
quant_thres=0.99999
thres=$(Rscript --vanilla ${path2script} ${PODxtx} ${quant_thres})

echo "Subset for XtX outliers only"
wdir=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/${macval}/combined
cd ${wdir}
mkdir -p outliers
XtXfile=myna_baypass_${macval}_combined_summary_pi_xtx_SNPs.out
XtXoutlierfile=outliers/${XtXfile//_summary_pi_xtx_SNPs.out/_xtx_outliers}_${quant_thres}.txt

# Subset table myna_*_pi_xtx.out based on $thres
#head -n 1 ${XtXfile} > ${XtXoutlierfile}
awk -v var="$thres" '$8>var' ${XtXfile} > ${XtXoutlierfile}


for conid in {1..7}
do

    echo "Get M_C2 threshold for contrast = "${conid}

    path2script=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/baypass_C2_threshold.R
    PODC2=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/${macval}/subset_001/G.btapods_IS_C2_GEA_summary_contrast.out
    quant_thres=0.99999
    thres=$(Rscript --vanilla ${path2script} ${PODC2} M_C2 ${conid} ${quant_thres})
    
    echo ${thres}

    echo "Subset for C2 outliers only, for contrast = "${conid}
    C2file=myna_baypass_${macval}_combined_IS_C2_GEA_summary_contrast_CON_00${conid}_SNPs.out
    C2outlierfile=outliers/${C2file//_SNPs.out/_outliers}_${quant_thres}.txt

    # Subset table myna_*_pi_xtx.out based on $thres
    # head -n 1 ${C2file} > ${C2outlierfile}
    awk -v var="$thres" '$7>var' ${C2file} > ${C2outlierfile}
done

wdir=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/
cd ${wdir}
mkdir IHS_outliers
IHSfile=WGS_IHS.txt
IHSoutlierfile=IHS_outliers/WGS_IHS_outliers_logpval_6.txt
awk '$4>6' ${IHSfile} > ${IHSoutlierfile}
# Rename key column names to match XtX and C2 for compatability
sed -e '1s/POSITION/POS/' -e '1s/CHR/chr/' -i ${IHSoutlierfile}

wdir=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/EHHS/
cd ${wdir}
mkdir IHS_outliers
IHSfile=WGS_IHS_chr8_SVremoved.txt
IHSoutlierfile=IHS_outliers/WGS_IHS_chr8_SVremoved_outliers_logpval_6.txt
awk '$4>6' ${IHSfile} > ${IHSoutlierfile}
# Rename key column names to match XtX and C2 for compatability
sed -e '1s/POSITION/POS/' -e '1s/CHR/chr' -i ${IHSoutlierfile}




