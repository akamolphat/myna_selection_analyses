#!/bin/bash -e
#SBATCH --job-name=07ai
#SBATCH --output=07_errors/07ai_%j.out
#SBATCH --error=07_errors/07ai_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=0:20:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4

module purge

BPDIR=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/
BPDIR1=mac10
combdir=combined
outcombfile=${BPDIR}${BPDIR1}/${combdir}/myna_baypass_${BPDIR1}_${combdir}_summary_pi_xtx_SNPs.out
mkdir -p ${BPDIR}${BPDIR1}/${combdir}

echo ${BPRUN}

cd ${BPDIR}${BPDIR1}

for runno in {1..100}
    do

    runno=$(printf %03d ${runno})
    BPRUN=subset_${runno}
    SNPout=${BPDIR}${BPDIR1}/${BPRUN}/myna_baypass_${BPDIR1}_${BPRUN}_SNPs.txt
    outfile=${BPDIR}${BPDIR1}/${BPRUN}/myna_baypass_${BPDIR1}_${BPRUN}_summary_pi_xtx.out
    outfileSNP=${BPDIR}${BPDIR1}/${BPRUN}/myna_baypass_${BPDIR1}_${BPRUN}_summary_pi_xtx_SNPs.out

    cd ${BPRUN}
    echo ${BPRUN}

    paste ${SNPout} ${outfile} | tr -s ' ' '\t' > ${outfileSNP}

    cd ${BPDIR}${BPDIR1}
    done

cd ${BPDIR}${BPDIR1}  # cd again to make sure that we are in the right place

echo "Combine files"
# echo ${outfile}
pwd
# covariate free
# First file
file1=${BPDIR}${BPDIR1}/subset_001/myna_baypass_${BPDIR1}_subset_001_summary_pi_xtx_SNPs.out

echo "Output header of covariate free XtX file"
head -n 1 ${file1} > ${outcombfile}
echo "Add the tails"
tail -n +2 -q subset_*/myna_baypass_mac10_subset_[0-9][0-9][0-9]_summary_pi_xtx_SNPs.out | sort -t $'\t' -nk4,4 >> ${outcombfile}
#wc -l ${outcombfile}

