#!/bin/bash -e
#SBATCH --job-name=07bi
#SBATCH --output=07_errors/07bi_%j.out
#SBATCH --error=07_errors/07bi_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4

module purge

BPDIR=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/
BPDIR1=mac10
combdir=combined
mkdir -p ${BPDIR}${BPDIR1}/${combdir}

echo ${BPRUN}

cd ${BPDIR}${BPDIR1}

for runno in {1..100}
do

    runno=$(printf %03d ${runno})
    BPRUN=subset_${runno}
    SNPout=${BPDIR}${BPDIR1}/${BPRUN}/myna_baypass_${BPDIR1}_${BPRUN}_SNPs.txt
    SNPrep=${BPDIR}${BPDIR1}/${BPRUN}/myna_baypass_${BPDIR1}_${BPRUN}_SNPs_contrast_rep.txt
    outfile=${BPDIR}${BPDIR1}/${BPRUN}/myna_baypass_${BPDIR1}_${BPRUN}_summary_contrast.out
    tempfileSNP=${BPDIR}${BPDIR1}/${BPRUN}/myna_baypass_${BPDIR1}_${BPRUN}_summary_contrast_SNPs_temp.out
    outfileSNP=${BPDIR}${BPDIR1}/${BPRUN}/myna_baypass_${BPDIR1}_${BPRUN}_summary_contrast_SNPs.out

    cd ${BPRUN}
    echo ${BPRUN}

    nlines=$(wc -l < ${outfile})
    nsnppop=$((nlines-1))
    nlines2=$(wc -l < ${SNPout})
    nsnp=$((nlines2-1))
    ncomp=$((nsnppop/nsnp))

    head -n 1 ${SNPout} > ${SNPrep}
    for ((i = 1; i <= $ncomp; i++))
    do 

        tail -n +2 -q ${SNPout} >> ${SNPrep}

    done

    paste ${SNPrep} ${outfile} > ${tempfileSNP}
    
    tr -s ' ' '\t' < ${tempfileSNP} > ${outfileSNP}

    rm ${tempfileSNP}

    cd ${BPDIR}${BPDIR1}
done

cd ${BPDIR}${BPDIR1}  # cd again to make sure that we are in the right place

echo "Combine files"
pwd

SNPout=${BPDIR}${BPDIR1}/subset_001/myna_baypass_${BPDIR1}_subset_001_SNPs.txt
file1=${BPDIR}${BPDIR1}/subset_001/myna_baypass_${BPDIR1}_subset_001_summary_contrast_SNPs.out

nl=$(wc -l < ${file1})
ns=$((nl-1))
nltwo=$(wc -l < ${SNPout})
nsnpone=$((nltwo-1))
ncmp=$((ns/nsnpone))

echo "Loop through each contrast"
for ((i = 1; i <= $ncmp; i++))
do
     compno=$(printf %03d ${i})
     outcombfile=${BPDIR}${BPDIR1}/${combdir}/myna_baypass_${BPDIR1}_${combdir}_summary_contrast_CON_${compno}_SNPs_int.out
     outcombfile2=${BPDIR}${BPDIR1}/${combdir}/myna_baypass_${BPDIR1}_${combdir}_summary_contrast_CON_${compno}_SNPs.out
     head -n 1 ${file1} > ${outcombfile}     

     cd ${BPDIR}${BPDIR1}
     
     echo "Loop through each subset"

     for runno in {1..100}
     do
     
         runno=$(printf %03d ${runno})
         BPRUN=subset_${runno}
         SNPout=${BPDIR}${BPDIR1}/${BPRUN}/myna_baypass_${BPDIR1}_${BPRUN}_SNPs.txt
         outfileSNP=${BPDIR}${BPDIR1}/${BPRUN}/myna_baypass_${BPDIR1}_${BPRUN}_summary_contrast_SNPs.out

         cd ${BPRUN}
         echo ${BPRUN}
     
         nlines=$(wc -l < ${outfileSNP})
         nsnppop=$((nlines-1))
         nlines2=$(wc -l < ${SNPout})
         nsnp=$((nlines2-1))
         ncomp=$((nsnppop/nsnp))
         nstart=$((((i-1)*nsnp)+1))
         tail -n +2 ${outfileSNP} | tail -n +${nstart} | head -n ${nsnp} >> ${outcombfile}

         cd ${BPDIR}${BPDIR1}

     done

     sort -t $'\t' -nk4,4 ${outcombfile} > ${outcombfile2}

     rm ${outcombfile}

done

