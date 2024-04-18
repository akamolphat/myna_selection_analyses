#!/bin/bash -e
#SBATCH --job-name=05ci
#SBATCH --output=05_errors/05ci_%j_%a.out
#SBATCH --error=05_errors/05ci_%j_%a.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=0:05:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4
#SBATCH --array=1-100

module purge

maxarrno=100
SNPDIR=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/
plinkfolder=plink_QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10/
SNPDIR=${SNPDIR}${plinkfolder}
BPDIR=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/
BPDIR1=mac10
runno=$(printf %03d $SLURM_ARRAY_TASK_ID)
if [[ $SLURM_ARRAY_TASK_ID -eq ${maxarrno} ]]
then
  divno=0
else
  divno=$SLURM_ARRAY_TASK_ID
fi

BPRUN=subset_${runno}

echo ${BPRUN}

SNPin=${SNPDIR}variant.baypass_SNPs.txt
cd ${BPDIR}${BPDIR1}

SNPout=${BPDIR}${BPDIR1}/${BPRUN}/myna_baypass_${BPDIR1}_${BPRUN}_SNPs.txt
echo $'chr\tPOS\tID\tMRKALL' > ${SNPout}
awk "NR % ${maxarrno} == $divno" ${SNPin} >> ${SNPout} 

