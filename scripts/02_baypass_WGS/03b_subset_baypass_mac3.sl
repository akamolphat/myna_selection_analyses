#!/bin/bash -e
#SBATCH --job-name=05bii
#SBATCH --output=05_errors/05bii_%j_%a.out
#SBATCH --error=05_errors/05bii_%j_%a.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=00:10:00
#SBATCH --mem=1G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4
#SBATCH --array=1-170

module purge

maxarrno=170
BPDIR=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/
BPDIR1=mac3
runno=$(printf %03d $SLURM_ARRAY_TASK_ID)
if [[ $SLURM_ARRAY_TASK_ID -eq ${maxarrno} ]]
then
  divno=0
else
  divno=$SLURM_ARRAY_TASK_ID
fi

BPRUN=subset_${runno}

echo ${BPRUN}

bayp=${BPDIR}${BPDIR1}/myna_baypass_${BPDIR1}.txt

cd ${BPDIR}${BPDIR1}

baypout=${BPDIR}${BPDIR1}/myna_baypass_${BPDIR1}_${BPRUN}.txt
awk "NR % ${maxarrno} == $divno" ${bayp} > ${baypout} 

mkdir -p ${BPDIR}${BPDIR1}/${BPRUN}
