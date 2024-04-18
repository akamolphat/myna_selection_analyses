#!/bin/bash -e
#SBATCH --job-name=02
#SBATCH --output=02_errors/02_%j.out
#SBATCH --error=02_errors/02_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --mem=6G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4

module purge
module load BayPass/2.31-intel-2022a



BPDIR=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/
BPDIR1=mac10
BPRUN=SV_only
seedno=5001
npop=11
bayp=myna_baypass_SV.txt
omega=subset_001/myna_baypass_${BPDIR1}_subset_001_mat_omega.out
contrastfile=${BPDIR}mac10/contrastfile.txt
cd ${BPDIR}${BPDIR1}
mkdir -p ${BPRUN}

i_baypass -seed ${seedno} -npop ${npop} -gfile ${bayp} -contrastfile ${contrastfile} -nthreads 4 -omegafile ${omega} -outprefix ${BPRUN}/myna_baypass_SV_core_C2

