#!/bin/bash -e
#SBATCH --job-name=06bii
#SBATCH --output=06_errors/06bii_%j.out
#SBATCH --error=06_errors/06bii_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --mem=16G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4

module purge
ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
ml BCFtools/1.16-GCC-11.3.0
module load PLINK/1.09b6.16
module load BayPass/2.31-intel-2022a
ml R/4.2.1-gimkl-2022a

seedno=5001
BPDIR=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/
BPDIR1=mac3
runno=001
BPRUN=subset_${runno}

bayp1=${BPDIR}${BPDIR1}/myna_baypass_${BPDIR1}.txt

echo "Simulate baypass"
omega=myna_baypass_${BPDIR1}_${BPRUN}_mat_omega.out
beta=myna_baypass_${BPDIR1}_${BPRUN}_summary_beta_params.out
wdir=${BPDIR}${BPDIR1}/${BPRUN}/
nsnp=20000

cd ${wdir}
echo "Start simulating PODs"

Rscript --vanilla /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/shared_scripts/baypass_simulate1.R ${omega} ${beta} ${bayp1} ${nsnp} ${wdir}

echo "Finish simulating PODs"

bayp=G.btapods

npop2=$(head -n 1 ${bayp} | awk '{print NF}')
npop=$((npop2/2))

echo "Run BayPass core model"
contrastfile=${BPDIR}${BPDIR1}/contrastfile.txt
i_baypass -seed ${seedno} -npop ${npop} -gfile ${bayp} -contrastfile ${contrastfile} -nthreads 4 -outprefix ${bayp}
