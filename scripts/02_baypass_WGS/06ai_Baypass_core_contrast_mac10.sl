#!/bin/bash -e
#SBATCH --job-name=06ai
#SBATCH --output=06_errors/06ai_%j.out
#SBATCH --error=06_errors/06ai_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --mem=6G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4
#SBATCH --array=1-100

module purge
ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
ml BCFtools/1.16-GCC-11.3.0
module load PLINK/1.09b6.16
module load BayPass/2.31-intel-2022a
ml R/4.2.1-gimkl-2022a

seedno=5001
BPDIR=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/
BPDIR1=mac10
runno=$(printf %03d $SLURM_ARRAY_TASK_ID)
BPRUN=subset_${runno}
DIR=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/
cd ${DIR}

namebase=plink_QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10
plinkfolder=${namebase}/
popind=${namebase//plink_/pop_ind.}.txt
popind=${DIR}${plinkfolder}${popind}
SNPtxtfile=${namebase//plink_/variant_calls.bial.}_SNPs.txt
SNPtxt=${DIR}${plinkfolder}${SNPtxtfile}

PLINK=${DIR}${plinkfolder}myna.plink.ped
PLINK2=${DIR}${plinkfolder}myna.plink

bayp=${BPDIR}${BPDIR1}/myna_baypass_${BPDIR1}_${BPRUN}.txt

npop=$(cut -f 1 ${popind} | sort | uniq | wc -l)
npop2=$((2*npop))

cd ${plinkfolder}

nsnp=$(wc -l < "myna.plink.map")
cd ${BPDIR}${BPDIR1}
mkdir -p ${BPRUN}

echo "Run BayPass core model"
i_baypass -seed ${seedno} -npop ${npop} -gfile ${bayp} -contrastfile ${contrastfile} -nthreads 4 -outprefix ${BPRUN}/myna_baypass_${BPDIR1}_${BPRUN}

