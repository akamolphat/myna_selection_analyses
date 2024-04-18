#!/bin/bash -e
#SBATCH --job-name=05aii
#SBATCH --output=05_errors/05aii_%j.out
#SBATCH --error=05_errors/05aii_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=02:00:00
#SBATCH --mem=32G
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

BPDIR=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/baypass/
BPDIR1=mac3
DIR=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/
cd ${DIR}

mkdir -p ${BPDIR}${BPDIR1}

namebase=plink_QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac3
plinkfolder=${namebase}/
popind=${namebase//plink_/pop_ind.}.txt
popind=${DIR}${plinkfolder}${popind}
SNPtxtfile=${namebase//plink_/variant.WGS.bial.}_SNPs.txt
SNPtxt=${DIR}${plinkfolder}${SNPtxtfile}

PLINK=${DIR}${plinkfolder}myna.plink.ped
PLINK2=${DIR}${plinkfolder}myna.plink

bayp=${BPDIR}${BPDIR1}/myna_baypass_${BPDIR1}.txt

npop=$(cut -f 1 ${popind} | sort | uniq | wc -l)
npop2=$((2*npop))

cd ${plinkfolder}

echo "calculate .frq from plink file"

plink --file ${PLINK2} --allow-extra-chr --freq counts --family --out myna.plink

echo "Finish calculating .frq"

echo "Reorder .frq.strat file"

awk 'FNR == NR { lineno[$3] = NR; next} {print lineno[$2], $0;}' ${SNPtxt} myna.plink.frq.strat | sort -k 1,1n | cut -d' ' -f2- > myna.plink.frq.strat.reordered

echo "Create baypass file"

tail -n +2 myna.plink.frq.strat.reordered | awk '{ $9 = $8 - $7 } 1' | awk '{print $7,$9}' | tr " " "\n" | awk -v pp=${npop2} '{if (NR % pp == 0){a=a $0"";print a; a=""} else a=a $0" "}' > ${bayp}

echo "Create a file with SNP list, generated from same reordered file"

awk -v npop=${npop} "NR % npop == 2" myna.plink.frq.strat.reordered > baypass_SNPs.txt

echo "Create file with SNP list with CHR and POS"

awk 'FNR == NR { lineno[$2] = NR; next} {print lineno[$3], $0;}' baypass_SNPs.txt ${SNPtxt} | sort -k 1,1n | cut -d' ' -f2- | cut -f1-3 | awk '{print $0"\t"NR}' > variant.baypass_SNPs.txt

