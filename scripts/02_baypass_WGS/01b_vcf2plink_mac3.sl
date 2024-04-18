#!/bin/bash -e
#SBATCH --job-name=04aii
#SBATCH --output=04_errors/04aii_%j.out
#SBATCH --error=04_errors/04aii_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task

module purge
ml VCFtools/0.1.15-GCC-9.2.0-Perl-5.30.1
ml BCFtools/1.16-GCC-11.3.0
module load PLINK/1.09b6.16

DIR=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/
cd ${DIR}

inputa=variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac3.snps.vcf.gz
outputa=myna.plink
outputb=myna.plink
outputSNPtxt=${inputa//.snps.vcf.gz/_SNPs.txt}

pfold=${inputa//variant.WGS.bial./plink_}
folderplink="$(basename ${pfold} .snps.vcf.gz)"/

mkdir -p ${folderplink}

# Convert VCF to PLINK
vcftools --gzvcf ${inputa} --plink --out ${folderplink}${outputa}

cd ${folderplink}
# Convert PLINK to BED
plink --file ${outputa} --make-bed --noweb --out ${outputb}

zgrep -v "^#" ${DIR}${inputa} | cut -f1-3 | awk '{print $0"\t"NR}' > ${outputSNPtxt}

# Update PED file as well with the correct populations

popmap1=${DIR}pop_map_WGS_ALL.txt
popmaptxt=${inputa//.snps.vcf.gz/.txt}
popmaptxt=${popmaptxt//variant.WGS.bial/pop_map}
popindtxt=${popmaptxt//pop_map/pop_ind}
popmap2=${DIR}${folderplink}${popmaptxt}
popind=${DIR}${folderplink}${popindtxt}

PLINK=${DIR}${folderplink}myna.plink.ped
PLINK2=${DIR}${folderplink}myna.plink

echo "update plink.ped file"
# Extract sample ID from plink
cut -f 2 $PLINK > x0.delete

cat x0.delete | tr -d "[:blank:]" > x1.delete

# subset pop_map_ALL with the samples in PLINK
grep -Fwf x1.delete ${popmap1} > ${popmap2}

# Order the samples in $popmap2 by x1.delete
awk '
    FNR==NR { a[NR]=$1; next }
    { b[$1] = (b[$1] ? b[$1]"\n" $0 : $0) }
    END{ for (i=1; i<=length(a); i++) print b[a[i]] }
' x1.delete ${popmap2} | awk '{print $2,"\t",$1}' > ${popind}

#remove first 2 columns from PLINK file
cut -f 3- $PLINK > x2.delete

# Replace first 2 columns from PLINK file with ${popind}
paste ${popind} x2.delete > ${PLINK}
rm x0.delete x1.delete x2.delete

npop=$(cut -f 1 ${popind} | sort | uniq | wc -l)
npop2=$((2*npop))

echo complete updating plink.ped file

