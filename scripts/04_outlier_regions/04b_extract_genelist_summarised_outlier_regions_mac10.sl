#!/bin/bash -e
#SBATCH --job-name=09fi
#SBATCH --output=09_errors/09fi_%j.out
#SBATCH --error=09_errors/09fi_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=0:30:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4

module purge
ml BEDTools/2.30.0-GCC-11.3.0

ml BEDTools/2.30.0-GCC-11.3.0

macval=mac10

bedfolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/results/combined_outliers/
gfffolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/VEP/gff/
gffin=${gfffolder}Superscaffold_cleaned_VEP.gff.gz

echo "n > 10 outlier regions"

bedin=${bedfolder}outliers_n10_summary.bed
genelist=${bedfolder}outliers_n10_summary_genelist.txt
bedtools intersect -wb -a ${bedin} -b ${gffin} | grep -oP '(?<=Preferred_name=).*?(?=;)' > temp.txt
bedtools intersect -wb -wa -a ${bedin} -b ${gffin} | grep Preferred_name | grep -v "Preferred_name=;" | cut -f 1,2,3,4,5,9,10 > temp2.txt

paste -d "\t" temp2.txt temp.txt > ${genelist}
rm temp.txt temp2.txt 

echo "n > 5 outlier regions"

bedin=${bedfolder}outliers_n5_summary.bed
genelist=${bedfolder}outliers_n5_summary_genelist.txt
bedtools intersect -wb -a ${bedin} -b ${gffin} | grep -oP '(?<=Preferred_name=).*?(?=;)' > temp.txt
bedtools intersect -wb -wa -a ${bedin} -b ${gffin} | grep Preferred_name | grep -v "Preferred_name=;" | cut -f 1,2,3,4,5,9,10 > temp2.txt

paste -d "\t" temp2.txt temp.txt > ${genelist}
rm temp.txt temp2.txt 
