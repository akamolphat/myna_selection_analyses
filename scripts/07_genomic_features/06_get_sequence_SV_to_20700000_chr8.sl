#!/bin/bash -e
#SBATCH --job-name=06
#SBATCH --output=06_errors/06_%j.out
#SBATCH --error=06_errors/06_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=01:00:00
#SBATCH --mem=32G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4

module purge

ml BEDTools/2.30.0-GCC-11.3.0

VEPfolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/VEP/

cd ${VEPfolder}
mkdir -p chr8_peak_genes/

cd chr8_peak_genes/

bed=chr8_region_after_rep_section_chr8.bed
fastafile=chr8_region_after_rep_section_chr8.fasta

echo -e 'Superscaffold_chr8\t20687525\t20700000' > ${bed}

inputfasta=${VEPfolder}ref/AcTris_vAus2.1.fasta

bedtools getfasta -fi ${inputfasta} -bed ${bed} -fo ${fastafile}


