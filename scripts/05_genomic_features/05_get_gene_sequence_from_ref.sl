#!/bin/bash -e
#SBATCH --job-name=05
#SBATCH --output=05_errors/05_%j.out
#SBATCH --error=05_errors/05_%j.err
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
mkdir chr8_peak_genes/

cd gff/

zcat Superscaffold_cleaned_VEP.gff.gz | grep Superscaffold_chr8 | awk '{if($3 == "gene" && $4 < 20750000 && $4 > 20600000) print $0}' | cut -f 1,4,5 > ${VEPfolder}chr8_peak_genes/chr8_peak_genes.bed

cd ${VEPfolder}ref/
zcat AcTris_vAus2.1.fasta.gz > AcTris_vAus2.1.fasta

cd ${VEPfolder}chr8_peak_genes/

inputfasta=${VEPfolder}ref/AcTris_vAus2.1.fasta
bedtools getfasta -fi ${inputfasta} -bed ${VEPfolder}chr8_peak_genes/chr8_peak_genes.bed -fo ${VEPfolder}chr8_peak_genes/chr8_peak_genes.fasta


