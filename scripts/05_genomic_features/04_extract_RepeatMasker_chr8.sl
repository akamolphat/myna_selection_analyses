#!/bin/bash -e
#SBATCH --job-name=04
#SBATCH --output=04_errors/04_%j.out
#SBATCH --error=04_errors/04_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=00:05:00
#SBATCH --mem=8G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4

module purge

ml BEDTools/2.30.0-GCC-11.3.0
ml VEP/107.0-GCC-11.3.0-Perl-5.34.1
ml tabix/0.2.6-GCCcore-9.2.0
ml HTSlib/1.19-GCC-11.3.0
ml BWA/0.7.17-GCC-11.3.0


folder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/RepeatMasker/

cd ${folder}
ln -s /nesi/nobackup/uoa02613/kstuart_projects/At1_MynaGenome/analysis/repeats/repeatmasker/AcTris_vAus2.0repeatlib_hardmask/AcTris_vAus2.0.fasta.out ./

grep "Superscaffold_chr8" AcTris_vAus2.0.fasta.out | awk '{if ($5 == "Superscaffold_chr8" && $6 < 20900000 && $7 > 20400000) print $0}' > Repeats_chr8_20400000_20900000.txt
