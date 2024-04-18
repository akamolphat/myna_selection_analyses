#!/bin/bash -e
#SBATCH --job-name=03
#SBATCH --output=03_errors/03_%j.out
#SBATCH --error=03_errors/03_%j.err
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


TEfolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/TE/

cd ${TEfolder}
ln -s /nesi/nobackup/uoa02613/kstuart_projects/At3_TEpopgen/data/mcclintock/mcclintock_out_genorescue2/mcclintock_te_sample_variants.vcf

grep Superscaffold_chr8 mcclintock_te_sample_variants.vcf | awk '{if ($2 > 20400000 && $2 < 20900000) print $0}' > TE_chr8_20400000_20900000.txt
