#!/bin/bash -e
#SBATCH --job-name=01
#SBATCH --output=01_errors/01_%j.out
#SBATCH --error=01_errors/01_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=06:00:00
#SBATCH --mem=50G
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

VEPfolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/VEP/
# Create softlinks to reference genome in Kat's folder -------------#####
mkdir -p ${VEPfolder}ref/
cd ${VEPfolder}ref
#ln -s /nesi/nobackup/uoa02613/kstuart_projects/At1_MynaGenome/analysis/curation/step4_scaffolding/ragtag_atris_synteny/renamed/AcTris_vAus2.1.fasta ./
refin=AcTris_vAus2.1.fasta
#bwa index ${refin}
#bgzip ${refin}

# Create softlinks to gff file in Kat's folder ---------------------#####
mkdir -p ${VEPfolder}gff/
GFFKat=/nesi/nobackup/uoa02613/kstuart_projects/At1_MynaGenome/annotation/eggnog/run1b/braker_gemoma_combined_katmanual_brakerforce_longestIsoform_transcripts.emapper.decorated.gff
cd ${VEPfolder}gff/
#ln -s ${GFFKat} ./

# Clean up GFF file ------------------------------------------------#####
GFF=braker_gemoma_combined_katmanual_brakerforce_longestIsoform_transcripts.emapper.decorated.gff
GFFout=Superscaffold_cleaned_VEP.gff.gz

## Subset GFF file for just the Superscaffolds and autosomes only, then add biotype=protein_coding to every transcript row, and then sort by chromosome 
#grep -v "#" ${GFF} | grep "Superscaffold" | grep -v "Superscaffold_chrW" | grep -v "Superscaffold_chrZ" | awk 'BEGIN{OFS=FS="\t"} $3=="transcript" {$9=$9";biotype=protein_coding"} {print}' | sort -k1,1 -k4,4n -k5,5n -t$'\t' | bgzip -c > ${GFFout}
#tabix -p gff ${GFFout}

# Clean up VCF file -------------------------------------------------#####
mkdir -p ${VEPfolder}vcf/
cd ${VEPfolder}vcf/
invcf=variant.WGS.bial.QUAL30_DP5-35.nosingledoubletons.filtind.5n.contigfilt1.mac10.snps.vcf.gz
#ln -s /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/BCFtools/vcftools_filtered/${invcf} ./

#tabix -h ${invcf} Superscaffold_chr8:20400000-20900000 > chr8_20400000_20900000.vcf 

# Run VEP -----------------------------------------------------------#####
cd ${VEPfolder}
## Full dataset -----------------------------------------------------
vep -i vcf/${invcf} -gff gff/${GFFout} -fasta ref/${refin}.gz -o output/Superscaffold_VEP_out

vep -i vcf/chr8_20400000_20900000.vcf -gff gff/${GFFout} -fasta ref/${refin}.gz -o output/chr8_20400000_20900000_VEP_out

grep missense output/chr8_20400000_20900000_VEP_out | cut -f 1 > output/chr8_20400000_20900000_VEP_out_missense.txt
