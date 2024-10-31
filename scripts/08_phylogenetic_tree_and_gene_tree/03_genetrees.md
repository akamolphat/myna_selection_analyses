# 03_genetrees.md

These commands were performed interactively so I have put it in an .md file

## Create directories to store outputs
```
cd /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/
mkdir AMY2A_sequences
```

## Extract exon positions of the two AMY2A copies and create bed file
``` 
cd /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/AMY2A_sequences

inputgffgz=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/VEP/gff/Superscaffold_cleaned_VEP.gff.gz
amy2acopy1bed="myna_AMY2A_copy1_exon.bed"

zgrep g_12137 $inputgffgz| grep exon | awk -v OFS='\t' '{print $1, $4-1, $5}' > ${amy2acopy1bed}

amy2acopy2bed="myna_AMY2A_copy2_exon.bed"
zgrep g_12138 $inputgffgz | grep exon | awk -v OFS='\t' '{print $1, $4-1, $5}' > ${amy2acopy2bed}

cd ../../AMY2A_sequences
```
## Extract exon sequences of the two AMY2A copies
```
ml BEDTools/2.30.0-GCC-11.3.0

VEPfolder=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/VEP/

fasta1file=myna_AMY2A_copy1_exon.fasta
amy2acopy1bed="myna_AMY2A_copy1_exon.bed"
amy2acopy2bed="myna_AMY2A_copy2_exon.bed"
inputfasta=${VEPfolder}ref/AcTris_vAus2.1.fasta

bedtools getfasta -fi ${inputfasta} -bed ${amy2acopy1bed} -fo ${fasta1file}

fasta2file=myna_AMY2A_copy2_exon.fasta
bedtools getfasta -fi ${inputfasta} -bed ${amy2acopy2bed} -fo ${fasta2file}
```

## Concatenate nucleotide sequences to get full gene sequence, perhaps for later use
```
grep -v ">" myna_AMY2A_copy1_exon.fasta | tr --delete '\n' > myna_AMY2A_copy1.fasta

grep -v ">" myna_AMY2A_copy2_exon.fasta | tr --delete '\n' > myna_AMY2A_copy2.fasta
# Create a multi-sequence fasta file for later use
echo ">myna_AMY2A_copy1" > myna_AMY2A.fasta
cat myna_AMY2A_copy1.fasta >> myna_AMY2A.fasta
echo >> myna_AMY2A.fasta
echo ">myna_AMY2A_copy2" >> myna_AMY2A.fasta
cat myna_AMY2A_copy2.fasta >> myna_AMY2A.fasta
```

Search of both sequences yield a similar AMY2A-like gene sequence in the common starling. These were the following accession number:
* [XM_014885285.1](https://www.ncbi.nlm.nih.gov/nuccore/XM_014885285.1/). The exons are to be downloaded manually from the [graphical viewer of the locus](https://www.ncbi.nlm.nih.gov/gene/106858541) or from the following [link](https://www.ncbi.nlm.nih.gov/projects/sviewer/sequence.cgi?netcache=0&id=XM_014885285.1&format=fasta&filename=XM_014885285.1.exons.fa&ranges=0-185,186-332,333-530,531-761,762-895,896-1018,1019-1118,1119-1240,1241-1366,1367-1680)
* [XM_014885297.1](https://www.ncbi.nlm.nih.gov/nuccore/XM_014885297.1/). The exons are to be downloaded manually from the [graphical viewer of the locus](https://www.ncbi.nlm.nih.gov/gene?cmd=retrieve&list_uids=106858550) or from the following [link](https://www.ncbi.nlm.nih.gov/projects/sviewer/sequence.cgi?netcache=0&id=XM_014885297.1&format=fasta&filename=XM_014885297.1.exons.fa&ranges=0-212,213-359,360-557,558-788,789-922,923-1045,1046-1145,1146-1267,1268-1393,1394-2090)

These genes are found in an un-named contig in the [North American common starling genome genbank accession GCF_001447265.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_001447265.1/).

The exon sequences of XM_014885285.1 and XM_014885297.1 were downloaded, and the exon sequences were blasted agains the [Australian common starling genome genbank accession GCA_023376015.1](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_023376015.1/) as follows:
```
queryfile=XM_014885285.1.exons.fa
subjectfile=/nesi/nobackup/uoa02613/kstuart_projects/At4_MynaStarling/data/resources/genomes/Svulgaris_vAU_1.0.fasta
blastn -query ${queryfile} -subject ${subjectfile} -outfmt 7 > XM_014885285_exon_blast_SvAU.txt

queryfile=XM_014885297.1.exons.fa
subjectfile=/nesi/nobackup/uoa02613/kstuart_projects/At4_MynaStarling/data/resources/genomes/Svulgaris_vAU_1.0.fasta
blastn -query ${queryfile} -subject ${subjectfile} -outfmt 7 > XM_014885297_exon_blast_SvAU.txt
```
The same was done against the [superb starling genome genbank accession GCA_015883425.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_015883425.2/):
```
# Align starling exons to superb starling
queryfile=XM_014885297.1.exons.fa
subjectfile=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/raw_data/Superb_starling/ncbi_dataset/data/GCA_015883425.2/GCA_015883425.2_CU_Lasu_v2_genomic.fna
blastn -query ${queryfile} -subject ${subjectfile} -outfmt 7 > XM_014885297_exon_blast_Lsuperbus.txt

queryfile=XM_014885285.1.exons.fa
subjectfile=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/raw_data/Superb_starling/ncbi_dataset/data/GCA_015883425.2/GCA_015883425.2_CU_Lasu_v2_genomic.fna
blastn -query ${queryfile} -subject ${subjectfile} -outfmt 7 > XM_014885285_exon_blast_Lsuperbus.txt
```
## Extract exon positions of the two AMY2A copies from the Australian common starling reference genome and create bed file

```
# Copy 1
grep -v "#" XM_014885297_exon_blast_SvAU.txt | sed -n '1p;3p;4p;5p;6p;8p;10p;12p;14p;16p' | awk -v OFS='\t' '{print $2, $9-1, $10}' > SvAU_AMY2A_copy1_exon.bed

# Copy 2
grep -v "#" XM_014885285_exon_blast_SvAU.txt | sed -n '1p;3p;4p;5p;6p;8p;10p;12p;14p;16p' | awk -v OFS='\t' '{print $2, $9-1, $10}' > SvAU_AMY2A_copy2_exon.bed
```
## Extract sequences from SvAU reference genome
```
reffile=/nesi/nobackup/uoa02613/kstuart_projects/At4_MynaStarling/data/resources/genomes/Svulgaris_vAU_1.0.fasta
# Copy 1
bedtools getfasta -fi ${reffile} -bed SvAU_AMY2A_copy1_exon.bed | grep -v ">" | tr --delete '\n' > SvAU_AMY2A_copy1.fasta
# Copy 2
bedtools getfasta -fi ${reffile} -bed SvAU_AMY2A_copy2_exon.bed | grep -v ">" | tr --delete '\n' > SvAU_AMY2A_copy2.fasta

# Concatenate exon sequences into one
# Copy 1
grep -v ">" XM_014885297.1.exons.fa | tr --delete '\n' > SvNA_AMY2A_copy1.fasta
# Copy 2
grep -v ">" XM_014885285.1.exons.fa | tr --delete '\n' > SvNA_AMY2A_copy2.fasta
```
##  Extract sequences of the two AMY2A copies from the superb starling reference genome
1. Extract exon positions of the two AMY2A copies from the superb starling reference genome and create bed file
2. Extract exon sequences 
3. Concatenate sequences
```
# Extract AMY2A Sequences from L. Superbus genome
grep -v "#" XM_014885297_exon_blast_Lsuperbus.txt | awk -v OFS='\t' '$10 >=11087272 && $10 <= 11093392' | awk -v OFS='\t' '{print $2, $10-1, $9}' > Lsuperbus_AMY2A_copy1_exon.bed

grep -v "#" XM_014885285_exon_blast_Lsuperbus.txt | awk -v OFS='\t' '$10 >= 11079747 && $10 <= 11083552' | awk -v OFS='\t' '{print $2, $10-1, $9}' > Lsuperbus_AMY2A_copy2_exon.bed

# Extract sequences from L. Superbus ref genome
reffile=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/raw_data/Superb_starling/ncbi_dataset/data/GCA_015883425.2/GCA_015883425.2_CU_Lasu_v2_genomic.fna
# Copy 1 # tac command reverse the rows and then concatenate. The strands are then complemented and reversed
bedtools getfasta -fi ${reffile} -bed Lsuperbus_AMY2A_copy1_exon.bed| grep -v ">" | tac | tr --delete '\n' | tr ACGTacgt TGCATGCA | rev > Lsuperbus_AMY2A_copy1.fasta
# Copy 2
bedtools getfasta -fi ${reffile} -bed Lsuperbus_AMY2A_copy2_exon.bed| grep -v ">" | tac | tr --delete '\n' | tr ACGTacgt TGCATGCA | rev > Lsuperbus_AMY2A_copy2.fasta
```
## Create single multiple sequence fasta file for MAFFT alignment
```
outfasta=AMY2A_sequences.fasta
echo ">A_tristis_AMY2A_copy1" > ${outfasta}
cat myna_AMY2A_copy1.fasta >> ${outfasta}
echo >> ${outfasta}
echo ">A_tristis_AMY2A_copy2" >> ${outfasta}
cat myna_AMY2A_copy2.fasta >> ${outfasta}
echo >> ${outfasta}
echo ">S_vulgaris_NA_AMY2A_copy1" >> ${outfasta}
cat SvNA_AMY2A_copy1.fasta >> ${outfasta}
echo >> ${outfasta}
echo ">S_vulgaris_NA_AMY2A_copy2" >> ${outfasta}
cat SvNA_AMY2A_copy2.fasta >> ${outfasta}
echo >> ${outfasta}
echo ">S_vulgaris_AU_AMY2A_copy1" >> ${outfasta}
cat SvAU_AMY2A_copy1.fasta >> ${outfasta}
echo >> ${outfasta}
echo ">S_vulgaris_AU_AMY2A_copy2" >> ${outfasta}
cat SvAU_AMY2A_copy2.fasta >> ${outfasta}
echo >> ${outfasta}
echo ">L_superbus_AMY2A_copy1" >> ${outfasta}
cat Lsuperbus_AMY2A_copy1.fasta >> ${outfasta}
echo >> ${outfasta}
echo ">L_superbus_AMY2A_copy2" >> ${outfasta}
cat Lsuperbus_AMY2A_copy2.fasta >> ${outfasta}
```

## MAFFT
```
ml MAFFT/7.505-gimkl-2022a-with-extensions
mkdir MAFFT
mafft --preservecase --maxiterate 1000 --localpair AMY2A_sequences.fasta > MAFFT/Aligned_AMY2A_sequences.fasta
```
## trimAl
Trim aligned sequences accordingly
```
ml trimAl/1.4.1-GCC-11.3.0 
mkdir trimal
trimal -in MAFFT/Aligned_AMY2A_sequences.fasta -out trimal/Aligned_AMY2A_sequences.trimal_gt0_5.fasta -gt .5 # sequence present in at least 50% of individuals
trimal -in MAFFT/Aligned_AMY2A_sequences.fasta -out trimal/Aligned_AMY2A_sequences.trimal_gt0.fasta -gt 0 # same as usual
trimal -in MAFFT/Aligned_AMY2A_sequences.fasta -out trimal/Aligned_AMY2A_sequences.trimal_gt1.fasta -gt 1 # sequence present in all individuals
```
After some visualisation of alignments in [msaviewer](https://www.ncbi.nlm.nih.gov/projects/msaviewer/), there is some speculation of the first 11 bases of the Aligned_AMY2A_sequences.trimal_gt1.fasta file and therefore the first 11 bases are trimmed as follows:

```
ml trimAl/1.4.1-GCC-11.3.0 
trimal -in trimal/Aligned_AMY2A_sequences.trimal_gt1.fasta(-selectcols { 0-10 }) -out trimal/Aligned_AMY2A_sequences.trimal_gt1_trimstart.fasta
```
## IQ-TREE
```
cd trimal
ml IQ-TREE/2.2.2.2-gimpi-2022a
iqtree -s Aligned_AMY2A_sequences.trimal_gt1.fasta
iqtree -s Aligned_AMY2A_sequences.trimal_gt0.fasta
iqtree -s Aligned_AMY2A_sequences.trimal_gt0_5.fasta
iqtree -s Aligned_AMY2A_sequences.trimal_gt1_trimstart.fasta
```
Resulting `.treefile` were visualised in [itol](https://itol.embl.de/upload.cgi)
