#!/bin/bash -e
#SBATCH --job-name=02a
#SBATCH --output=02_errors/02a_%j.out
#SBATCH --error=02_errors/02a_%j.err
#SBATCH --mail-user=kats326@aucklanduni.ac.nz
#SBATCH --mail-type=ALL
#SBATCH --time=1:00:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --account=uoa02613
#SBATCH --profile task
#SBATCH --cpus-per-task=4

module purge
module load Miniconda3/22.11.1-1

source activate /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/programs/miniconda/envs/mamba 
#source activate conda activate /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/programs/miniconda/envs/mamba 
#conda activate /nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/programs/miniconda/envs/mamba

wdir=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/Potential_repeats/
mapdir=/nesi/nobackup/uoa02613/A_selection_analyses/selection_analyses/WGS/data/processed/mapped_reads/
cd ${wdir}
mkdir -p repeat_ID_482137

chr=chr8
#outputfolder=repeat_ID_482137/Madhya_Pradesh
#samplot plot -o ${outputfolder} -c ${chr} -s 20678727 -e 20687524 -t DEL -d 100 -n J761 J762 J763 J764 J766 J767 J768 -b ${mapdir}J761.sorted.dup.reheader.bam ${mapdir}J762.sorted.dup.reheader.bam ${mapdir}J763.sorted.dup.reheader.bam ${mapdir}J764.sorted.dup.reheader.bam ${mapdir}J766.sorted.dup.reheader.bam ${mapdir}J767.sorted.dup.reheader.bam ${mapdir}J768.sorted.dup.reheader.bam

#outputfolder=repeat_ID_482137/Tamil_Nadu
#samplot plot -o ${outputfolder} -c ${chr} -s 20678727 -e 20687524 -t DEL -d 100 -n J753 J754 J755 J756 J757 J758 J759 J760 -b ${mapdir}J753.sorted.dup.reheader.bam ${mapdir}J754.sorted.dup.reheader.bam ${mapdir}J755.sorted.dup.reheader.bam ${mapdir}J756.sorted.dup.reheader.bam ${mapdir}J757.sorted.dup.reheader.bam ${mapdir}J758.sorted.dup.reheader.bam ${mapdir}J759.sorted.dup.reheader.bam ${mapdir}J760.sorted.dup.reheader.bam

#popid=(J769 J770 J771 J773 J774 J775 J776)
#outputfolder=repeat_ID_482137/Cairns  
#samplot plot -o ${outputfolder} -c ${chr} -s 20678727 -e 20687524 -t DEL -d 100 -n ${popid[0]} ${popid[1]} ${popid[2]} ${popid[3]} ${popid[4]} ${popid[5]} ${popid[6]} -b ${mapdir}${popid[0]}.sorted.dup.reheader.bam ${mapdir}${popid[1]}.sorted.dup.reheader.bam ${mapdir}${popid[2]}.sorted.dup.reheader.bam ${mapdir}${popid[3]}.sorted.dup.reheader.bam ${mapdir}${popid[4]}.sorted.dup.reheader.bam ${mapdir}${popid[5]}.sorted.dup.reheader.bam ${mapdir}${popid[6]}.sorted.dup.reheader.bam 

#popid=(J777 J778 J779 J780 J781 J782 J783 J784)
#outputfolder=repeat_ID_482137/Melbourne  
#samplot plot -o ${outputfolder} -c ${chr} -s 20678727 -e 20687524 -t DEL -d 100 -n ${popid[0]} ${popid[1]} ${popid[2]} ${popid[3]} ${popid[4]} ${popid[5]} ${popid[6]} ${popid[7]} -b ${mapdir}${popid[0]}.sorted.dup.reheader.bam  ${mapdir}${popid[1]}.sorted.dup.reheader.bam ${mapdir}${popid[2]}.sorted.dup.reheader.bam ${mapdir}${popid[3]}.sorted.dup.reheader.bam ${mapdir}${popid[4]}.sorted.dup.reheader.bam ${mapdir}${popid[5]}.sorted.dup.reheader.bam ${mapdir}${popid[6]}.sorted.dup.reheader.bam ${mapdir}${popid[7]}.sorted.dup.reheader.bam 

popid=(J785 J786 J787 J788 J789 J790 J791 J792)
outputfolder=repeat_ID_482137/South_Africa  
samplot plot -o ${outputfolder} -c ${chr} -s 20678727 -e 20687524 -t DEL -d 100 -n ${popid[0]} ${popid[1]} ${popid[2]} ${popid[3]} ${popid[4]} ${popid[5]} ${popid[6]} ${popid[7]} -b ${mapdir}${popid[0]}.sorted.dup.reheader.bam  ${mapdir}${popid[1]}.sorted.dup.reheader.bam ${mapdir}${popid[2]}.sorted.dup.reheader.bam ${mapdir}${popid[3]}.sorted.dup.reheader.bam ${mapdir}${popid[4]}.sorted.dup.reheader.bam ${mapdir}${popid[5]}.sorted.dup.reheader.bam ${mapdir}${popid[6]}.sorted.dup.reheader.bam ${mapdir}${popid[7]}.sorted.dup.reheader.bam

popid=(J793 J794 J795 J796 J797 J798 J799)
outputfolder=repeat_ID_482137/Fiji  
samplot plot -o ${outputfolder} -c ${chr} -s 20678727 -e 20687524 -t DEL -d 100 -n ${popid[0]} ${popid[1]} ${popid[2]} ${popid[3]} ${popid[4]} ${popid[5]} ${popid[6]} -b ${mapdir}${popid[0]}.sorted.dup.reheader.bam ${mapdir}${popid[1]}.sorted.dup.reheader.bam ${mapdir}${popid[2]}.sorted.dup.reheader.bam ${mapdir}${popid[3]}.sorted.dup.reheader.bam ${mapdir}${popid[4]}.sorted.dup.reheader.bam ${mapdir}${popid[5]}.sorted.dup.reheader.bam ${mapdir}${popid[6]}.sorted.dup.reheader.bam

popid=(J718 J719 J720 J721 J722 J723 J724 J725)
outputfolder=repeat_ID_482137/Leigh  
samplot plot -o ${outputfolder} -c ${chr} -s 20678727 -e 20687524 -t DEL -d 100 -n ${popid[0]} ${popid[1]} ${popid[2]} ${popid[3]} ${popid[4]} ${popid[5]} ${popid[6]} ${popid[7]} -b ${mapdir}${popid[0]}.sorted.dup.reheader.bam  ${mapdir}${popid[1]}.sorted.dup.reheader.bam ${mapdir}${popid[2]}.sorted.dup.reheader.bam ${mapdir}${popid[3]}.sorted.dup.reheader.bam ${mapdir}${popid[4]}.sorted.dup.reheader.bam ${mapdir}${popid[5]}.sorted.dup.reheader.bam ${mapdir}${popid[6]}.sorted.dup.reheader.bam ${mapdir}${popid[7]}.sorted.dup.reheader.bam

popid=(J726 J727 J728 J729 J730 J731 J732 J733)
outputfolder=repeat_ID_482137/Napier  
samplot plot -o ${outputfolder} -c ${chr} -s 20678727 -e 20687524 -t DEL -d 100 -n ${popid[0]} ${popid[1]} ${popid[2]} ${popid[3]} ${popid[4]} ${popid[5]} ${popid[6]} ${popid[7]} -b ${mapdir}${popid[0]}.sorted.dup.reheader.bam  ${mapdir}${popid[1]}.sorted.dup.reheader.bam ${mapdir}${popid[2]}.sorted.dup.reheader.bam ${mapdir}${popid[3]}.sorted.dup.reheader.bam ${mapdir}${popid[4]}.sorted.dup.reheader.bam ${mapdir}${popid[5]}.sorted.dup.reheader.bam ${mapdir}${popid[6]}.sorted.dup.reheader.bam ${mapdir}${popid[7]}.sorted.dup.reheader.bam 

popid=(J734 J735 J736 J737 J738 J739 J740 J741)
outputfolder=repeat_ID_482137/Great_Barrier_Island  
samplot plot -o ${outputfolder} -c ${chr} -s 20678727 -e 20687524 -t DEL -d 100 -n ${popid[0]} ${popid[1]} ${popid[2]} ${popid[3]} ${popid[4]} ${popid[5]} ${popid[6]} ${popid[7]} -b ${mapdir}${popid[0]}.sorted.dup.reheader.bam  ${mapdir}${popid[1]}.sorted.dup.reheader.bam ${mapdir}${popid[2]}.sorted.dup.reheader.bam ${mapdir}${popid[3]}.sorted.dup.reheader.bam ${mapdir}${popid[4]}.sorted.dup.reheader.bam ${mapdir}${popid[5]}.sorted.dup.reheader.bam ${mapdir}${popid[6]}.sorted.dup.reheader.bam ${mapdir}${popid[7]}.sorted.dup.reheader.bam 

popid=(J742 J743 J744 J745 J746 J747)
outputfolder=repeat_ID_482137/Maharashtra_subpop_A  
samplot plot -o ${outputfolder} -c ${chr} -s 20678727 -e 20687524 -t DEL -d 100 -n ${popid[0]} ${popid[1]} ${popid[2]} ${popid[3]} ${popid[4]} ${popid[5]} -b ${mapdir}${popid[0]}.sorted.dup.reheader.bam ${mapdir}${popid[1]}.sorted.dup.reheader.bam ${mapdir}${popid[2]}.sorted.dup.reheader.bam ${mapdir}${popid[3]}.sorted.dup.reheader.bam ${mapdir}${popid[4]}.sorted.dup.reheader.bam ${mapdir}${popid[5]}.sorted.dup.reheader.bam 

popid=(J748 J749 J750 J751 J752)
outputfolder=repeat_ID_482137/Maharashtra  
samplot plot -o ${outputfolder} -c ${chr} -s 20678727 -e 20687524 -t DEL -d 100 -n ${popid[0]} ${popid[1]} ${popid[2]} ${popid[3]} ${popid[4]} -b ${mapdir}${popid[0]}.sorted.dup.reheader.bam ${mapdir}${popid[1]}.sorted.dup.reheader.bam ${mapdir}${popid[2]}.sorted.dup.reheader.bam ${mapdir}${popid[3]}.sorted.dup.reheader.bam ${mapdir}${popid[4]}.sorted.dup.reheader.bam 
