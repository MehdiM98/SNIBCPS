#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=0-12:00
#SBATCH -c 6
#SBATCH --mem-per-cpu=8000
###SBATCH --mem 16G
#SBATCH --account=def-pcampeau
#SBATCH -o Star.out
#SBATCH -e Star.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mehdi.manjoura@umontreal.ca


#Loading modules needed for the script to be run
 
module load gcc/9.3.0 star/2.7.9a
module load r/4.2.2
module load python/3.7.7

# Generating the reference genome index


STAR --runThreadN 8 --runMode genomeGenerate --genomeDir ./indexGenRef --genomeFastaFiles mm10.fa --sjdbGTFfile mm10.refGene.gtf --sjdbOverhang 99



#Aligning reads (.fastq) to reference genome (.fa)


#Before, please make sure you create N directories depending on the number of samples, in our case poolA poolB,..,poolF

#TO create directories: mkdir poolA poolB poolC poolD poolE poolF

#Move your fastq files respectively to directories (2 files per directory, forward and reverse), to move: mv file1_R1.fastq.gz ./nameofDirectory1 mv file1_R2.fastq.gz ./nameofDirectory1

 
cd ./poolA
gunzip *.gz
cd ..
STAR --genomeDir ./indexGenRef/ --runThreadN 8 --readFilesIn ./poolA/poolA_R1.fastq ./poolA/poolA_R2.fastq --outFileNamePrefix poolA --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100
cd ./poolB
gunzip *.gz
cd ..
STAR --genomeDir ./indexGenRef/ --runThreadN 8 --readFilesIn ./poolB/poolB_R1.fastq ./poolB/poolB_R2.fastq --outFileNamePrefix poolB --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100
cd ./poolC
gunzip *.gz
cd ..
STAR --genomeDir ./indexGenRef/ --runThreadN 8 --readFilesIn ./poolC/poolC_R1.fastq ./poolC/poolC_R2.fastq --outFileNamePrefix poolC --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100
cd ./poolD
gunzip *.gz
cd ..
STAR --genomeDir ./indexGenRef/ --runThreadN 8 --readFilesIn ./poolD/poolD_R1.fastq ./poolD/poolD_R2.fastq --outFileNamePrefix poolD --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100
cd ./poolE
gunzip *.gz
cd ..
STAR --genomeDir ./indexGenRef/ --runThreadN 8 --readFilesIn ./poolE/poolE_R1.fastq ./poolE/poolE_R2.fastq --outFileNamePrefix poolE --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100
cd ./poolF
gunzip *.gz
cd ..
STAR --genomeDir ./indexGenRef/ --runThreadN 8 --readFilesIn ./poolF/poolF_R1.fastq ./poolF/poolF_R2.fastq --outFileNamePrefix poolF --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100

# TEtranscripts and DESeq2

cd ..
source ENV/bin/activate
cd ENV

Rscript des.R

TEtranscripts --sortByPos --mode multi -t poolBAligned.sortedByCoord.out.bam poolDAligned.sortedByCoord.out.bam poolFAligned.sortedByCoord.out.bam -c poolAAligned.sortedByCoord.out.bam poolCAligned.sortedByCoord.out.bam poolEAligned.sortedByCoord.out.bam --GTF mm10.refGene.gtf --TE mm10_rmsk_TE.gtf



