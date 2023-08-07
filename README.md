# SNIBCPS

RNA-seq processing: $sbatch deg.sh
---------------------
As we are working in ComputeCanada's nodes, jobs have to be slurmed. The only thing to modify is the mail, put yours instead.
 (all files have to be decompressed once in your directory except fastq.gz files that will be automatically decompressed once in their directory)
 
1- Download mm10.fa.gz and mm10.refGene.gtf from UCSC browser (reference genome to be indexed, and gene annotation file respectively).
2- Download your fastq files.
3- Create a folder called indexGenRef
   command: mkdir indexGenRef
4- Create N number of directories depending on the number of the samples you have and name them as follow:  poolA, poolB, ......,poolN
   command: mkdir poolA poolB ......poolN
5- In our case we had pair-ended files, rename and move each 2 files per sample respectively to poolA, poolB, poolC ........ poolN
   command:  
   mv NS1111111111_R1.fastq poolA_R1.fastq.gz (renaming)
   mv NS1111111111_R2.fastq poolA_R2.fastq.gz (renaming)
   mv poolA_R1.fastq ./poolA (moving)
   mv poolA_R2.fastq ./poolA (moving)

As you see below, we are aligning 6 samples, if you have 4 or 8 or N number, you will have to align each samples, using the right directory, with the right files, see below:

STAR --genomeDir ./indexGenRef/ --runThreadN 8 --readFilesIn ./poolA/poolA_R1.fastq ./poolA/poolA_R2.fastq --outFileNamePrefix poolA --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100

STAR --genomeDir ./indexGenRef/ --runThreadN 8 --readFilesIn ./poolB/poolB_R1.fastq ./poolB/poolB_R2.fastq --outFileNamePrefix poolB --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100

STAR --genomeDir ./indexGenRef/ --runThreadN 8 --readFilesIn ./poolC/poolC_R1.fastq ./poolC/poolC_R2.fastq --outFileNamePrefix poolC --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100

STAR --genomeDir ./indexGenRef/ --runThreadN 8 --readFilesIn ./poolD/poolD_R1.fastq ./poolD/poolD_R2.fastq --outFileNamePrefix poolD --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100

STAR --genomeDir ./indexGenRef/ --runThreadN 8 --readFilesIn ./poolE/poolE_R1.fastq ./poolE/poolE_R2.fastq --outFileNamePrefix poolE --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100

STAR --genomeDir ./indexGenRef/ --runThreadN 8 --readFilesIn ./poolF/poolF_R1.fastq ./poolF/poolF_R2.fastq --outFileNamePrefix poolF --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFilterMultimapNmax 100 --winAnchorMultimapNmax 100

-- The lines will produce 6 sorted bam files:

poolBAligned.sortedByCoord.out.bam 
poolDAligned.sortedByCoord.out.bam 
poolFAligned.sortedByCoord.out.bam 
poolAAligned.sortedByCoord.out.bam
poolCAligned.sortedByCoord.out.bam
poolEAligned.sortedByCoord.out.bam

First 3 are mutant/treatment
Last 3 are control/wildtype

# TEtranscripts and DESeq2

cd ..
source ENV/bin/activate
cd ENV

Rscript des.R

TEtranscripts --sortByPos --mode multi -t poolBAligned.sortedByCoord.out.bam poolDAligned.sortedByCoord.out.bam poolFAligned.sortedByCoord.out.bam -c poolAAligned.sortedByCoord.out.bam poolCAligned.sortedByCoord.out.bam poolEAligned.sortedByCoord.out.bam --GTF mm10.refGene.gtf --TE mm10_rmsk_TE.gtf
