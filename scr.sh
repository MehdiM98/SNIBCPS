#!/bin/bash
#SBATCH --time=0-00:30:00                                  # Runtime is sam output/hour, see below, format --time=day-hours:minutes
###SBATCH --nodes=1                                     # In case all the tasks need to be on the same node
#SBATCH --ntasks=4                                      # Increase per sample would require other bowtie settings!
#SBATCH --mem-per-cpu 8000                              # Samtools uses ?, GB comes with 4 or 8GB per CPU
#SBATCH --account=def-pcampeau                  # Not your own account
#SBATCH -o GEAanalysisJune14H3.3_FreshHemic_E15vsIgG1_E15.out                  # See https://genome.ucsc.edu/goldenPath/help/bam.html for sam>>
#SBATCH -e GEAanalysisJune14H3.3_FreshHemic_E15vsIgG1.err                  # Give names to output prefix (here: toolsjobFeb4 and toolsFeb4>
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mehdi.manjoura@umontreal.ca



module load gcc/9.3.0
module load python/3.7.7
module load r/4.2.2
module load bedtools/2.30.0

mkdir Module1Files Module1Graphs Module2Files Module2Graphs Module3Files Module3Graphs Module4Files Module4Graphs Module5Files Module5Graphs Module6Files Module6Graphs

bash SEACR_1.3.sh ./BedGrphs/CHD3_2_E15.BEDgraph ./BedGrphs/IgG2_E15.BEDgraph non stringent output

mv output.stringent.bed sampleX.csv


Rscript ChIP.R
python updatedCode.py



#scp mehdim98@narval.computecanada.ca:scratch/codes/pipe/Module1Graphs/* /home/mehdim/Desktop/Directoryname/
