#!/bin/bash
#SBATCH --partition=popgenom
#SBATCH --cpus-per-task=1
#SBATCH --array=1-24:1 #only 24 of the samples had more than 20M paired reads, the others don't need to be subsampled
#SBATCH --output=seqtk.out
#SBATCH --error=seqtk.err
#SBATCH --mem=50G

inputR1=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' qual24paths_to_be_subsampled.txt | awk '{print $2}')
inputR2=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' qual24paths_to_be_subsampled.txt | awk '{print $3}')
SampleName=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' qual24paths_to_be_subsampled.txt | awk '{print $1}')


seqtk sample -s 100 $inputR1 13000000 > "$SampleName".R1.sub.fq
seqtk sample -s 100 $inputR2 13000000 > "$SampleName".R2.sub.fq