#!/bin/bash
#SBATCH --partition=popgenom	### Partition
#SBATCH --job-name=CntMatr      ### Job Name
#SBATCH --output=CntMatr.out	### File in which to store job output
#SBATCH --error=CntMatr.err		### File in which to store job error messages
#SBATCH --time=7-00:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlcoffin@ksu.edu


python prepDE.py -g /data/kelley/projects/tarcreek/data/gene_expression/Gaffinis_Gene_Matrix.csv -t /data/kelley/projects/tarcreek/data/gene_expression/Gaffinis_Transcript_Matrix.csv

#done