#!/bin/bash
#SBATCH --partition=kamiak	### Partition
#SBATCH --job-name=TrimQ24      	### Job Name
#SBATCH --output=TrimQ24.out       ### File in which to store job output
#SBATCH --error=TrimQ24.err        ### File in which to store job error messages
#SBATCH --time=7-00:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --array=1-53:1			### Submit separate jobs for each sample in the array
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlcoffin@ksu.edu

###Runs Trim Galore on the files specified in the "qual0samplepaths.txt".  Note: this is a tab delimited file.

SampleName=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' qual0samplepaths.txt | awk '{print $1}')
read1=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' qual0samplepaths.txt | awk '{print $2}')
read2=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' qual0samplepaths.txt | awk '{print $3}')
                
       trim_galore --path_to_cutadapt /data/kelley/projects/programs/cutadapt/bin/cutadapt --quality 24 --fastqc_args "--noextract --nogroup" --stringency 6 -e 0.1 --gzip --length 50 --output_dir /data/kelley/projects/tarcreek/data/trimmed/qual24 --paired $read1 $read2


