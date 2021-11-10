#!/bin/bash
#SBATCH --partition=popgenom	### Partition
#SBATCH --job-name=strtie1      ### Job Name
#SBATCH --output=strtie1.out       ### File in which to store job output
#SBATCH --error=strtie1.err        ### File in which to store job error messages
#SBATCH --array=1-53:1			### Submit separate jobs for each sample in the array
#SBATCH --time=7-00:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --ntasks-per-node=2     ### Number of tasks to be launched per Node

##Run StringTie to to extract putatively expressed regions for each individual (reference guided with the .gff file)
##Version 1.3.2d
## -e only estimate the abundance of given reference transcripts (requires -G)
## -B enable output of Ballgown table files which will be created in the same directory as the output GTF (requires -G, -o recommended)
## -G reference annotation to use for guiding the assembly process (GTF/GFF3)


input=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' 01_mapped_bam_paths.txt | awk '{print $2}')
SampleName=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' 01_mapped_bam_paths.txt | awk '{print $1}')

	output=/data/kelley/projects/tarcreek/data/gene_expression/ballgown/$SampleName/$SampleName.gtf
	echo $SampleName

mkdir /data/kelley/projects/tarcreek/data/gene_expression/ballgown/$SampleName
cd /data/kelley/projects/tarcreek/data/gene_expression
stringtie -e -B -G /data/kelley/projects/tarcreek/data/reference/Xmaculatus_reference/GCF_002775205.1_X_maculatus-5.0-male_genomic.gff -o $output $input

	echo $SampleName

