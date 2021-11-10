#!/bin/bash
#SBATCH --partition=kamiak	### Partition
#SBATCH --job-name=alignXmac      ### Job Name
#SBATCH --output=alignXmac.out       ### File in which to store job output
#SBATCH --error=alignXmac.err        ### File in which to store job error messages
#SBATCH --time=7-00:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --array=1-53:1			### Submit separate jobs for each sample in the array
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --ntasks-per-node=2     ### Number of tasks to be launched per Node
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlcoffin@ksu.edu

###Be sure to change your reference genome location, we don't all want to map to Pmex

##"bwa mem" aligns reads from the files specified in the "qual24samplepaths.txt" to a reference genome.
#-t Number of threads to use on the cluster
#-R Change read group header line such as '@RG\tID:foo\tSM:bar' 
#-o specifies where to output the mapped reads (in sam format)

SampleName=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' qual24paths.txt | awk '{print $1}')
read1=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' qual24paths.txt | awk '{print $2}')
read2=$(sed -n ''$SLURM_ARRAY_TASK_ID'p' qual24paths.txt | awk '{print $3}')
		echo $SampleName
		output=/data/kelley/projects/tarcreek/data/mapping/X_maculatus_mapping/01_mapped_sam_files/$SampleName
		echo $output
       
       bwa mem -t 4 -R "@RG\tID:$SampleName\tSM:$SampleName\tLB:$SampleName\tPL:HiSeq2500" /data/kelley/projects/tarcreek/data/reference/Xmaculatus_reference/GCF_002775205.1_X_maculatus-5.0-male_genomic.fna $read1 $read2 -o $output
