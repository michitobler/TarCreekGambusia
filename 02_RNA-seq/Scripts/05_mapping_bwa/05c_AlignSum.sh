#!/bin/bash
#SBATCH --partition=kamiak	### Partition
#SBATCH --job-name=alignsum      ### Job Name
#SBATCH --output=alignsum.out       ### File in which to store job output
#SBATCH --error=alignsum.err        ### File in which to store job error messages
#SBATCH --time=7-00:00:00       ### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlcoffin@ksu.edu

for i in {1..53}

do

        reference=/data/kelley/projects/tarcreek/data/reference/Xmaculatus_reference/GCF_002775205.1_X_maculatus-5.0-male_genomic.fna
        Input=$(cat 03_bam_paths.txt | sed -n ''$i'p' | awk '{print $2}')
        SampleName=$(cat 03_bam_paths.txt | sed -n ''$i'p' | awk '{print $1}')
        echo $SampleName
       output=/data/kelley/projects/tarcreek/data/mapping/X_maculatus_mapping/03_alignment_summary/$SampleName.txt
       echo $output


java -jar /data/kelley/projects/programs/picard-tools-1.141/picard.jar CollectAlignmentSummaryMetrics \
        R=$reference \
        I=$Input \
        O=$output \
        
done