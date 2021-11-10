#!/bin/bash
#SBATCH --output=TrimQ0_last27.out	  ### File in which to store job output
#SBATCH --error=TrimQ0_last27.err        ### File in which to store job error messages
#SBATCH --time=7-00:00:00	### Wall clock time limit in Days-HH:MM:SS
#SBATCH --nodes=1               ### Node count required for the job
#SBATCH --ntasks-per-node=1     ### Number of tasks to be launched per Node
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlcoffin@ksu.edu

# Runs Trim Galore on the files specified in the "Limia_samplepaths.txt".  Note: this is a tab delimited file.

# The for loop is hard coded to the number to lines in the "samplepaths.txt" file. Right now I have 12
#     entries in the "Limia_samplepaths.txt" table so the loop is {1..12}. If you only had 10 sequencing file pairs
#     you'd need to change the loop to {1..10}.

# "cat" opens the file.
# "sed" is grabbing the ith line of the file.
# "awk is grabbing the 1st, 2nd, or 3rd column from the line specified.
# "mikdir" is creating a directory with the name supplied in "samplepaths.txt"

# The end result of running this script is trimmed files that match and a fastqc report that
# matches the directory structure of your pre-processed read data.

##--quality trims reads with a quality score of 0, which shouldn't be many reads. So this round of trimming is mostly just removing adapters and really bad quality reads
##--stringency anything with 6bp overlap with the adapters will be trimmed
##-e ERROR RATE: Maximum allowed error rate (no. of errors divided by the length of the matchingregion) (default: 0.1)
##--gzip ouputs a gzipped file
##--length Discard reads that became shorter than specified length (50) because of either quality or adapter trimming
##--paired location of the paired reads that will be mapped together


##--fastqc_args
#--noextract outputs a compressed file
#--nogroup  Disables grouping of bases for reads >50bp


for i in {1..27}

do
        read1=$(cat affinis_raw_paths_last27.txt | sed -n ''$i'p' | awk '{print $2}')
        read2=$(cat affinis_raw_paths_last27.txt | sed -n ''$i'p' | awk '{print $3}')
        SampleName=$(cat affinis_raw_paths_last27.txt | sed -n ''$i'p' | awk '{print $1}')
        echo $SampleName

      trim_galore --path_to_cutadapt /data/kelley/projects/programs/cutadapt/bin/cutadapt --quality 0 --fastqc_args "--noextract --nogroup" --stringency 6 -e 0.1 --gzip --length 50 --output_dir /data/kelley/projects/tarcreek/data/trimmed/qual0 --paired $read1 $read2
done

