#!/bin/bash
#SBATCH --output=count_out.txt
#SBATCH --error=count_error.txt
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlcoffin@ksu.edu

# Script locates all fastq files from current directory to all subfolders and goes through each fastq file and counts the number of reads in each file.  This script only counts the number of reads in the read 1 file because the number of reads in the read 2 file should be the same.
####################################################################################################################################

# Find all read 1 fastq files from current directory and all subfolders. Output full path (`pwd` needed otherwise you'd only get the path from the current directory). -type specifies a normal file. -name allows us to look for things based on file name.

fileList=($(find `pwd` -type f -name '*?????R1*.fastq.gz'))

# Initialize variables for use below. Change the output file name here.
Output_file_name="reads_per_fastq.csv"
Total_Number_Of_Reads=0
Array_Length=${#fileList[@]}
Loop_Number=0

# Message to user
echo "Counting all the things..."


# Header for the output CSV file
echo "File_Name,number_of_reads" > $Output_file_name

# Loops through all elements in the array. Note that "i" holds the value for the element in the array during each iteration. In this case "i" holds the path to the fastq file.
for i in ${fileList[@]}
do

	Loop_Number=`expr $Loop_Number + 1`

	File_Name=`echo $i | awk -F "/" '{print $NF}'`
	Number_of_Lines=`gunzip -c $i | wc -l`
	Number_Of_Reads=`expr $Number_of_Lines / 4`
	
	Total_Number_Of_Reads=`expr $Total_Number_Of_Reads + $Number_Of_Reads`

	echo "$File_Name,$Number_Of_Reads" >> $Output_file_name

	# print message about what file has been completed
	echo "Finished $Loop_Number/$Array_Length"

done

# Print whitespace line and total read to output file.
echo "" >> $Output_file_name
echo "Total_Number_Of_Reads,$Total_Number_Of_Reads" >> $Output_file_name

# Message to user
echo "Finished counting all the things!!!"