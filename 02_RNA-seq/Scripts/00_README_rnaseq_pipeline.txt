#Date Created: Feb 27, 2020
#Author: John L Coffin
#Collaborating Authors: Joanna L Kelley, Punidan D Jeyasingh, & Michael Tobler

For this experiment, we were interested in observing gene expression differences that were 
commonly shared between populations of Gambusia affinis in a polluted section of the Tar Creek 
Superfund site in NE Oklahoma and nearby uncontaminated streams. We captured 6 female 
Gambusia affinis at each of 3 sites (Tar Creek (polluted), Coal Creek (unpolluted), and 
Little Elm (unpolluted)), and excised gill, liver, and brain tissue from each individual 
for RNA-sequencing. Raw reads were downloaded on the Kamiak high performance computing cluster 
at Washington State U at /data/kelley/projects/tarcreek/data/raw/gz_files.

Step 1. 01_Reads_And_Index_Per_Fastq.sh
The first step in the analysis of these reads was to count the number of raw reads in each 
file to see if there were any issues that need correcting. 


Step 2. 02a_trimgalore_qual0.sh and then 02b_trimgalore_qual24.sh
Next, put the raw reads through two rounds of trimming based on the quality scores of each 
read in the fastq file. First round is quality score of 0, then use quality score of 24.


Step 3. 03_subsampling_excessively_large_files_seqtk.sh
After trimming, you can re-run a version of the 01_Reads_And_Index_Per_Fastq.sh script to 
count the number of reads in each file. After doing this, we realized that some samples had 
way more reads than the other samples. So to normalize this, we subsampled those samples 
with excessive reads with seqtk.


Steps 4 and 5. all scripts in the 04_mapping_bwa_subsamples folder and 05_mapping_bwa folder
We mapped all of the reads with BWA. In the 04 folder, we did this on the subsampled files, 
and then in the 05 folder we ran it on the rest of the samples. The first script 
(04a or 05a) maps the high-quality (qual 24+) reads to the Xiphophorus maculatus reference 
genome. The second script (04b or 05b) contains a pipeline to prep and convert the aligned 
.sam files into .bam files. It first removes soft-clipped reads, ensures that the read-pairs 
contain the correct information about their matching mate read, then converts sam to bam 
format and gives the new .bam file a header. Finally, the .bam file gets sorted in 
coordinate order.


Step 6. 06_stringtie.sh
To quantify and compare which transcripts are expressed, we used Stringtie, which extracts 
putatively expressed regions for each sample. 


Step 7. 07_GenerateCountsMatrices.sh and prepDE.py
In the 07_GenerateCountsMatrices.sh script, we used the prepDE.py script to generate a 
count matrix of each gene identified by Stringtie.


Step 8. 08_TarCreekGambusia_RNAseq_analyses.R
This script was then used to compare the gene count matrix for each tissue between 
individuals from Tar Creek compared to those from both Coal and Little Elm Creeks.


Step 9. 09_TarCreekGambusia_DAPC_visualize_expression_variation.R
This script enables visualization of gene expression variation among tissues and between sites 
using a technique called discriminant analysis of principal components (DAPC).

