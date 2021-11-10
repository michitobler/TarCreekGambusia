#!/bin/bash
#SBATCH --partition=kamiak	
#SBATCH --job-name=concat_align_sum     
#SBATCH --output=concat_align_sum.out	
#SBATCH --error=concat_align_sum.err       
#SBATCH --time=7-00:00:00	
#SBATCH --nodes=1              
#SBATCH --ntasks-per-node=1    
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jlcoffin@ksu.edu

awk '{print FILENAME, $0}' 18105FL*.txt | cat > Mapping_Stats_ALL-INFO-3.txt
grep -e ".txt PAIR*" Mapping_Stats_ALL-INFO-3.txt > Mapping_Stats-3.txt

