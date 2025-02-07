#!/bin/bash
#SBATCH -J read_stats # A single job name for the array
#SBATCH --partition=shared # Partition
#SBATCH -n 4 # Number of Nodes required
#SBATCH --mem=24G # Memory request 
#SBATCH -t 4-00:00 # Maximum execution time (D-HH:MM)

[ -e COMMUNITIES_STATS.tsv ] && rm COMMUNITIES_STATS.tsv

while read sample; do

    python get_communities_stats.py $sample >> COMMUNITIES_STATS.tsv

done < community_fai_indexes.txt
