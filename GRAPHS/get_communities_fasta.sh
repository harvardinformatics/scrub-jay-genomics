#!/bin/bash
#SBATCH -J getfasta # A single job name for the array
#SBATCH --partition=holy-smokes # Partition
#SBATCH -n 4 # Number of Nodes required
#SBATCH --mem=32G # Memory request 
#SBATCH -t 2-0:00 # Maximum execution time (D-HH:MM)

source ~/anaconda3/bin/activate samtools

seq 1 946 | while read i; do 
#    echo "$i"
    variable=$(sed -n "${i}"p no_ref_community_list.txt | awk 'BEGIN{FS=" "} {print $2}')
#    echo "$variable"
    samtools faidx combined_assemblies.fasta $(cat communities/combined_assemblies.partition.paf.edges.weights.txt.community.${variable}.txt)  > no_ref_communities/allbird_community.${variable}.fa
    samtools faidx no_ref_communities/allbird_community.${variable}.fa
done
