#!/bin/bash
#SBATCH -J read_stats # A single job name for the array
#SBATCH --partition=shared # Partition
#SBATCH -n 4 # Number of Nodes required
#SBATCH --mem=24G # Memory request 
#SBATCH -t 4-00:00 # Maximum execution time (D-HH:MM)

[ -e results/assembly_qc/READS_STATS.tsv ] && rm results/assembly_qc/READS_STATS.tsv

while read sample; do

    ID=$(echo "$sample" | cut -f1)
    READS=$(echo "$sample" | cut -f3 | sed 's/,/ /g')

    cat $READS | seqkit stats -T -i $ID >> results/assembly_qc/READS_STATS.tsv

#    printf "%s\t%s\n" "$ID" "$READS"

done < sample_sheet.txt
