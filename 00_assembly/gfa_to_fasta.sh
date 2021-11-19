#!/bin/bash
#SBATCH --partition=holy-smokes # Partition
#SBATCH -n 4 # Number of Nodes required
#SBATCH --mem=24G # Memory request 
#SBATCH -t 4-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o gfa_to_fasta.out # Standard output
#SBATCH -e gfa_to_fasta.err # Standard error

ls /n/holylfs05/LABS/informatics/Everyone/scrubjay/01_assembly/hifiasm_output/*.p_ctg.gfa > gfa_fofn.txt

while read p; do 

    sample=$(basename $p .gfa)
    awk '/^S/{print ">"$2;print $3}' $p | gzip > assemblies/${sample}.fa.gz

done < gfa_fofn.txt
