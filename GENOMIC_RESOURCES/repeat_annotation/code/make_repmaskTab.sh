#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -t 1-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=16G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o fix.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e fix.err  # File to which STDERR will be written, %j inserts jobid

while read p; do
	dir=$(dirname $p)
	sample=$(basename $p)

	tail -n +4 ${dir}/${sample} | sed 's/^[ \t]*//' | sed 's/ \+/\t/g' > ${dir}/${sample}.tab

done < repmask_out_fofn.txt
