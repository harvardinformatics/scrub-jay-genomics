#!/bin/bash
#SBATCH -n 1                # Number of cores
#SBATCH -t 3-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=16G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o count.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e count.err  # File to which STDERR will be written, %j inserts jobid

python counts_repeats.py -i repMaskGff_fofn.txt > scrubjay_repeat_counts.tsv
