#!/bin/bash
#SBATCH -n 8                # Number of cores
#SBATCH -N 1
#SBATCH -t 7-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared,serial_requeue   # Partition to submit to
#SBATCH --mem=40G           # Memory pool for all cores (see also --mem-per-cpu)


source ~/mambaforge/bin/activate stringtie

stringtie -p 8 --merge gtf_list.txt -o ./scrubjay_stringtie_merge.gtf
