#!/bin/bash
#SBATCH -n 12                # Number of cores
#SBATCH -N 1
#SBATCH -t 7-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared   # Partition to submit to
#SBATCH --mem=64G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o repmask.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e repmask.err  # File to which STDERR will be written, %j inserts jobid


/n/holyscratch01/informatics/dkhost/holylfs_shortcut/dfam-tetools.sh --container /n/holyscratch01/informatics/dkhost/holylfs_shortcut/dfam-tetools-latest.sif -- RepeatMasker -xsmall -pa 9 -gff -no_is -lib /n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/HIC_REFERENCE_ASSEMS/AW_365336_combined_repeats_v2.fasta ../AW_365336_FSJragtag.v1.fasta
