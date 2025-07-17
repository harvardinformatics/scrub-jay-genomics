#!/bin/bash
#SBATCH -n 12                # Number of cores
#SBATCH -N 1
#SBATCH -t 14-0:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p holy-cow   # Partition to submit to
#SBATCH --mem=120G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -o repmod.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e repmod.err  # File to which STDERR will be written, %j inserts jobid

/n/holylfs05/LABS/informatics/Users/dkhost/dfam-tetools.sh --container /n/holylfs05/LABS/informatics/Users/dkhost/dfam-tetools-latest.sif -- BuildDatabase -name AW_365336_FSJragtag.v1 ../AW_365336_FSJragtag.v1.fasta

/n/holylfs05/LABS/informatics/Users/dkhost/dfam-tetools.sh --container /n/holylfs05/LABS/informatics/Users/dkhost/dfam-tetools-latest.sif -- RepeatModeler -database AW_365336_FSJragtag.v1 -pa 12 -LTRStruct
