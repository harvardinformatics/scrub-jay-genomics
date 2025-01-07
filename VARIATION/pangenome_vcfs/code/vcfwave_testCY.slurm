#!/bin/bash
#SBATCH -J vcfwave # A single job name for the array
#SBATCH -a 1-46
#SBATCH --partition=shared # Partition
#SBATCH -n 8 # Number of Nodes required
#SBATCH --mem=64G # Memory request 
#SBATCH -t 3-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o vcfwave_%A_%a.out # Standard output
#SBATCH -e vcfwave_%A_%a.err # Standard error

source ~/mambaforge/bin/activate vcflib-manual

gzip testCY/community21_testCY.vcf

singularity exec /n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/pggb_native/pggb_cannon.sif vcfbub -l 0 -a 100000 --input testCY/community21_testCY.vcf.gz > testCY/community21_testCY.bub.vcf
/n/holylfs05/LABS/informatics/Users/dkhost/vcflib/build/vcfwave -I 1000 -t 8 testCY/community21_testCY.bub.vcf > testCY/community21_testCY.bub.wave.vcf
