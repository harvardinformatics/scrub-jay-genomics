#!/bin/bash
#SBATCH -J isoquant # A single job name for the array
#SBATCH --partition=shared # Partition
#SBATCH -n 8 # Number of Nodes required
#SBATCH -N 1
#SBATCH --mem=72G # Memory request 
#SBATCH -t 7-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o isoquant.out # Standard output
#SBATCH -e isoquant.err # Standard error

source ~/anaconda3/bin/activate isoquant 

REF="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/HIC_REFERENCE_ASSEMS/AW_365336_FSJragtag.v1.ALTLABEL.fasta"
FASTQ1="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/data/hifi_annotation/m64190e_220824_192345.hifi_reads.fastq.gz"
FASTQ2="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/data/hifi_annotation/m64190e_220909_182409.hifi_reads.fastq.gz"

#isoquant.py --threads 12 --reference $REF --fastq $FASTQ1 $FASTQ2 --data_type nanopore -o /n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/GenomeAnnotation/IsoQuant/

source ~/anaconda3/bin/activate agat
agat_sp_merge_annotations.pl -f 00_m64190e_220824_192345.hifi_reads/00_m64190e_220824_192345.hifi_reads.fastq.transcript_models.gtf -f 01_m64190e_220909_182409.hifi_reads/01_m64190e_220909_182409.hifi_reads.fastq.transcript_models.gtf --out isoquant_merge.gtf
