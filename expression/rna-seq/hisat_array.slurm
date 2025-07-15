#!/bin/bash
#SBATCH -n 8                # Number of cores
#SBATCH -N 1
#SBATCH -t 7-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared,serial_requeue   # Partition to submit to
#SBATCH --mem=40G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -a 1-62
#SBATCH -o hisat.%A.%a.out # File to which STDOUT will be written
#SBATCH -e hisat.%A.%a.err # File to which STDERR will be written


export SLURM_TMPDIR=/n/holyscratch01/informatics/dkhost/tmp
export TMPDIR=/n/holyscratch01/informatics/dkhost/tmp

SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" sample_list.txt)

INDEX="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/HIC_REFERENCE_ASSEMS/AW_365336_FSJragtag.v1.ALTLABEL"
FASTQDIR="/n/holylfs05/LABS/informatics/Everyone/scrubjay/RNAseq/fastq/Scrub_Jay_"
READ1="${FASTQDIR}/${SAMPLE}_R1_001.fastq.gz"
READ2="${FASTQDIR}/${SAMPLE}_R2_001.fastq.gz"
OUTSAM="${SAMPLE}_aligned.sam"
OUTBAM="${SAMPLE}_aligned.sorted.bam"

#echo $SAMPLE

mkdir -p "${SAMPLE}_hisat"

source ~/mambaforge/bin/activate hisat2

hisat2 -p 8 -x $INDEX -1 $READ1 -2 $READ2 -S ${SAMPLE}_hisat/$OUTSAM

source ~/mambaforge/bin/activate samtools

samtools sort -@ 8 -T $(mktemp -d) -O bam -o ${SAMPLE}_hisat/$OUTBAM ${SAMPLE}_hisat/$OUTSAM 
