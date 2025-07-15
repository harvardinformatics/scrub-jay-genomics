#!/bin/bash
#SBATCH -n 4                # Number of cores
#SBATCH -N 1
#SBATCH -t 3-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared,serial_requeue   # Partition to submit to
#SBATCH --mem=32G           # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH -a 1-24
#SBATCH -o kallisto.%J.%a.out # File to which STDOUT will be written
#SBATCH -e kallisto.%J.%a.err # File to which STDERR will be written



SAMPLE=$(sed "${SLURM_ARRAY_TASK_ID}q;d" sample_list.txt)

INDEX="/n/holyscratch01/jlewis_lab/pnair/humanTranscriptome.kalIdx"
FASTQDIR="/n/holyscratch01/jlewis_lab/pnair/RNAseq_dir"
READ1="${FASTQDIR}/${SAMPLE}_R1_001.fastq.gz"
READ2="${FASTQDIR}/${SAMPLE}_R2_001.fastq.gz"

#echo $SAMPLE

mkdir -p "/n/holyscratch01/jlewis_lab/pnair/${SAMPLE}_kallisto"

OUTDIR="/n/holyscratch01/jlewis_lab/pnair/${SAMPLE}_kallisto"

source ~/mambaforge/bin/activate kallisto

kallisto quant -t 4 -b 30 -i $INDEX -o $OUTDIR $READ1 $READ2
