#!/bin/bash
#SBATCH -J hic-map # A single job name for the array
#SBATCH -a 1-3
#SBATCH --partition=holy-smokes # Partition
#SBATCH -n 16 # Number of Nodes required
#SBATCH -N 1
#SBATCH --mem=96G # Memory request 
#SBATCH -t 14-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o hic-map_%A_%a.out # Standard output
#SBATCH -e hic-map_%A_%a.err # Standard error

module load samtools
module load bwa/0.7.17-fasrc01
module load bedtools2/2.26.0-fasrc01

export SLURM_TMPDIR=/n/holyscratch01/informatics/dkhost/tmp
export TMPDIR=/n/holyscratch01/informatics/dkhost/tmp

#Assign assembly and read files from spec sheet
INFO=$(sed "${SLURM_ARRAY_TASK_ID}q;d" hic_fofn.txt)

SAMPLE=$(echo $INFO | awk 'BEGIN{FS=" "} {print $1}')

ASSEM="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/assemblies/${SAMPLE}.p_ctg.fa"

READ1=$(echo $INFO | awk 'BEGIN{FS=" "} {split($4,a,",")} {print a[1]}')
READ2=$(echo $INFO | awk 'BEGIN{FS=" "} {split($4,a,",")} {print a[2]}')

MAP_DIR="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/HiC_mapping/"


#Paths to scripts from Arima and Esrice pipelines
FILTER="/n/holyscratch01/informatics/dkhost/holylfs_shortcut/mapping_pipeline/filter_five_end.pl"
COMBINE="/n/holyscratch01/informatics/dkhost/holylfs_shortcut/esrice-hic-pipeline/combine_ends.py"


#Index draft genome
#bwa index $ASSEM


#Map reads
bwa mem -t 16 $ASSEM $READ1 | samtools view -F 4 -Sb - > ${MAP_DIR}/${SAMPLE}_R1.bam

bwa mem -t 16 $ASSEM $READ2 | samtools view -F 4 -Sb - > ${MAP_DIR}/${SAMPLE}_R2.bam


#Filter
samtools view -h ${MAP_DIR}/${SAMPLE}_R1.bam | perl $FILTER | samtools view -Sb - > ${MAP_DIR}/${SAMPLE}_filt_R1.bam

samtools view -h ${MAP_DIR}/${SAMPLE}_R2.bam | perl $FILTER | samtools view -Sb - > ${MAP_DIR}/${SAMPLE}_filt_R2.bam


#Combine and de-dup
samtools faidx $ASSEM

python $COMBINE ${MAP_DIR}/${SAMPLE}_filt_R1.bam ${MAP_DIR}/${SAMPLE}_filt_R2.bam | samtools fixmate -m - - | samtools sort -@ 16 -T /n/holyscratch01/informatics/dkhost/tmp - | samtools markdup -@ 16 -r - ${MAP_DIR}/${SAMPLE}_combined.bam


#Convert to BED
bamToBed -i ${MAP_DIR}/${SAMPLE}_combined.bam | sort -k 4 > ${MAP_DIR}/${SAMPLE}_combined.bed
