#!/bin/bash
#SBATCH -J yahs_array # A single job name for the array
#SBATCH -a 1-2
#SBATCH --partition=holy-smokes # Partition
#SBATCH -n 12 # Number of Nodes required
#SBATCH -N 1
#SBATCH --mem=92G # Memory request 
#SBATCH -t 7-0:00 # Maximum execution time (D-HH:MM)
#SBATCH -o yahs_%A_%a.out # Standard output
#SBATCH -e yahs_%A_%a.err # Standard error

module load Java/1.8
module load samtools/1.5-fasrc02

YAHS="/n/holyscratch01/informatics/dkhost/holylfs_shortcut/yahs/yahs"
JUICER_PRE="/n/holyscratch01/informatics/dkhost/holylfs_shortcut/yahs/juicer_pre"
JUICER_JAR="/n/holyscratch01/informatics/dkhost/holylfs_shortcut/juicer/juicer_tools_1.22.01.jar"

#Assign assembly and read files from spec sheet
INFO=$(sed "${SLURM_ARRAY_TASK_ID}q;d" hic_fofn.txt)

SAMPLE=$(echo $INFO | awk 'BEGIN{FS=" "} {print $1}')

ASSEM="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/assemblies/${SAMPLE}.p_ctg.fa"
INDEX="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/assemblies/${SAMPLE}.p_ctg.fa.fai"
BED="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/HiC_mapping/${SAMPLE}_combined.bed"

#mkdir /n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/HiC_scaffolding_YaHS/${SAMPLE}
OUTDIR="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/HiC_scaffolding_YaHS/${SAMPLE}"

#Commands
#$YAHS -o $OUTDIR/${SAMPLE}_yahs_scaff $ASSEM $BED

cd $OUTDIR
#samtools faidx ${SAMPLE}_yahs_scaff_scaffolds_final.fa
SCAFFINDEX="${SAMPLE}_yahs_scaff_scaffolds_final.fa.fai"

#cut -f1,2 $SCAFFINDEX > scaffolds_final.chrom.sizes

#($JUICER_PRE ${SAMPLE}_yahs_scaff.bin ${SAMPLE}_yahs_scaff_scaffolds_final.agp $INDEX | sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G | awk 'NF' > alignments_sorted.txt.part) && (mv alignments_sorted.txt.part alignments_sorted.txt)

(java -Xms60g -Xmx90g -jar ${JUICER_JAR} pre --threads 12 alignments_sorted.txt out.hic.part scaffolds_final.chrom.sizes) && (mv out.hic.part out.hic) && (rm alignments_sorted.txt)
