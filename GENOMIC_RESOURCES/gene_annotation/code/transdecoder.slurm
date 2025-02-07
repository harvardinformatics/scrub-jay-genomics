#!/bin/bash
#SBATCH -n 6               # Number of cores
#SBATCH -N 1
#SBATCH -t 7-00:00          # Runtime in D-HH:MM, minimum of 10 minutes
#SBATCH -p shared,serial_requeue   # Partition to submit to
#SBATCH --mem=40G           # Memory pool for all cores (see also --mem-per-cpu)

GTF2FASTA="/n/holylfs05/LABS/informatics/Users/dkhost/TransDecoder-TransDecoder-v5.7.0/util/gtf_genome_to_cdna_fasta.pl"
GTF2GFF3="/n/holylfs05/LABS/informatics/Users/dkhost/TransDecoder-TransDecoder-v5.7.0/util/gtf_to_alignment_gff3.pl"
TDLONGORFS="/n/holylfs05/LABS/informatics/Users/dkhost/TransDecoder-TransDecoder-v5.7.0/TransDecoder.LongOrfs"
TDPREDICT="/n/holylfs05/LABS/informatics/Users/dkhost/TransDecoder-TransDecoder-v5.7.0/TransDecoder.Predict"
CDNA2GENOME="/n/holylfs05/LABS/informatics/Users/dkhost/TransDecoder-TransDecoder-v5.7.0/util/cdna_alignment_orf_to_genome_orf.pl"

GTF="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/gene_annotation/stringtie_RNAseq/scrubjay_stringtie_merge.gtf"
REF="/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/HIC_REFERENCE_ASSEMS/AW_365336_FSJragtag.v1.ALTLABEL.fasta"

$GTF2FASTA $GTF $REF > scrubjay_stringtie_transcripts.fa

$GTF2GFF3 $GTF > scrubjay_stringtie_merged.gff3

$TDLONGORFS -t scrubjay_stringtie_transcripts.fa
$TDPREDICT -t scrubjay_stringtie_transcripts.fa

$CDNA2GENOME scrubjay_stringtie_transcripts.fa.transdecoder.gff3 scrubjay_stringtie_merged.gff3 scrubjay_stringtie_transcripts.fa > scrubjay_stringtie_transcripts.fa.transdecoder.genome.gff3

#scrubjay_stringtie_transcripts.fa.transdecoder_dir/longest_orfs.pep

#FastaSplitter.py -f scrubjay_stringtie_transcripts.fa.transdecoder_dir/longest_orfs.pep -maxn 1000 -counter blastp/orf_splitfile_list.txt
