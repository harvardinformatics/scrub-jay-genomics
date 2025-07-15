#!/usr/bin/env bash
# ==============================================================================
#  Scrub‑Jay Dipcall + PSMC Pipeline
# ------------------------------------------------------------------------------
#  Purpose
#  -------
#  * Call diploid variants with Dipcall against the AW reference genome
#  * Produce phased VCFs and summary statistics per sample
#  * Convert Dipcall output to PSMC input and infer demographic history
#  * Generate PSMC plots (effective population size through time)
# ==============================================================================

# -----------------------------------------------------------------------------
# 0) Environment setup
# -----------------------------------------------------------------------------
cd /n/holyscratch01/edwards_lab/bfang/Scrub_jay/Dipcall

# -- Add Dipcall toolkit -------------------------------------------------------
export PATH=$PATH:/n/home00/bfang/programs/dipcall/dipcall.kit

# -- Reference & assemblies ----------------------------------------------------
PATH_asm=/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/assemblies
Ref_scrubjay=/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/HIC_REFERENCE_ASSEMS/AW_365336_FSJragtag.v1.fasta

# -- Required software modules (FASRC) ----------------------------------------
module load seqtk/1.2-fasrc01
module load gnuplot/5.2.6-fasrc01
module load ghostscript/9.16-fasrc01

module load python
source activate fasrc   # activates conda env containing vcf2snp.pl, etc.

# -- Add PSMC binaries ---------------------------------------------------------
export PATH=$PATH:/n/home00/bfang/programs/PSMC/psmc/
export PATH=$PATH:/n/home00/bfang/programs/PSMC/psmc/utils

# -----------------------------------------------------------------------------
# 1) Sample‑specific pipeline (example: "AA_SRR23446543")
#    ‑ Duplicate this block for each sample, adjusting IDs/paths accordingly.
# -----------------------------------------------------------------------------

###############################################
#### for example | sample "AA_SRR23446543" ####
###############################################

###############
### Dipcall ###
###############
/n/home00/bfang/programs/dipcall/dipcall.kit/run-dipcall AA_SRR23446543 \
  $Ref_scrubjay \
  $PATH_asm/AA_SRR23446543.hap1.p_ctg.fa \
  $PATH_asm/AA_SRR23446543.hap2.p_ctg.fa \
  > /n/holyscratch01/edwards_lab/bfang/Scrub_jay/Dipcall_results/AA_SRR23446543/AA_SRR23446543.mak

cd /n/holyscratch01/edwards_lab/bfang/Scrub_jay/Dipcall_results/AA_SRR23446543/
make -j2 -f AA_SRR23446543.mak 2>&1 | tee AA_SRR23446543.log

## -- Generate phased VCF ------------------------------------------------------
/n/home00/bfang/programs/dipcall/dipcall.kit/k8 \
  /n/home00/bfang/programs/dipcall/dipcall.kit/dipcall-aux.js vcfpair -a AA_SRR23446543.pair.vcf.gz \
  | /n/home00/bfang/programs/dipcall/dipcall.kit/htsbox bgzip \
  > AA_SRR23446543.dip.vcf.gz

## -- Remove sex chromosomes & summarise --------------------------------------
grep -v '^ScYP8k310HRSCAF43chZ_RagTag\|^SUPER_W_RagTag' AA_SRR23446543.dip.bed \
  > AA_SRR23446543.dip_nosex.bed  # remove sex‑linked regions

/n/home00/bfang/programs/dipcall/dipcall.kit/dipcall-aux.js dipsum \
  AA_SRR23446543.dip_nosex.bed AA_SRR23446543.dip.vcf.gz \
  > Statistics_AA_SRR23446543.txt

############
### PSMC ###
############
cd /n/holyscratch01/edwards_lab/bfang/Scrub_jay/Dipcall_results/AA_SRR23446543

# (Optional) sample‑specific reference; update if needed
Ref_scrubjay=/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/HIC_REFERENCE_ASSEMS/AA_SRR23446543_FSJragtag.v1.fasta

## -- Convert to PSMC input ----------------------------------------------------
seqtk mutfa $Ref_scrubjay <(gzip -dc AA_SRR23446543.dip.vcf.gz | vcf2snp.pl -) \
  | seqtk seq -cM AA_SRR23446543.dip_nosex.bed -l80 \
  | fq2psmcfa - \
  > AA_SRR23446543.psmcfa

## -- Run PSMC -----------------------------------------------------------------
psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o AA_SRR23446543.psmc AA_SRR23446543.psmcfa

## -- Plot Ne trajectory -------------------------------------------------------
psmc_plot.pl -p -T AA_SRR23446543 AA_SRR23446543 AA_SRR23446543.psmc
