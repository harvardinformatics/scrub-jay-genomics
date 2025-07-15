# ==============================================================================
#  Scrub-Jay Variant Analysis Pipeline (SFS, Neutrality Tests, DFE)
# ------------------------------------------------------------------------------
#  Description:
#    ‑ Computes summary statistics from PGGB‐derived VCF tables
#    ‑ Builds Site Frequency Spectra (SFS) for three scrub‐jay species
#    ‑ Performs neutrality tests (DoS / α) by genomic feature and variant class
#    ‑ Prepares input tables for downstream DFE inference (fastDFE / anavar)
# ==============================================================================

# -----------------------------------------------------------------------------
# 1.  Set‑up & Libraries
# -----------------------------------------------------------------------------

library(data.table)   # Fast I/O + efficient data manipulation
library(ggplot2)      # Plotting
library(cowplot)      # Multi‑panel figure assembly
library(parallel)     # Multi‑core helpers (mclapply)
library(bedtoolsr)    # BEDTools wrapper for R (not used directly below)

setwd("/n/netscratch/edwards_lab/Lab/bohao/Scrub_jay/SFS_2024Sep")
CORES <- 30           # Number of CPU cores for parallel string ops

# -----------------------------------------------------------------------------
# 2.  Reference & Chromosome Lists
# -----------------------------------------------------------------------------

## (a) Complete reference .fai (31 chromosomes)
fai <- fread("/n/holylfs05/LABS/informatics/Everyone/scrubjay/scrub-jay-genomics/workflow/results/HIC_REFERENCE_ASSEMS/AW_365336_FSJragtag.v1.fasta.fai")

## (b) Chromosomes retained in the PGGB graph
chrs_pggb <- gsub("aphWoo1#REF#", "", fread("/n/netscratch/edwards_lab/Lab/bohao/Scrub_jay/SFS_Nov2023/Chrs_Comparative_Nov7_2023.txt", header = FALSE)$V1)[!chrs_pggb %like% "ch17"]
chrs_pggb_auto <- chrs_pggb[!(chrs_pggb %like% "_W_" | chrs_pggb %like% "chZ")]  # autosomes only

# -----------------------------------------------------------------------------
# 3.  Load & Clean PGGB Variant Table
# -----------------------------------------------------------------------------

tab <- fread("/n/holyscratch01/edwards_lab/Users/bfang/Assemblies_SJ/VCF_table_Jan2024/variation_final_v2.tab")

## Strip reference prefix from chrom column
tab[, chrom := gsub("aphWoo1#REF#", "", chrom)]

## --- Harmonise variant subtype labels ---------------------------------------

tab[subtype == "Complex_NV" & type == "MNP", subtype := "MNP"]  # treat complex NVs as MNPs

## --- Fix column name typos ---------------------------------------------------

tab <- tab[, setnames(.SD, old = c("inv", "polarized", "aa"), new = c("inv_1", "polarized_1", "aa_1"))]
tab <- tab[, setnames(.SD, old = c("inv_1", "polarized_1", "aa_1"), new = c("polarized", "aa", "inv"))]

# -----------------------------------------------------------------------------
# 4.  Genomic Feature Annotation (cnee, exon, intron, intergenic)
# -----------------------------------------------------------------------------

tab[, overlap_new := overlap]                                         # copy original column

tab[overlap_new %like% "cnee",      overlap_new := "cnee"]
tab[overlap_new %like% "exon",      overlap_new := "exon"]
tab[overlap_new %like% "intron",    overlap_new := "intron"]
tab[overlap_new %like% ",",        overlap_new := "intergenic"]
tab[overlap_new ==     "ce",       overlap_new := "intergenic"]

# -----------------------------------------------------------------------------
# 5.  Allele Counts Per Population (haplotype‑aware)
# -----------------------------------------------------------------------------

## Column slices for each species
geno_AW <- tab[, 56:70]
geno_AC <- tab[, 27:40]
geno_AI <- tab[, 41:55]

### Helper: split diploid genotypes into haplotype vectors --------------------

FUN_get_hap_geno <- function(geno) {
  geno <- geno[, mclapply(.SD, function(x) gsub("/", "|", x), mc.cores = CORES)]
  geno <- cbind(
    # Hap1 columns
    geno[, mclapply(.SD, function(x) sapply(strsplit(x, "\\|"), "[", 1), mc.cores = CORES)][, setnames(.SD, new = paste0(names(.SD), "_hap1"))],
    # Hap2 columns
    geno[, mclapply(.SD, function(x) sapply(strsplit(x, "\\|"), "[", 2), mc.cores = CORES)][, setnames(.SD, new = paste0(names(.SD), "_hap2"))]
  )
  geno
}

## Expand to haplotypes
geno_AW_haps <- FUN_get_hap_geno(geno_AW)
geno_AC_haps <- FUN_get_hap_geno(geno_AC)
geno_AI_haps <- FUN_get_hap_geno(geno_AI)

### Allele count (per site, per species) -----------------------------------

tab$allele_count_AW <- apply(geno_AW_haps, 1, function(x) length(unique(x[x != "."])))
tab$allele_count_AC <- apply(geno_AC_haps, 1, function(x) length(unique(x[x != "."])))
tab$allele_count_AI <- apply(geno_AI_haps, 1, function(x) length(unique(x[x != "."])))

# -----------------------------------------------------------------------------
# 6.  Variant Counts by Species & Allelicity (bi‑ and multi‑allelic)
# -----------------------------------------------------------------------------

## Unique key (chrom + coordinates + subtype + per‑species DAC)

tab_1 <- tab[, unique(.SD, by = c("chrom", "bedStart", "bedEnd", "subtype", "AW_DAC", "AC_DAC", "AI_DAC"))]

### (a) Bi‑allelic sites -------------------------------------------------------

tab_3_species <- merge(
  merge(
    tab_1[allele_count == 2][AW_DAC != 0 & AW_DAC != AW_AN][, .(AW = .N), by = subtype],
    tab_1[allele_count == 2][AC_DAC != 0 & AC_DAC != AC_AN][, .(AC = .N), by = subtype],
    all = TRUE
  ),
  tab_1[allele_count == 2][AI_DAC != 0 & AI_DAC != AI_AN][, .(AI = .N), by = subtype],
  all = TRUE
)

# Grand total across populations

tab_3_species_sum <- rbind(
  tab_1[allele_count == 2][AW_DAC != 0 & AW_DAC != AW_AN],
  tab_1[allele_count == 2][AC_DAC != 0 & AC_DAC != AC_AN],
  tab_1[allele_count == 2][AI_DAC != 0 & AI_DAC != AI_AN]
)[order(chrom, bedStart, bedEnd)][, .(All = .N), by = subtype]

tab_all_biallelic <- merge(tab_3_species_sum, tab_3_species, all = TRUE)[order(nchar(subtype))][, allele_count := "Bi‑allelic"]

### (b) Multi‑allelic sites ----------------------------------------------------

tab_3_species_multiallelic <- merge(
  merge(
    tab_1[allele_count_AW > 2][, .(AW = .N), by = subtype],
    tab_1[allele_count_AC > 2][, .(AC = .N), by = subtype],
    all = TRUE
  ),
  tab_1[allele_count_AI > 2][, .(AI = .N), by = subtype],
  all = TRUE
)

tab_3_species_multiallelic_sum <- rbind(
  tab_1[allele_count_AW > 2],
  tab_1[allele_count_AC > 2],
  tab_1[allele_count_AI > 2]
)[order(chrom, bedStart, bedEnd)][, .(All = .N), by = subtype]

tab_all_multiallelic <- merge(tab_3_species_multiallelic_sum, tab_3_species_multiallelic, all = TRUE)[order(nchar(subtype))][, allele_count := "Multi‑allelic"]

tab_all <- rbind(tab_all_biallelic, tab_all_multiallelic)[, .(allele_count, subtype, All, AW, AC, AI)]

# fwrite(tab_all, file = "Sum_Count_PGGB_3sps.txt", sep = "\t")

# -----------------------------------------------------------------------------
# 7.  Site Frequency Spectrum (SFS)
# -----------------------------------------------------------------------------

# Keep bi‑allelic, polarized variants; remove complex SV labels
DAC <- tab_1[allele_count == 2 & polarized == TRUE][!subtype %like% "omplex"][, All_DAC := AW_DAC + AC_DAC + AI_DAC]

## Coarse variant classes ------------------------------------------------------
DAC[subtype == "SNP",              type := "SNP"]
DAC[subtype == "MNP",              type := "MNP"]
DAC[subtype %in% c("INS", "DEL"), type := "INDEL"]
DAC[subtype %like% "SV",            type := "SV"]

### SFS per population / per type --------------------------------------------

sfs_AW <- DAC[(AW_AN) >= 30][, .(frequency = .N), by = .(type, overlap_new, AW_DAC)][, setnames(.SD, "AW_DAC", "DAC")]
sfs_AC <- DAC[(AC_AN) >= 28][, .(frequency = .N), by = .(type, overlap_new, AC_DAC)][, setnames(.SD, "AC_DAC", "DAC")]
sfs_AI <- DAC[(AI_AN) >= 30][, .(frequency = .N), by = .(type, overlap_new, AI_DAC)][, setnames(.SD, "AI_DAC", "DAC")]

## Convert to proportional spectra ---------------

sfs_AW_prop <- rbindlist(lapply(split(sfs_AW[!DAC %in% c(0, 30)], by = c("type", "overlap_new")),
                                function(x) x[, proportion := frequency / sum(frequency)][order(type, overlap_new, DAC)][, population := "AW"]))

sfs_AC_prop <- rbindlist(lapply(split(sfs_AC[!DAC %in% c(0, 28)], by = c("type", "overlap_new")),
                                function(x) x[, proportion := frequency / sum(frequency)][order(type, overlap_new, DAC)][, population := "AC"]))

sfs_AI_prop <- rbindlist(lapply(split(sfs_AI[!DAC %in% c(0, 30)], by = c("type", "overlap_new")),
                                function(x) x[, proportion := frequency / sum(frequency)][order(type, overlap_new, DAC)][, population := "AI"]))

sfs <- rbind(sfs_AW_prop, sfs_AC_prop, sfs_AI_prop)
sfs <- sfs[, c(6, 1:5)][order(population, type, overlap_new, DAC)]

# fwrite(sfs, "SFS.txt", sep = "\t")

### SFS per population (species) / per subtype -----------------------------------------

sfs_AW_subtype <- DAC[(AW_AN) >= 30][, .(frequency = .N), by = .(subtype, overlap_new, AW_DAC)][, setnames(.SD, "AW_DAC", "DAC")]
sfs_AC_subtype <- DAC[(AC_AN) >= 28][, .(frequency = .N), by = .(subtype, overlap_new, AC_DAC)][, setnames(.SD, "AC_DAC", "DAC")]
sfs_AI_subtype <- DAC[(AI_AN) >= 30][, .(frequency = .N), by = .(subtype, overlap_new, AI_DAC)][, setnames(.SD, "AI_DAC", "DAC")]

sfs_AW_prop_subtype <- rbindlist(lapply(split(sfs_AW_subtype[!DAC %in% c(0, 30)], by = c("subtype", "overlap_new")),
                                        function(x) x[, proportion := frequency / sum(frequency)][order(subtype, overlap_new, DAC)][, population := "AW"]))

sfs_AC_prop_subtype <- rbindlist(lapply(split(sfs_AC_subtype[!DAC %in% c(0, 28)], by = c("subtype", "overlap_new")),
                                        function(x) x[, proportion := frequency / sum(frequency)][order(subtype, overlap_new, DAC)][, population := "AC"]))

sfs_AI_prop_subtype <- rbindlist(lapply(split(sfs_AI_subtype[!DAC %in% c(0, 30)], by = c("subtype", "overlap_new")),
                                        function(x) x[, proportion := frequency / sum(frequency)][order(subtype, overlap_new, DAC)][, population := "AI"]))

sfs_subtype <- rbind(sfs_AW_prop_subtype, sfs_AC_prop_subtype, sfs_AI_prop_subtype)
sfs_subtype <- sfs_subtype[, c(6, 1:5)][order(population, subtype, overlap_new, DAC)]

# -----------------------------------------------------------------------------
# 8.  Neutrality Tests (DoS & α) – Species‑specific variants
# -----------------------------------------------------------------------------

# alpha = 1-((Ds*Pn)/(Dn*Ps)) # exclude low DAC <= 4
# DoS = (Dn/(Dn + Ds))-(Pn/(Pn + Ps))

# Dn = fixed_variant
# Ds = fixed_SNP
# Pn = poly_variant
# Ps = poly_SNP

# Filters: retain sites segregating in only one species; exclude DAC ≤ 4

sfs_AW <- DAC[(AW_AN) >= 30][AI_DAC == 0][AC_DAC == 0][, .(frequency = .N), by = .(type, overlap_new, AW_DAC)][, setnames(.SD, "AW_DAC", "DAC")]
sfs_AC <- DAC[(AC_AN) >= 28][AW_DAC == 0][AI_DAC == 0][, .(frequency = .N), by = .(type, overlap_new, AC_DAC)][, setnames(.SD, "AC_DAC", "DAC")]
sfs_AI <- DAC[(AI_AN) >= 30][AW_DAC == 0][AC_DAC == 0][, .(frequency = .N), by = .(type, overlap_new, AI_DAC)][, setnames(.SD, "AI_DAC", "DAC")]

sfs_fastDFE_all <- rbind(
  sfs_AW[, species := "AW"],
  sfs_AC[, species := "AC"],
  sfs_AI[, species := "AI"]
)

## Flag polymorphic classes ----------------------------------------------------

tab_filtered <- sfs_fastDFE_all[DAC != 0][
  , seg_type := fifelse(species == "AW" & DAC == 30, "fixed",
                        fifelse(species == "AC" & DAC == 28, "fixed",
                                fifelse(species == "AI" & DAC == 30, "fixed",
                                        fifelse(DAC <= 4, "low_freq", "poly"))))][]

tab_summary <- tab_filtered[, .(count = sum(frequency)), by = .(species, overlap_new, type, seg_type)]

tab_wider <- dcast(tab_summary, species + overlap_new ~ seg_type + type, value.var = "count", fill = 0)

## Ensure numeric columns ------------------------------------------------------

tab_wider[, (3:ncol(tab_wider)) := lapply(.SD, as.numeric), .SDcols = 3:ncol(tab_wider)]

## Bring in intergenic SNP counts (serves as neutral reference) ---------------

tab_wider_1 <- rbindlist(lapply(split(tab_wider, by = "species"), function(x) {
  data.table(x, x[overlap_new == "intergenic", .(fixed_SNP_intergenic = fixed_SNP,
                                                 poly_SNP_intergenic  = poly_SNP)])
}))

## DoS & α calculations --------------------------------------------------------

neutrality_test_subtype_overlap <- tab_wider_1[, .(
  alpha_SNP   = 1 - ((fixed_SNP_intergenic * poly_SNP)   / (poly_SNP_intergenic * fixed_SNP)),
  alpha_INDEL = 1 - ((fixed_SNP_intergenic * poly_INDEL) / (poly_SNP_intergenic * fixed_INDEL)),
  alpha_SV    = 1 - ((fixed_SNP_intergenic * poly_SV)    / (poly_SNP_intergenic * fixed_SV)),
  
  DoS_SNP   = (fixed_SNP   / (fixed_SNP   + fixed_SNP_intergenic)) - (poly_SNP   / (poly_SNP   + poly_SNP_intergenic)),
  DoS_INDEL = (fixed_INDEL / (fixed_INDEL + fixed_SNP_intergenic)) - (poly_INDEL / (poly_INDEL + poly_SNP_intergenic)),
  DoS_SV    = (fixed_SV    / (fixed_SV    + fixed_SNP_intergenic)) - (poly_SV    / (poly_SV    + poly_SNP_intergenic))
), by = .(species, overlap_new)]

## Round for neat output

neutrality_test_subtype_overlap[, (3:ncol(neutrality_test_subtype_overlap)) := lapply(.SD, round, 5), .SDcols = 3:ncol(neutrality_test_subtype_overlap)]

fwrite(neutrality_test_subtype_overlap,
       "Neutral_test_Type_2025May30_species_specific_sites.txt",
       sep = "\t", col.names = TRUE)

# -----------------------------------------------------------------------------
# 9.  Plotting: DoS & α by Variant Class
# -----------------------------------------------------------------------------

DOS_plot <- melt(neutrality_test_subtype_overlap, id.vars = c("species", "overlap_new"))
DOS_plot[, `:=`(test    = sapply(strsplit(variable, "_"), "[", 1),
                variant = sapply(strsplit(variable, "_"), function(x) paste(x[-1], collapse = "_")))]

DOS_plot[, variant := factor(variant, levels = DOS_plot[, .N, by = variant]$variant)]
DOS_plot[, overlap_new := factor(overlap_new, levels = c("cnee", "exon", "intron", "intergenic", "promoter"))]

## (a) DoS ---------------------------------------------------------------------

p_DOS <- ggplot(DOS_plot[test == "DoS" & overlap_new != "promoter"],
                aes(overlap_new, value, colour = species)) +
  facet_wrap(~ variant, scales = "free", nrow = 3) +
  geom_point(size = 4, alpha = 0.8) +
  geom_line(aes(group = species)) +
  theme_minimal(base_size = 16) +
  theme(aspect.ratio = 1 / 1.618,
        legend.position = "top",
        plot.title       = element_text(face = "bold"),
        axis.line        = element_line(linewidth = 0.3, colour = "#333333"),
        axis.ticks       = element_line(linewidth = 0.3, colour = "#333333"),
        axis.text        = element_text(colour = "#333333"),
        strip.background = element_blank(),
        panel.background = element_blank(),
        plot.background  = element_blank(),
        panel.spacing    = unit(1, "lines"),
        panel.grid       = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title.x      = element_text(margin = margin(t = 10, unit = "pt")),
        axis.title.y      = element_text(margin = margin(r = 10, unit = "pt"))) +
  scale_colour_manual(values = c("AC" = "#a9adb8", "AI" = "#114ef2", "AW" = "#8aa5ee")) +
  labs(x = NULL, y = "Direction of selection", colour = NULL)

## (b) α ----------------------------------------------------------------------

p_alpha <- ggplot(DOS_plot[test == "alpha" & overlap_new != "promoter"],
                  aes(overlap_new, value, colour = species)) +
  facet_wrap(~ variant, scales = "free", nrow = 3) +
  geom_point(size = 4, alpha = 0.8) +
  geom_line(aes(group = species)) +
  theme_minimal(base_size = 16) +
  theme(aspect.ratio = 1 / 1.618,
        legend.position = "top",
        plot.title       = element_text(face = "bold"),
        axis.line        = element_line(linewidth = 0.3, colour = "#333333"),
        axis.ticks       = element_line(linewidth = 0.3, colour = "#333333"),
        axis.text        = element_text(colour = "#333333"),
        strip.background = element_blank(),
        panel.background = element_blank(),
        plot.background  = element_blank(),
        panel.spacing    = unit(1, "lines"),
        panel.grid       = element_blank(),
        axis.ticks.length = unit(0.2, "cm"),
        axis.title.x      = element_text(margin = margin(t = 10, unit = "pt")),
        axis.title.y      = element_text(margin = margin(r = 10, unit = "pt"))) +
  scale_colour_manual(values = c("AC" = "#a9adb8", "AI" = "#114ef2", "AW" = "#8aa5ee")) +
  labs(x = NULL, y = "alpha", colour = NULL)

## Combine panels & save -------------------------------------------------------

p_test <- plot_grid(p_alpha, p_DOS, align = "hv", ncol = 2)

ggsave(plot = p_test,
       filename = "Figures_Jan2025/Neutrality_test_Type_DAC4_2025May30_segregating_in_single_species.pdf",
       width = 8, height = 8)

# -----------------------------------------------------------------------------
# 10.  Distribution of Fitness Effects (DFE)
# -----------------------------------------------------------------------------

# For full DFE workflow (fastDFE + anavar), please see:
#   https://github.com/fangbohao/house-finch-pangenome/tree/main/5_Distribution_of_fitness_effects
# -----------------------------------------------------------------------------
