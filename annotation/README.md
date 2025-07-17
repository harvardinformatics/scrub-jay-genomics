Code to annotation genes and repeats is contained here. 

(1) Repeats

We used RepeatModeler to de novo identify interspersed repetitive elements in the scrub jay reference genome (`repeatmodeler.slurm`). We used KMC and Satellite Repeat Finder (SRF) to identify tandem repeats in the *combined* assemblies (i.e. all assemblies in the graph concatenated together) (`srf_combined_assem.slurm`). The consensus sequences from RepeatModeler and SRF were combined into one library, which was then provided to RepeatMasker to annotate repeats in the reference genome as well as the haplotype assemblies for all other samples (`repeatmasker.slurm` and `repeatmasker_array.slurm`). We also used the library to annotate repeats in SV regions (`repeatmasker_SVs.slurm`). To quantify satellite abundances in the individual assemblies and the raw reads themselves we mapped the satellites from SRF (`map_sats.slurm` and `map_sats_vs_reads.slurm`)

(2) Genes

We used TOGA to project gene models from chicken GRC7b and zebra finch 1.4 onto the reference scrub jay assembly. 
Subsequently, custom scripts were used to merge redundant models and produce a final, filtered gene set (see: genes/postprocess ): 
- filtered TOGA output BED files to retain only projections that are Intact (I), Partially Intact (PI) or Uncertain Loss (UL) (`filter_toga_beds.sh`)
- convert BED to GTF
- add gene IDs (`add_geneids_to_gtf.py`)
- merged chicken and zebra finch annotations and removed redundant annotations (`gffcompare.slurm`)
- cleaned up labels and ported over CDS and gene annotations from source TOGA annotations (`fixLabels_gffutils.py` and `add_CDS_gene_annots_V2.py`)