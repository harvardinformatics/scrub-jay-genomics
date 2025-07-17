This directory contains code for expression analysis:

(1) pangene-exp.Rmd contains the final expression analysis code, using the data in `exp-final-files`
(2) rna-seq contains code to run QC (rseqc), as well as kallisto to generate count tables for each RNA-seq sample. For kallisto, we used the merged, non-redundant annotation from TOGA (`rna-seq/kallisto_array_redo.slurm`). We then merged expression info from all samples together into a single table (`collate_kallisto.sh`)