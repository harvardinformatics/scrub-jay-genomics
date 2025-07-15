# Scrub Jay (Pan)Genomics

Code from the manuscript: [Comparative population pangenomes reveal unexpected complexity and fitness effects of structural variants](https://www.biorxiv.org/content/10.1101/2025.02.11.637762v1)

Structural variants (SVs) are widespread in vertebrate genomes, yet their evolutionary dynamics remain poorly understood. Using 45 long-read de novo genome assemblies and pangenome tools, we analyze SVs within three closely related species of North American jays (Aphelocoma, scrub-jays) displaying a 60-fold range in effective population size. We find rapid evolution of genome architecture, including ∼100 Mb variation in genome size driven by dynamic satellite landscapes with unexpectedly long (> 10 kb) repeat units and widespread variation in gene content, influencing gene expression. SVs exhibit slightly deleterious dynamics modulated by variant length and population size, with strong evidence of adaptive fixation only in large populations. Our results demonstrate how population size shapes the distribution of SVs and the importance of pangenomes to characterizing genomic diversity.

This repository contains analysis code for documentation and reproducibility purposes, it does not contain code that is intended to be reused. 

The sample data is stored in scrubjay-samples.csv in the top-level directory. The rest of the code base is divided into sections:
- assembly: code for producing assemblies of each individual from raw HiFi data
- annotation: code for producing gene and repeat annotations for the reference individual
- graphs: code for producing the pggb graph and associated files, as well as other graph analysis code
- variation: code for producing a final vcf from the pggb graph
- popgen: code for population genetics analysis
- expression: code for expression analysis

