# SigGenes
  
SigGenes automates pairwise differential testing between groups and extraction of signatures.
Documentation will follow. Use at own risk, no warranty.

## Installation

```r
# Install Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install()

# Install limma and circtools
to_install <- c("edgeR", "limma", "Matrix", "S4Vectors", "scuttle", "SummarizedExperiment", "SingleCellExperiment", "atpoint/SigGenes")
BiocManager::install(to_install)
```
