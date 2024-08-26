# SigGenes
  
![CI](https://github.com/ATpoint/SigGenes/actions/workflows/ci.yml/badge.svg)
  
## Installation

```r
# Install Bioconductor
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install()

# Install limma and circtools
to_install <- c("edgeR", "limma", "Matrix", "S4Vectors", "scuttle", "SummarizedExperiment", "SingleCellExperiment", "atpoint/SigGenes")
BiocManager::install(to_install)
```
