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

## Modes

- Default: Uses raw counts and runs limma-voom. If `sizeFactors(x)` are present will use these. Else, 
will use `x$norm.factors` if present, and if none present, will run default edgeR normalization.

- If pseudobulk aggregation, will run calcNormFactors after applying min.pct for at least one group.

- If use_directly, will use the provided assay directly.

CAVE:

- Prefilter
- All ine one model or separeete
