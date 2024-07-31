#library(SigGenes)

# Load RNA-seq data for CD4T-, CD8T and naive B cells from Haemopedia:
#counts <- readRDS(paste0(
#            system.file("extdata",package="CreateGeneSignatures"),
#            "/haemopedia_subset.rds"))

counts <- readRDS("inst/extdata/haemopedia_subset.rds")

library(SingleCellExperiment)
coldata <- DataFrame(group = gsub("\\..", "", colnames(counts)), row.names = colnames(counts), check.names = FALSE)
sce <- SingleCellExperiment(assays=list(counts = counts), colData = coldata)

# Pairwise Mode
de <- de_limma(x = sce, use_assay = "counts", aggregate_by = NULL,
               main_covariate = "group", min_fc = 1.5, use_existing_sf = FALSE, return_object = TRUE,
               other_covariates = NULL, block = NULL,
               min_pct = 0, use_weights = FALSE,
               mode = "pairwise", delim = "_vs_",
               use_filterByExpr = FALSE)

signatures <- create_signatures(de, min_prop = 1, n = Inf)

lengths(signatures)

library(pheatmap)
logcpm <- log2(edgeR::cpm(de$DGEList,log=FALSE)+1)

col_order <- unlist(lapply(names(signatures), function(x) which(de$DGEList$samples$group %in% x)))

logcpmZ <- t(scale(t(logcpm[unique(unlist(signatures)),])))
pheatmap(mat=logcpmZ[,col_order],
         show_rownames=FALSE, cluster_rows=FALSE, cluster_cols=FALSE)

# Average Mode
de2 <- de_voom(x = sce, use_assay = "counts", aggregate_by = NULL,
               main_covariate = "group", min_fc = 1.5, use_existing_sf = FALSE, return_object = FALSE,
               other_covariates = NULL, block = NULL,
               min_pct = 0, use_weights = FALSE,
               mode = "average", delim = "_vs_",
               use_filterByExpr = FALSE)

signatures2 <- sapply(de2, function(x){

  xy <- x[x$logFC > 0 & x$adj.P.Val < 0.05,]
  xy <- xy[order(xy$t, decreasing = TRUE),]
  head(rownames(xy), 100)

}, simplify = FALSE)

lengths(signatures2)

logcpmZ2 <- t(scale(t(logcpm[unique(unlist(signatures2)),])))

pheatmap(mat=logcpmZ2[,col_order],
         show_rownames=FALSE, cluster_rows=FALSE, cluster_cols=FALSE)
