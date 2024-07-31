#' Differential expression analysis on single-cell or pseudobulk level using limma-voom.
#'
#' Function accepts a SingleCellExperiment and expects an assay "counts" representing raw counts.
#' It then performs DE analysis forming all possible contrasts with limma-voom \(Law et al. 2014\), using the \code{voomLmFit} function from the `edgeR` package.
#' It is similar to running \code{voom} followed by \(optionally\) \code{duplicateCorrelation} and \code{lmFit}, for details see the details section in \code{edgeR::voomLmFit}.
#' DE testing is done by either contrasting all groupwise combinations, or by testing each group level vs the average of all other group levels.
#' The function can optionally aggregate single-cell data into pseudobulks. The function expects one group
#' to be the main covariate for testing and allows an arbitrary number of other covariates to be adjusted for. Alternatively, the user can specify a blocking
#' factor, for example to account for repeated measures.
#'
#' The output is a list with DE stats per contrast. All possible contrasts are displayed, for example for groups A/B/C it's A-B, A-C, B-C, B-A, C-A, C-B. The
#' individual results for identical contrasts e.g. A-B / B-A can differ if \code{min_pct} is above 0, as each list is filtered that the non-reference expresses
#' the gene with more than \code{min_pct} of cells. If {min_pct} is 0 then these two lists would be the same.
#'
#'
#' @param x A SingleCellExperiment
#' @param use_assay Name of the assay representing the raw counts
#' @param aggregate_by Vector with names of the colData columns to aggregate by. If NULL then
#'   no aggregation is done prior to DE analysis.
#' @param use_existing_sf Logical, whether to use the size factors from \code{sizeFactors(x)} for normalization.
#'   Only applies if no pseudobulk aggregation is performed. If FALSE and no pseudobulk aggregation is performed then
#'   normalization is done via the \code{calcNormFactors} function from `edgeR`.
#'   The default is TRUE as one usually does DE analysis on single-cell level after preprocessing and cluster, so normalization
#'   usually has already been performed. If TRUE and no size factors can be found will run library size normalization with the
#'   \code{librarySizeFactors} function from `scuttle`.
#' @param main_covariate Name of the main covariate for the testing.
#' @param other_covariates Vector with other covariate names to adjust for in DE analysis.
#' @param block Name of a factor to block for.
#' @param min_pct The minimum percentage of cells in at least one of the two group levels per contrast expressing the gene (expressed means counts above 0).
#' @param min_fc Minimum fold change to test against via \code{limma::treat()}.
#' @param use_weights Logical, whether to use empirical sample weights in \code{edgeR::voomLmFit}
#' @param mode Either "pairwise" to run all pairwise comparisons between `main_covariate` levels,
#'   or "average" to test each level against the average of all other levels.
#' @param use_filterByExpr Logical, whether to run \code{filterByExpr} from `edgeR` for prefiltering of counts.
#'   This makes sense when testing for example bulk RNA-seq data.
#' @param return_object Logical, whether to return the DGEList in the output
#' @param delim A delimiter string used in the names of the output list.
#'
#' @examples
#' set.seed(1)
#' sce <- scuttle::mockSCE(ncells=1000, ngenes=1000, nspikes=0)
#' sce$group <- rep(LETTERS[1:4], each=250)
#' sce$donor <- rep(LETTERS[1:4], 250)
#' res_pairwise <- de_limma(x = sce, aggregate_by = c("group", "donor"),
#'                          main_covariate = "group", other_covariates = c("donor"))
#' res_pairwise_block <- de_limma(x = sce, aggregate_by = c("group", "donor"),
#'                                main_covariate = "group", block = "donor")
#' res_average <- de_limma(x = sce, aggregate_by = c("group", "donor"),
#'                         main_covariate = "group", other_covariates = c("donor"),
#'                         mode="average")
#'
#' @author Alexander Bender
#'
#' @references
#' Law, CW, Chen, Y, Shi, W, Smyth, GK (2014).
#' Voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology 15, R29.
#' \url{doi:10.1186/gb-2014-15-2-r29}
#'
#' @seealso \code{\link{voom}}
#'
#' @importFrom SingleCellExperiment colData sizeFactors sizeFactors<-
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom scuttle aggregateAcrossCells librarySizeFactors
#' @importFrom edgeR calcNormFactors DGEList filterByExpr voomLmFit
#' @importFrom stats model.matrix
#' @importFrom limma contrasts.fit is.fullrank makeContrasts treat topTreat
#' @importFrom utils combn
#' @importFrom S4Vectors SimpleList
#' @importFrom methods as is
#'
#' @export
de_limma <- function(x, use_assay = "counts", aggregate_by=NULL, use_existing_sf = TRUE,
                     main_covariate, other_covariates=NULL, block=NULL,
                     min_pct=0, min_fc = 1.0, use_weights = FALSE,
                     mode=c("pairwise", "average"),
                     use_filterByExpr = FALSE, return_object = FALSE, delim="_vs_"){

  # Checks
  is_sce <- is(x, "SingleCellExperiment")
  if(!is_sce) stop("x must be a SingleCellExperiment")

  if(!use_assay %in% assayNames(x))
    stop(paste("x contains no assay named", use_assay))

  if(!is.null(aggregate_by)){
    is_agg <- aggregate_by %in% colnames(colData(x))
    if(!sum(is_agg) == length(aggregate_by)) stop("Not all entries in <aggregate_by> are in colData(x)")
  }

  if(!is.null(block)){
    if(!block %in% colnames(colData(x)))
      stop("Name of block not in colData(x)")
  }

  if(!is(use_existing_sf, "logical"))
    stop("use_existing_sf must be logical")

  is_main <- main_covariate %in% colnames(colData(x))
  if(!is_main) stop("Not all entries in <main_covariate> are in colData(x)")

  if(!is.null(other_covariates)){
    is_other <- other_covariates %in% colnames(colData(x))
    if(!is_other) stop("Not all entries in <other_covariates> are in colData(x)")
  }

  is_pct <- min_pct >= 0 & is.numeric(min_pct)
  if(!is_pct) stop("min_pct must be numeric and not negative")

  is_fc <- min_fc >= 1.0 & is.numeric(min_fc)
  if(!is_fc) stop("is_fc must be numeric and >= 1.0")

  if(!is.logical(use_weights)) stop("use_weights must be logical")

  invisible(match.arg(mode, c("pairwise", "average")))
  mode <- match.arg(mode)

  if(!is(delim, "character")) stop("delim must be a charcter")

  # Aggregate to pseudobulk
  if(!is.null(aggregate_by)){

    pb <- aggregateAcrossCells(x, id=colData(x)[,aggregate_by], use.assay.type="counts")
    pb <- DGEList(counts = assay(pb), samples = data.frame(colData(pb), check.names = FALSE))
    run_norm <- TRUE

  } else {

    w <- which(assayNames(x) == "counts")
    pb <- DGEList(counts = assay(x, use_assay), samples = data.frame(colData(x), check.names = FALSE))
    sf <- suppressWarnings(sizeFactors(x))
    # scran::convertTo()
    if (!is.null(sf)) {
      nf <- log(sf/y$samples$lib.size)
      nf <- exp(nf - mean(nf))
      y$samples$norm.factors <- nf
    }

    if(use_existing_sf){

      if(is.null(sizeFactors(x))){

        sf <- librarySizeFactors(x)
        sizeFactors(x) <- sf
        run_norm <- FALSE

      }

      run_norm <- FALSE

    } else {

      run_norm <- TRUE

    }

  }

  pb$samples <- droplevels.data.frame(pb$samples)
  pb$genes <- NULL

  # Require replication in at least one group
  max_n <- max(table(pb$samples[[main_covariate]]))
  if(max_n == 1) stop("No replication in any group -- DE analysis not possible!")

  # Design, allowing any additional covariates
  others <- if(!is.null(other_covariates)) paste0(" + ", other_covariates) else ""
  f <- as.formula(paste0("~ 0 + ", main_covariate, others))
  design <- model.matrix(f, pb$samples)
  colnames(design) <- gsub(paste0("^", main_covariate), "", colnames(design))

  is_fr <- is.fullrank(design)
  if(!is_fr) stop("Design is not full-rank, meaning you cannot adjust for these covariates!")

  # Optionally prefilter with edgeR
  if(use_filterByExpr){

    keep_fbe <- filterByExpr(pb, group = pb$samples[[main_covariate]])
    pb <- pb[keep_fbe,]


  }

  # Prefilter to genes that in at least one group is expressed by, then normalize
  pexp <- get_pexpr(assay(x, use_assay), group = x[[main_covariate]])

  if(min_pct > 0){

    keep <- apply(pexp >= min_pct, 1, sum) > 0
    pb <- pb[keep,]

  }

  rm(x)

  if(nrow(pb) == 0) return(NULL)
  if(nrow(pb) < 10){
    warning("Fewer than 10 genes (arbitrary threshold to trigger this warning) left after expression filter.")
  }

  # Only run normalization if doing pseudobulk analysis
  if(run_norm){
    pb <- calcNormFactors(pb, method = "TMMwsp")
  }

  # For now we do voomLmFit
  if(is.null(block)){
    blocker <- NULL
  } else {
    blocker <- pb$samples[[block]]
  }

  v <- voomLmFit(counts = pb, design = design, sample.weights = use_weights, block = blocker)

  # Make all pairwise comparisons and make the contrasts.fit
  ux <- pb$samples[[main_covariate]]
  if(is(ux, "factor")){
    u <- levels(droplevels(ux))
  } else {
    u <- unique(as.character(ux))
  }

  # -- Pairwise -- # -----------------------------------------------------------
  if(mode=="pairwise"){

    contrasts <- .make_contrasts_pairwise(u, design, delim)
    v <- contrasts.fit(fit = v, contrasts = contrasts)
    v <- treat(v, fc=min_fc)

    # Extract results per contrast
    iter <- colnames(contrasts)
    de_results <- lapply(iter, function(i){

      s <- strsplit(i, delim)[[1]]
      first <- s[1]
      second <- s[2]

      tt <- topTreat(fit = v, coef = i, number=Inf)
      tt <- tt[order(tt$t, decreasing = TRUE),]
      pe <- pexp[rownames(tt), c(first, second)]
      colnames(pe) <- c("pct.1", "pct.2")
      tt <- cbind(tt, pe)

      tt <- tt[tt$pct.1 >= min_pct | tt$pct.2 >= min_pct,]
      tt$adj.P.Val <- p.adjust(tt$P.Value, "BH")
      tt

    })

    names(de_results) <- iter

  }

  if(mode=="average"){

    contrasts <- .make_contrasts_average(u, design)
    v <- contrasts.fit(fit = v, contrasts = contrasts)
    v <- treat(v, fc=min_fc)

    # Extract significant genes per contrast
    iter <- colnames(contrasts)
    de_results <- lapply(iter, function(i){

      tt <- topTreat(fit = v, coef = i, number=Inf)
      tt <- tt[order(tt$t, decreasing = TRUE),]
      pct.1 <- round(as.numeric(pexp[rownames(tt), i, drop=TRUE]), 2)
      pct.2 <- round(as.numeric(apply(pexp[rownames(tt), setdiff(colnames(pexp), i), drop=FALSE], 1, mean)), 2)
      tt$pct.1 <- pct.1
      tt$pct.2 <- pct.2

      tt <- tt[tt$pct.1 >= min_pct,]
      tt$adj.P.Val <- p.adjust(tt$P.Value, "BH")
      tt

    })

    names(de_results) <- iter

  }

  de_results <- sapply(de_results, function(x) data.frame(Gene = rownames(x), x), simplify = FALSE)
  de_results <- as(de_results, "SimpleList")

  to_return <- SimpleList(results = de_results)

  if(return_object){

    to_return[["DGEList"]] <- pb

  }

  to_return[["params"]] <- list()
  to_return[["params"]]$mode <- mode
  to_return[["params"]]$delim <- delim
  to_return[["params"]]$min_pct <- min_pct

  return(to_return)

}
