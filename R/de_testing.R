#' Differential expression analysis with limma
#'
#' Performs differential expression analysis on a `SingleCellExperiment` using
#' limma. The input may be aggregated to pseudobulk samples, converted from
#' counts to log2 counts per million, or used directly when the selected assay
#' already contains prefiltered log2-scale expression values.
#'
#' Pairwise comparisons test every level of `main_covariate` against every other
#' level. Average comparisons test each level against the average of all
#' remaining levels.
#'
#' @param x A `SingleCellExperiment`.
#' @param use_assay Character scalar specifying the assay to use. Defaults to
#'   `"counts"`.
#' @param aggregate_by Optional character vector naming columns in `colData(x)`
#'   used to aggregate cells into pseudobulk samples. All model covariates should
#'   be retained by the aggregation.
#' @param use_directly Logical scalar. If `TRUE`, the selected assay is used
#'   directly for testing and is assumed to be prefiltered and on a log2 scale.
#'   Normalization and prefiltering are skipped.
#' @param limma_method Moderation method passed to `limma::treat`. Either
#'   `"trend"` to use intensity-dependent empirical Bayes moderation or
#'   `"ordinary"` for ordinary empirical Bayes moderation.
#' @param prefilter_method Gene-filtering method. One of
#'   `"percent_expressed"`, `"filterByExpr"`, or `"none"`.
#' @param mode Contrast mode. `"pairwise"` performs all pairwise comparisons.
#'   `"average"` compares each group against the average of all remaining
#'   groups.
#' @param main_covariate Character scalar naming the primary grouping variable
#'   in `colData(x)`.
#' @param other_covariates Optional character vector naming additional model
#'   covariates in `colData(x)`.
#' @param min_pct Numeric scalar greater than or equal to zero. Minimum
#'   percentage of expressing cells required in at least one group when
#'   `prefilter_method = "percent_expressed"`.
#' @param min_expr Numeric scalar greater than or equal to zero. Expression
#'   threshold passed to `get_pexpr()`.
#' @param min_fc Numeric scalar greater than or equal to one. Minimum fold-change
#'   threshold on the non-log scale passed to `limma::treat`.
#' @param use_weights Logical scalar. If `TRUE`, sample weights are estimated
#'   with `limma::arrayWeights`.
#' @param delim Character scalar used to separate group names in pairwise
#'   contrast names. Defaults to `"_vs_"`.
#' @param return_object Logical scalar. If `TRUE`, the processed
#'   `SingleCellExperiment` is included in the returned object.
#' @param verbose Logical scalar. If `TRUE`, progress messages are printed.
#'
#' @return A `SimpleList` containing:
#' \itemize{
#'   \item `results`: a `SimpleList` of differential-expression tables, one per
#'     contrast.
#'   \item `params`: parameters and model information used for the analysis.
#'   \item `sce`: the processed `SingleCellExperiment`, included only when
#'     `return_object = TRUE`.
#' }
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' set.seed(1)
#' sce <- mockSCE(ncells = 1000, ngenes = 1000, nspikes = 0)
#' sce$group <- rep(LETTERS[1:4], each = 250)
#' sce$donor <- rep(LETTERS[1:4], 250)
#'
#' res <- de_limma(
#'   x = sce,
#'   aggregate_by = c("group", "donor"),
#'   main_covariate = "group",
#'   other_covariates = "donor",
#'   limma_method = "trend",
#'   mode = "pairwise"
#' )
#'
#' @importFrom SingleCellExperiment sizeFactors
#' @importFrom SummarizedExperiment assay assay<- assayNames assayNames<- colData
#' @importFrom edgeR calcNormFactors cpm filterByExpr
#' @importFrom stats p.adjust
#' @importFrom limma arrayWeights contrasts.fit lmFit treat topTreat
#' @importFrom S4Vectors SimpleList
#' @importFrom methods as is
#' @importFrom Matrix colSums
#'
#' @export
de_limma <- function(
    x, use_assay = "counts", aggregate_by = NULL, use_directly = FALSE,
    limma_method = c("trend", "ordinary"),
    prefilter_method = c("percent_expressed", "filterByExpr", "none"),
    mode = c("pairwise", "average"), main_covariate, other_covariates = NULL,
    min_pct = 0, min_expr = 0, min_fc = 1, use_weights = FALSE, delim = "_vs_",
    return_object = FALSE, verbose = FALSE) {

  if (!is(x, "SingleCellExperiment"))
    stop("x must be a SingleCellExperiment")

  if (!is.character(use_assay) || length(use_assay) != 1L || is.na(use_assay) || !use_assay %in% assayNames(x))
    stop("use_assay must name an assay in x")

  if (!is.null(aggregate_by) && (!is.character(aggregate_by) || length(aggregate_by) == 0L || anyNA(aggregate_by) || !all(aggregate_by %in% colnames(colData(x)))))
    stop("aggregate_by must contain column names from colData(x)")

  if (!is.character(main_covariate) || length(main_covariate) != 1L || is.na(main_covariate) || !main_covariate %in% colnames(colData(x)))
    stop("main_covariate must name one column in colData(x)")

  if (!is.null(other_covariates) && (!is.character(other_covariates) || length(other_covariates) == 0L || anyNA(other_covariates) || !all(other_covariates %in% colnames(colData(x)))))
    stop("other_covariates must contain column names from colData(x)")

  if (!is.numeric(min_pct) || length(min_pct) != 1L || is.na(min_pct) || min_pct < 0)
    stop("min_pct must be a numeric scalar greater than or equal to zero")

  if (!is.numeric(min_expr) || length(min_expr) != 1L || is.na(min_expr) || min_expr < 0)
    stop("min_expr must be a numeric scalar greater than or equal to zero")

  if (!is.numeric(min_fc) || length(min_fc) != 1L || is.na(min_fc) || min_fc < 1)
    stop("min_fc must be a numeric scalar greater than or equal to one")

  if (!is.logical(use_directly) || length(use_directly) != 1L || is.na(use_directly))
    stop("use_directly must be TRUE or FALSE")

  if (!is.logical(use_weights) || length(use_weights) != 1L || is.na(use_weights))
    stop("use_weights must be TRUE or FALSE")

  if (!is.logical(return_object) || length(return_object) != 1L || is.na(return_object))
    stop("return_object must be TRUE or FALSE")

  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose))
    stop("verbose must be TRUE or FALSE")

  if (!is.character(delim) || length(delim) != 1L || is.na(delim) || !nzchar(delim))
    stop("delim must be a non-empty character scalar")

  limma_method <- match.arg(limma_method)
  prefilter_method <- match.arg(prefilter_method)
  mode <- match.arg(mode)

  if (use_directly)
    prefilter_method <- "none"

  mode <- match.arg(mode)

  co <- assay(x, use_assay)
  md <- droplevels(data.frame(colData(x), check.names = FALSE))

  # Aggregate to pseudobulk
  if (!is.null(aggregate_by)) {

    if (verbose) message("Aggregating to pseudobulk")

    y <- aggregate_to_pseudobulk(
      count_matrix = co,
      metadata = md,
      aggregate_by = aggregate_by,
      format = "sce"
    )

    assayNames(y) <- "main"

  } else {

    y <- x
    assay(y, "main") <- assay(y, use_assay)

  }

  rm(co, md)

  design <- make_design(main_covariate, other_covariates, colData(y))

  # Normalization
  y$lib.size <- colSums(assay(y, "main"))

  need_run_norm <- FALSE

  if(!is.null(aggregate_by))
    need_run_norm <- TRUE

  if(is.null(sizeFactors(y)) & is.null(y$norm.factors) & !use_directly)
    need_run_norm <- TRUE

  # Size factors have priority over norm.factors
  if(!is.null(sizeFactors(y)) & !use_directly){

    y$norm.factors <- .sf2nf(sf = sizeFactors(y), ls = y$lib.size)

  }

  # Prefiltering before running norm
  keep <- 1:nrow(y)

  if (prefilter_method == "filterByExpr") {
    if(verbose) message("Run prefiltering")
    keep <- filterByExpr(assay(y, "main"), design = design)
  }

  if (prefilter_method == "percent_expressed" | min_pct > 0) {
    if(verbose) message("Run prefiltering")
    pexp <- get_pexpr(data = assay(x, use_assay), group = colData(x)[[main_covariate]], threshold = min_expr, digits = 2)
    keep <- rowSums(pexp >= min_pct) > 0
  }

  if(sum(keep) < 50)
    warning("Fewer than 50 genes (arbitrary threshold to trigger this message) genes left after prefiltering")

  y <- y[keep,]
  y$lib.size <- colSums(assay(y, "main"))

  # Calculate norm.factors if necessary
  if(use_directly)
    need_run_norm <- FALSE

  if(need_run_norm){
    if(verbose) message("Run normLibSizes")
    y$norm.factors <- as.numeric(calcNormFactors(assay(y, "main")))
  }

  # Now calculate CPMs
  if(!use_directly){

    if(verbose) message("Calculating logcpms")
    assay(y, "main") <- cpm(y = assay(y, "main"), lib.size = y$norm.factors * y$lib.size, log = TRUE)

  }

  # Optionally use sample weights
  if(use_weights){

    if(verbose) message("Calculating sample weights")
    y$sample_weights <- arrayWeights(assay(y, "main"), design = design)

  }

  # Now do testing
  if(verbose) message("Running testing")
  fit <- lmFit(assay(y, "main"), design = design, weights = y$sample_weights)

  group <- colData(y)[[main_covariate]]
  groups <- if (is.factor(group)) levels(droplevels(group)) else unique(as.character(group))

  if (mode == "pairwise") {
    pairs <- combn(groups, 2)
    contrasts <- .make_contrasts_pairwise(groups, design, delim)
  } else {
    contrasts <- .make_contrasts_average(groups, design)
  }

  fit <- contrasts.fit(fit, contrasts)
  fit <- treat(fit, fc = min_fc, trend = (limma_method == "trend"))

  results <- lapply(seq_len(ncol(contrasts)), function(i) {

    tt <- topTreat(fit, coef = i, number = Inf, confint = TRUE)
    tt <- tt[order(tt$t, decreasing = TRUE), , drop = FALSE]

    if (prefilter_method == "percent_expressed" | min_pct > 0) {

      if (mode == "pairwise") {

        pct <- pexp[rownames(tt), pairs[, i], drop = FALSE]
        colnames(pct) <- c("pct.1", "pct.2")

      } else {

        group_i <- colnames(contrasts)[i]

        pct <- cbind(
          pct.1 = pexp[rownames(tt), group_i],
          pct.2 = rowMeans(pexp[rownames(tt), setdiff(groups, group_i), drop = FALSE])
        )

      }

      tt <- cbind(tt, pct)
      tt <- tt[tt$pct.1 >= min_pct | tt$pct.2 >= min_pct, , drop = FALSE]
      tt$adj.P.Val <- p.adjust(tt$P.Value, method = "BH")
    }

    data.frame(Gene = rownames(tt), tt, row.names = NULL)

  })

  names(results) <- colnames(contrasts)

  out <- SimpleList(results = as(results, "SimpleList"))
  out$params <- list(
    mode = mode, delim = delim, min_pct = min_pct, min_expr = min_expr,
    formula = attr(design, "formula"), limma_method = limma_method,
    filtering_method = prefilter_method
  )

  if (return_object) {
    out$sce <- y
  }

  return(out)
}
