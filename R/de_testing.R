#' Automated differential expression via the limma package
#'
#' @param x A `SingleCellExperiment`
#' @param use_assay Name of the assay used for testing, "counts" by default
#' @param aggregate_by A character vector of column names in `colData(x)` used for pseudobulk aggregation.
#' @param use_directly Logical, if `TRUE`, use provided assay directly for testing. Assumes data are prefiltered and on log2 scale.
#' @param limma_method Either "trend" or "voom". See details
#' @param prefilter_method The prefiltering strategy, see details
#' @param mode Comparison mode, either "pairwise" (all-vs-all) or `"average"` (one-vs-avg.of.rest)
#' @param main_covariate The main grouping variable in `colData(x)` used for contrasts
#' @param other_covariates Optional character vector of additional covariates to include in the model
#' @param block Optional blocking variable (e.g., donor ID) from `colData(x)` for limma's duplicate correlation
#' @param min_pct Numeric between 0 and 100, minimum percent expression per group required to keep a gene (used with "percent_expressed" method).
#' @param min_fc Minimum absolute fold change threshold (non-log scale) used in `limma::treat`
#' @param use_weights Logical, if `TRUE`, will estimate empirical sample weights
#' @param delim Delimiter used to name pairwise contrasts (default: "_vs_").
#' @param return_object Logical, if `TRUE`, return `DGEList` object and design/contrast info as attributes. Will be the pseudobulked one if that option was chosen
#' @param verbose Logical, if `TRUE`, print progress messages

#' @examples
#'
#' # Testing on pseudobulk level with limma-voom
#' library(SingleCellExperiment)
#' set.seed(1)
#' sce <- mockSCE(ncells = 1000, ngenes = 1000, nspikes = 0)
#' sce$group <- rep(LETTERS[1:4], each = 250)
#' sce$donor <- rep(LETTERS[1:4], 250)
#' sizeFactors(sce) <- rnorm(ncol(sce), mean = 1, sd = .2)
#'
#' res_pairwise <- de_limma(
#'   x = sce, aggregate_by = c("group", "donor"), mode = "pairwise",
#'   main_covariate = "group", other_covariates = c("donor"), limma_method = "trend"
#' )
#'
#' @author Alexander Bender
#'
#' @importFrom SingleCellExperiment colData sizeFactors sizeFactors<-
#' @importFrom SummarizedExperiment assay assayNames colData colData<-
#' @importFrom edgeR calcNormFactors cpm DGEList filterByExpr voomLmFit
#' @importFrom stats median p.adjust
#' @importFrom limma arrayWeights contrasts.fit duplicateCorrelation is.fullrank makeContrasts lmFit treat topTreat
#' @importFrom S4Vectors SimpleList
#' @importFrom methods as is
#'
#' @export
de_limma <- function(
    x, use_assay = "counts", aggregate_by = NULL, use_directly = FALSE,
    limma_method = c("trend", "voom"), prefilter_method = c("percent_expressed", "filterByExpr", "none"),
    mode = c("pairwise", "average"), main_covariate, other_covariates = NULL, block = NULL,
    min_pct = 0, min_fc = 1.0, use_weights = FALSE, delim = "_vs_", return_object = FALSE, verbose = FALSE) {
  # Checks
  is_sce <- is(x, "SingleCellExperiment")
  if (!is_sce) stop("x must be a SingleCellExperiment")

  if (!use_assay %in% assayNames(x)) {
    stop(paste("x contains no assay named", use_assay))
  }

  if (!is.null(aggregate_by)) {
    is_agg <- aggregate_by %in% colnames(colData(x))
    if (!sum(is_agg) == length(aggregate_by)) stop("Not all entries in aggregate_by are in colData(x)")
  }

  if (!is(use_directly, "logical")) stop("use_directly must be logical")

  is_too_long_main <- length(main_covariate)
  if (is_too_long_main > 1) stop("main_covariate must be length 1")

  if (!main_covariate %in% colnames(colData(x))) {
    stop("main_covariate not in colData(x)")
  }

  if (!is.null(other_covariates)) {
    if (!other_covariates %in% colnames(colData(x))) {
      stop("other_covariates not in colData(x)")
    }
  }

  if (!is.null(block)) {
    if (!block %in% colnames(colData(x))) {
      stop("Name of block not in colData(x)")
    }
  }

  if (!is.numeric(min_pct) | min_pct < 0 | min_pct > 100) {
    stop("min_pct must be numeric and between 0 and 100")
  }

  if (min_fc < 1.0) stop("min_fc must be >= 1, it's not logFC!")

  if (!is.logical(use_weights)) stop("use_weights must be logical")

  limma_method <- match.arg(limma_method)
  limma_method <- match.arg(limma_method, choices = c("voom", "trend"))

  trend <- limma_method == "trend"

  mode <- match.arg(mode)
  mode <- match.arg(mode, c("pairwise", "average"))

  prefilter_method <- match.arg(prefilter_method)
  prefilter_method <- match.arg(prefilter_method, c("filterByExpr", "percent_expressed", "none"))

  if (!is.null(aggregate_by) & use_directly) stop("aggregate_by and use_directly cannot be used both at the same time")

  colData(x) <- droplevels.data.frame(colData(x))

  max_n <- max(table(x[[main_covariate]]))
  if (max_n == 1) stop("No replication in any group -- DE analysis not possible!")

  run_calc_norm_factors <- FALSE

  if (!is.null(aggregate_by)) {
    if (verbose) message("Running pseudobulk aggregation")

    tmp <- SigGenes:::.aggreg_pseudobulk(count_matrix = assay(x, use_assay), metadata = colData(x), aggregate_by = aggregate_by, sep = ".")
    y <- DGEList(counts = tmp$counts, samples = data.frame(tmp$metadata, check.names = FALSE))
    rm(tmp)
    run_calc_norm_factors <- TRUE
  }

  if (use_directly) {
    if (verbose) message("Using the assay directly for differential testing")
    y <- DGEList(counts = assay(x, use_assay), samples = data.frame(colData(x), check.names = FALSE))
    limma_method <- "trend"
  }

  # From here on we go with the DGEList format
  if (!use_directly & is.null(aggregate_by)) {
    y <- DGEList(counts = assay(x, use_assay), samples = data.frame(colData(x), check.names = FALSE))

    # Presence of size factors has priority
    y$samples$sizeFactor <- NULL

    if (!is.null(sizeFactors(x))) {
      if (verbose) message("Using existing size factors for normalization")
      y$samples$size_factors <- as.numeric(sizeFactors(x))
      y$samples$norm.factors <- .sf2nf(y$samples$size_factors, y$samples$lib.size)
    } else {
      # If no size factors but norm.factors are present use these
      nf1 <- y$samples$norm.factors.1
      if (!is.null(nf1)) {
        if (verbose) message("Using existing normalization factors + lib.size for normalization")
        y$samples$norm.factors <- nf1
        y$samples$norm.factors.1 <- NULL
      } else {
        run_calc_norm_factors <- TRUE
      }
    }
  }

  y$genes <- NULL

  # Make the design, always putting the main covariate first without intercept
  design <- SigGenes:::.make_design(main_covariate, other_covariates, y$samples)

  # Prefiltering choices
  if (prefilter_method == "filterByExpr") {
    if (verbose) message("Using filterByExpr for prefiltering")
    keep <- filterByExpr(y, design = design)
  }

  pexp <- NULL
  if (prefilter_method == "percent_expressed") {
    if (verbose) message("Using percent expression for prefiltering")
    pexp <- get_pexpr(assay(x, use_assay), group = colData(x)[[main_covariate]])

    if (min_pct > 0) {
      keep <- apply(pexp >= min_pct, 1, sum) > 0
    } else {
      keep <- 1:nrow(y)
    }
  }

  if (prefilter_method == "none") {
    if (verbose) message("Using no prefiltering")
    keep <- 1:nrow(y)
  }

  y <- y[keep, ]

  # Normalize with default method if necessary
  if (run_calc_norm_factors) {
    if (verbose) message("Running default edgeR normalization")
    y <- calcNormFactors(object = y)
  }

  if (nrow(y) == 0) {
    warning("No genes left after prefilter! -- Returning empty list")
    empty_list <- SimpleList()
    return(empty_list)
  }

  if (nrow(y) < 10) {
    warning("Fewer than 10 rows (arbitrary threshold to trigger this warning) left after prefilter")
  }

  # For duplicate correlation
  if (is.null(block)) {
    blocker <- NULL
  } else {
    blocker <- y$samples[[block]]
  }

  # Use voom
  if (limma_method == "voom") {
    if (verbose) message("Running voomLmFit")
    fit <- voomLmFit(counts = y, design = design, sample.weights = use_weights, block = blocker)
  } else {
    aw <- NULL

    if (use_weights) {
      if (verbose) message("Estimating sample weights")
      aw <- arrayWeights(object = y$counts, design = design)
    }

    # For trend, convert to logcpm
    if (!use_directly) {
      if (verbose) message("Converting counts to logCPM with edgeR")
      lcpm <- cpm(y, log = TRUE)
    } else {
      # If using assay directly, the input assay is used
      lcpm <- y$counts
    }

    dconsensus <- NULL
    if (!is.null(block)) {
      if (nrow(lcpm) > 5000) {
        set.seed(1)
        spl <- sample(x = 1:nrow(lcpm), size = 5000, replace = FALSE)
      } else {
        spl <- 1:nrow(lcpm)
      }

      if (verbose) message("Running duplicateCorrelation")
      dcor <- duplicateCorrelation(object = lcpm[spl, ], design = design, block = blocker, weights = aw)
      dconsensus <- dcor$consensus
    }

    if (verbose) message("Running lmFit")
    fit <- lmFit(object = lcpm, design = design, weights = aw, correlation = dconsensus)
  }

  # Make all contrasts and fit it
  ux <- y$samples[[main_covariate]]
  if (is(ux, "factor")) {
    u <- levels(droplevels(ux))
  } else {
    u <- unique(as.character(ux))
  }

  # -- Pairwise -- # -----------------------------------------------------------
  if (mode == "pairwise") {
    if (verbose) message("Running all-vs-all pairwise testing")
    contrasts <- .make_contrasts_pairwise(u, design, delim)
    fit <- contrasts.fit(fit = fit, contrasts = contrasts)
    fit <- treat(fit, fc = min_fc, trend = trend)

    # Extract results per contrast
    iter <- colnames(contrasts)
    de_results <- lapply(iter, function(i) {
      s <- strsplit(i, delim)[[1]]
      first <- s[1]
      second <- s[2]

      tt <- topTreat(fit = fit, coef = i, number = Inf, confint = TRUE)
      tt <- tt[order(tt$t, decreasing = TRUE), ]

      if (!is.null(pexp)) {
        pe <- pexp[rownames(tt), c(first, second)]
        colnames(pe) <- c("pct.1", "pct.2")
        tt <- cbind(tt, pe)
        tt <- tt[tt$pct.1 >= min_pct | tt$pct.2 >= min_pct, ]
        tt$adj.P.Val <- p.adjust(tt$P.Value, "BH")
      }

      tt
    })

    names(de_results) <- iter
  }

  if (mode == "average") {
    if (verbose) message("Running one-vs-average testing")
    contrasts <- .make_contrasts_average(u, design)
    fit <- contrasts.fit(fit = fit, contrasts = contrasts)
    fit <- treat(fit, fc = min_fc)

    # Extract significant genes per contrast
    iter <- colnames(contrasts)
    de_results <- lapply(iter, function(i) {
      tt <- topTreat(fit = fit, coef = i, number = Inf)
      tt <- tt[order(tt$t, decreasing = TRUE), ]

      if (!is.null(pexp)) {
        pct.1 <- round(as.numeric(pexp[rownames(tt), i, drop = TRUE]), 2)
        pct.2 <- round(as.numeric(apply(pexp[rownames(tt), setdiff(colnames(pexp), i), drop = FALSE], 1, mean)), 2)
        tt$pct.1 <- pct.1
        tt$pct.2 <- pct.2
        tt <- tt[tt$pct.1 >= min_pct, ]
        tt$adj.P.Val <- p.adjust(tt$P.Value, "BH")
      }

      tt
    })

    names(de_results) <- iter
  }

  de_results <- sapply(de_results, function(x) data.frame(Gene = rownames(x), x), simplify = FALSE)
  de_results <- as(de_results, "SimpleList")

  to_return <- SimpleList(results = de_results)

  if (return_object) {
    y$samples$weights <- aw
    y$design <- design
    y$contrasts <- contrasts
    to_return[["DGEList"]] <- y
  }

  # Return params as attributes
  to_return[["params"]] <- list()
  to_return[["params"]]$mode <- mode
  to_return[["params"]]$delim <- delim
  to_return[["params"]]$min_pct <- min_pct
  to_return[["params"]]$formula <- attr(design, "formula")

  if (verbose) message("Done!")
  return(to_return)
}
