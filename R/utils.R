#' Percentage of expressing cells per group
#'
#' Calculates, for each gene or observation, the percentage of cells or samples
#' in each group that meet the expression threshold.
#'
#' @param data Numeric matrix-like object with genes or observations in rows and
#'   cells or samples in columns. Dense matrices, data frames, and sparse
#'   \code{Matrix} objects are supported.
#' @param group Vector of length \code{ncol(data)} defining the group of each
#'   column.
#' @param threshold Non-negative expression threshold. When \code{threshold = 0},
#'   values greater than zero are considered expressed. Otherwise, values greater
#'   than or equal to \code{threshold} are considered expressed.
#' @param digits Number of decimal places used to round percentages.
#'
#' @return A numeric matrix with genes or observations in rows, groups in
#'   columns, and percentages as values.
#'
#' @examples
#' ngenes <- 10
#' nsamples <- 1000
#'
#' set.seed(1)
#' data <- matrix(
#'   rnbinom(ngenes * nsamples, size = 1, mu = 0.5),
#'   nrow = ngenes
#' )
#' rownames(data) <- paste0("gene", seq_len(ngenes))
#' group <- rep(LETTERS[1:10], each = nsamples / 10)
#'
#' get_pexpr(data, group)
#'
#' # Sparse matrices are also supported
#' data_sparse <- Matrix::Matrix(data, sparse = TRUE)
#' get_pexpr(data_sparse, group)
#'
#' @author Alexander Bender
#'
#' @importFrom Matrix Matrix sparse.model.matrix t
#'
#' @export
get_pexpr <- function(data, group, threshold = 0, digits = 2) {
  if (ncol(data) != length(group))
    stop("ncol(data) != length(group)")

  if (!is.numeric(threshold) || length(threshold) != 1L ||
      is.na(threshold) || threshold < 0)
    stop("threshold must be a single non-negative number")

  if(anyNA(group))
    stop("group must not contain missing values")

  group <- factor(group)
  detected <- if (threshold == 0) data > 0 else data >= threshold

  if (inherits(data, "sparseMatrix")) {
    design <- sparse.model.matrix(~ group - 1)
    out <- detected %*% design
    out <- t(t(out) / tabulate(group, nbins = nlevels(group)))
    colnames(out) <- levels(group)
  } else {
    out <- rowsum(t(detected * 1), group, reorder = FALSE)
    out <- t(out / tabulate(group, nbins = nlevels(group)))
  }

  ro <- round(100 * out, digits)
  ro <- data.frame(as.matrix(ro), check.names = FALSE, check.rows = FALSE)

  return(ro)

}

#' Makes all pairwise contrasts based on an intercept-less design
#' @importFrom utils combn
#' @importFrom limma makeContrasts
#' @keywords internal
.make_contrasts_pairwise <- function(u, design, delim) {
  all_combinations <- combn(u, 2)

  all_contrasts <- lapply(1:ncol(all_combinations), function(cc) {
    one <- all_combinations[1, cc, drop = TRUE]
    two <- all_combinations[2, cc, drop = TRUE]

    p_raw <- paste0(one, "-", two)
    p_valid <- paste0(make.names(one), "-", make.names(two))
    col_valid <- make.names(colnames(design))

    co <- suppressWarnings(makeContrasts(contrasts = p_valid, levels = col_valid))
    rownames(co) <- colnames(design)
    colnames(co) <- paste0(one, delim, two)
    co
  })

  contrasts <- do.call(cbind, all_contrasts)
  colnames(contrasts) <- gsub("   -   ", delim, colnames(contrasts))

  return(contrasts)
}

#' Makes all averaging contrasts based on an intercept-less design
#' @importFrom limma makeContrasts
#' @keywords internal
.make_contrasts_average <- function(u, design) {
  all_contrasts <- lapply(u, function(cc) {
    # potentially invalid
    one <- cc
    two <- setdiff(u, cc)
    two2 <- paste0("(", paste(two, collapse = "+"), ") / ", length(two))

    # valid
    onev <- make.names(one)
    twov <- make.names(two)
    two2v <- paste0("(", paste(twov, collapse = "+"), ") / ", length(twov))
    cnv <- make.names(colnames(design))
    cov <- paste0(onev, " - ", two2v)
    con <- suppressWarnings(makeContrasts(contrasts = cov, levels = cnv))
    con

    rownames(con) <- colnames(design)
    colnames(con) <- cc
    con
  })

  contrasts <- do.call(cbind, all_contrasts)
  contrasts

  return(contrasts)
}

#' Modified from Kolde et al in \code{RobustRankAggreg::rankMatrix()}
#' @param x list with ranked elements
.rankmatrix <- function(x) {
  u <- unique(c(x, recursive = TRUE))
  N <- length(u)

  rmat <- matrix(NA,
    nrow = length(u), ncol = length(x),
    dimnames = list(u, names(x))
  )

  for (i in names(x)) {
    rmat[x[[i]], i] <- seq(1, length(x[[i]]))
  }

  return(rmat)
}

# Convert size factors to edgeR-style normalization factors
.sf2nf <- function(sf, ls) {
  nf <- log(sf / ls)
  nf <- exp(nf - mean(nf))
  nf
}

#' Aggregate counts to pseudobulk
#'
#' Aggregates columns of a count matrix by combinations of metadata variables,
#' returning summed counts and pseudobulk-level metadata.
#'
#' @param count_matrix Numeric matrix-like object with features in rows and cells
#'   in columns.
#' @param metadata Data frame-like object with one row per cell.
#' @param aggregate_by Character vector naming metadata columns used to define
#'   pseudobulk samples.
#' @param sep Character string separating metadata values in pseudobulk names.
#' @param format EIther list or sce. If list will return a SimpleList with the count matrix and metadata,
#'   and if sce then returns a SingleCellExperiment with first assay being the counts and metadata as colData.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{counts}}{Feature-by-pseudobulk matrix of summed counts.}
#'   \item{\code{metadata}}{Pseudobulk metadata including the number of cells.}
#' }
#'
#' @examples
#' sce <- mockSCE()
#' sce$group <- rep(c("A", "B"), ncol(sce) / 2)
#' sce$group2 <- rep(c("A", "B", "C", "D"), each = ncol(sce) / 4)
#'
#' pb <- aggregate_to_pseudobulk(
#'   count_matrix = SummarizedExperiment::assay(sce),
#'   metadata = SummarizedExperiment::colData(sce),
#'   aggregate_by = c("group", "group2")
#' )
#'
#' @importFrom Matrix sparse.model.matrix t
#' @importFrom S4Vectors SimpleList DataFrame
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
aggregate_to_pseudobulk <- function(count_matrix, metadata, aggregate_by, sep = ".", format = c("list", "sce")) {

  if (ncol(count_matrix) != nrow(metadata))
    stop("ncol(count_matrix) != nrow(metadata)")

  aggreg <- data.frame(metadata[, aggregate_by, drop = FALSE])
  if (anyNA(aggreg))
    stop("Aggregation variables must not contain missing values. Check metadata!")

  format <- match.arg(format, choices = c("list", "sce"), several.ok = FALSE)

  ids <- do.call(paste, c(aggreg, sep = sep))
  group <- factor(ids, levels = unique(ids))

  if (inherits(count_matrix, "sparseMatrix")) {
    design <- sparse.model.matrix(~ group - 1)
    counts <- count_matrix %*% design
    } else {
    counts <- t(rowsum(t(count_matrix), group, reorder = FALSE))
  }

  colnames(counts) <- levels(group)
  meta <- aggreg[match(levels(group), ids), , drop = FALSE]
  rownames(meta) <- levels(group)
  meta$ncells <- tabulate(group)

  if(format == "list")
    out <- SimpleList(counts = counts, metadata = meta)

  if(format == "sce"){

    out <- SingleCellExperiment(
      assays = list(counts = counts), colData = DataFrame(meta, check.names = FALSE)
    )

  }

  return(out)
}

#' Create an intercept-free design matrix
#'
#' Builds a model matrix for a main covariate and optional adjustment
#' covariates, then verifies that the resulting design is full rank and has
#' unique coefficient names.
#'
#' @param main_covariate Character string naming the primary covariate.
#' @param other_covariates Optional character vector naming additional
#'   covariates to include in the design.
#' @param metadata Data frame-like object containing the covariates.
#'
#' @return A numeric design matrix. The model formula is stored in the
#'   \code{"formula"} attribute.
#'
#' @examples
#' metadata <- data.frame(
#'   group = factor(rep(c("control", "treated"), each = 4)),
#'   batch = factor(rep(c("A", "B"), 4))
#' )
#'
#' design <- make_design(
#'   main_covariate = "group",
#'   other_covariates = "batch",
#'   metadata = metadata
#' )
#'
#' @importFrom stats as.formula model.matrix
#' @importFrom limma is.fullrank
#'
#' @export
make_design <- function(main_covariate, other_covariates = NULL, metadata) {
  covariates <- c(main_covariate, other_covariates)

  if (!length(main_covariate) || length(main_covariate) != 1L)
    stop("main_covariate must be a single column name")

  if (!all(covariates %in% colnames(metadata)))
    stop("Not all covariates are present in metadata")

  formula <- as.formula(
    paste("~ 0 +", paste(covariates, collapse = " + "))
  )

  design <- model.matrix(formula, metadata)
  attr(design, "formula") <- paste(deparse(formula), collapse = "")

  main_columns <- startsWith(colnames(design), main_covariate)
  colnames(design)[main_columns] <- substring(
    colnames(design)[main_columns],
    nchar(main_covariate) + 1L
  )

  if (anyDuplicated(colnames(design)))
    stop("Design matrix contains duplicated coefficient names")

  if (!is.fullrank(design))
    stop("Design matrix is not full rank")

  design
}

# An exact copy of scuttle::mockSCE so we don't need it as dependency
#' @importFrom stats rnbinom runif
#' @importFrom S4Vectors DataFrame
#' @importFrom SummarizedExperiment colData<- colData
#' @importFrom SingleCellExperiment SingleCellExperiment altExp<-
#' @export
mockSCE <- function(ncells = 200, ngenes = 2000, nspikes = 100) {
  spike.means <- 2^runif(nspikes, 3, 8)
  spike.disp <- 100 / spike.means + 0.5
  spike.data <- matrix(rnbinom(nspikes * ncells,
    mu = spike.means,
    size = 1 / spike.disp
  ), ncol = ncells)
  rownames(spike.data) <- sprintf("Spike_%s", formatC(seq_len(nspikes),
    width = 4, flag = 0
  ))
  cell.means <- 2^runif(ngenes, 2, 10)
  cell.disp <- 100 / cell.means + 0.5
  cell.data <- matrix(rnbinom(ngenes * ncells,
    mu = cell.means,
    size = 1 / cell.disp
  ), ncol = ncells)
  rownames(cell.data) <- sprintf("Gene_%s", formatC(seq_len(ngenes),
    width = 4, flag = 0
  ))
  colnames(cell.data) <- sprintf("Cell_%s", formatC(seq_len(ncells),
    width = 3, flag = 0
  ))
  sce <- SingleCellExperiment(list(counts = cell.data))
  colData(sce) <- cbind(colData(sce), DataFrame(Mutation_Status = sample(c(
    "positive",
    "negative"
  ), ncells, replace = TRUE), Cell_Cycle = sample(c(
    "S",
    "G0", "G1", "G2M"
  ), ncells, replace = TRUE), Treatment = sample(c(
    "treat1",
    "treat2"
  ), ncells, replace = TRUE)))
  colnames(spike.data) <- colnames(sce)
  altExp(sce, "Spikes") <- SingleCellExperiment(list(counts = spike.data))
  sce
}
