#' Calculate percentage of expression per group for each gene.
#'
#'
#' @param data Numeric data.frame, matrix or similar with rows being genes/observations and columns being cells/samples.
#' @param group Vector of length \code{ncol(data)} with group information
#' @param threshold Consider counts above this threshold as "expressed"
#' @param digits Round percentage to this digit
#'
#' @examples
#' ngenes <- 10
#' nsamples <- 1000
#' set.seed(1)
#' data <- matrix(rnbinom(ngenes * nsamples, mu = 0.5, nsamples), ncol = nsamples, byrow = TRUE)
#' rownames(data) <- paste0("gene", 1:nrow(data))
#' group <- rep(LETTERS[1:10], each = nsamples / 10)
#' get_pexpr(data, group)
#'
#' # Also works on Csparse matrix
#' data_sparse <- as(data, "CsparseMatrix")
#' get_pexpr(data_sparse, group)
#'
#' @author Alexander Bender
#'
#' @importFrom Matrix t
#'
#' @export
get_pexpr <- function(data, group, threshold = 0, digits = 2) {
  if (ncol(data) != length(group)) stop("ncol(data) != length(group)")
  if (!is.numeric(threshold) | threshold < 0) stop("threshold must be numeric and > 0")

  datar <- (data >= threshold) * 1
  a <- rowsum(x = t(datar), group = group)
  b <- as.numeric(table(group)[rownames(a)])
  f <- round(100 * t(apply(a, 2, function(x) x / b)), digits = digits)
  f
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

#' Aggregate a SingleCellExperiment assay to pseudobulk sums using only base R
#' @param count_matrix Matrix of (raw) counts to aggregate
#' @param metadata Data.frame with annotations
#' @param aggregate_by Character vector indicating which columns of metadata to use
#' @param sep String used to delimit colnames in aggregated data.frame
#' @examples
#' sce <- mockSCE()
#' sce$group <- rep(c("A", "B"), ncol(sce)/2)
#' sce$group2 <- rep(c("A", "B", "C", "D"), each = ncol(sce)/4)
#' pb <- .aggreg_pseudobulk(count_matrix = assay(sce), metadata = colData(sce), aggregate_by = c("group", "group2"), sep = ".")
#' @export
.aggreg_pseudobulk <- function(count_matrix, metadata, aggregate_by, sep = ".") {
  aggreg <- data.frame(metadata[, aggregate_by, drop = FALSE])

  aggreg_groups <- do.call(rbind, lapply(1:nrow(aggreg), function(z) {
    data.frame(
      concat = paste(aggreg[z, , drop = TRUE], collapse = sep),
      aggreg[z,, drop = FALSE ]
    )
  }))

  aggregated <- t(rowsum(t(count_matrix), group = aggreg_groups$concat))
  aggreg_groups_unique <- unique(aggreg_groups)
  rownames(aggreg_groups_unique) <- aggreg_groups_unique$concat
  tbl <- table(aggreg_groups$concat)
  aggreg_groups_unique$ncells <- as.numeric(tbl[aggreg_groups_unique$concat])
  aggreg_groups_unique$concat <- NULL
  new_ids <- aggreg_groups_unique[colnames(aggregated), , drop = FALSE]
  l <- list(counts = aggregated, metadata = new_ids)

  return(l)
}

#' @importFrom stats as.formula model.matrix
#' @importFrom limma is.fullrank
#' @export
.make_design <- function(main_covariate, other_covariates, metadata) {
  others <- if (!is.null(other_covariates)) paste(other_covariates, collapse = " + ") else ""
  strg <- paste("~ 0", main_covariate, others, sep = " + ")
  f <- as.formula(gsub(" \\+ $", "", strg))
  design <- model.matrix(f, metadata)
  attr(design, "formula") <- paste(f, collapse = "")
  colnames(design) <- gsub(paste0("^", main_covariate), "", colnames(design))

  is_fr <- is.fullrank(design)
  if (!is_fr) stop("Design is not full-rank, meaning you cannot adjust for these covariates!")

  return(design)
}

# An exact copy of scuttle::mockSCE
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
