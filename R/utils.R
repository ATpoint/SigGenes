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
#' data <- matrix(rnbinom(ngenes*nsamples, mu=0.5, nsamples), ncol=nsamples, byrow=TRUE)
#' rownames(data) <- paste0("gene", 1:nrow(data))
#' group <- rep(LETTERS[1:10], each=nsamples / 10)
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
get_pexpr <- function(data, group, threshold=0, digits=2){

  if(ncol(data)!=length(group)) stop("ncol(data) != length(group)")
  if(!is.numeric(threshold) | threshold < 0) stop("threshold must be numeric and > 0")

  datar <- (data>threshold) * 1
  a <- rowsum(x=t(datar), group=group)
  b <- as.numeric(table(group)[rownames(a)])
  f <- round(100*t(apply(a, 2, function(x) x/b)), digits=digits)
  f

}

#' Makes all pairwise contrasts based on an intercept-less design
#' @importFrom utils combn
#' @importFrom limma makeContrasts
#' @keywords internal
.make_contrasts_pairwise <- function(u, design, delim){

  all_combinations <- combn(u, 2)

  all_contrasts <- lapply(1:ncol(all_combinations), function(cc){

    one <- all_combinations[1,cc,drop=TRUE]
    two <- all_combinations[2,cc,drop=TRUE]

    p_raw <- paste0(one, "-", two)
    p_valid <- paste0(make.names(one), "-", make.names(two))
    col_valid <- make.names(colnames(design))

    co <- suppressWarnings(makeContrasts(contrasts = p_valid, levels=col_valid))
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
.make_contrasts_average <- function(u, design){

  all_contrasts <- lapply(u, function(cc){

    # potentially invalid
    one <- cc
    two <- setdiff(u, cc)
    two2 <- paste0("(", paste(two, collapse="+"), ") / ", length(two))

    # valid
    onev <- make.names(one)
    twov <- make.names(two)
    two2v <- paste0("(", paste(twov, collapse="+"), ") / ", length(twov))
    cnv <- make.names(colnames(design))
    cov <-paste0(onev, " - ", two2v)
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
.rankmatrix <- function(x){
  u = unique(c(x, recursive = TRUE))
  N = length(u)

  rmat = matrix(NA, nrow=length(u), ncol=length(x),
                dimnames=list(u, names(x)))

  for (i in names(x)) {

    rmat[x[[i]], i] = seq(1, length(x[[i]]))

  }

  return(rmat)

}
