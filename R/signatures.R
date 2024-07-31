#' Internal function to filter and rank DEG results
#'
#' @author Alexander Bender
#'
#' @param res The output of \code{de_limma} run with \code{mode="pairwise"}
#'
#' @inheritParams create_signatures
#'
#' @importFrom methods as is
#'
#' @keywords internal
#'
rank_degs <- function(res, delim,
                      signif_column, signif_threshold,
                      effect_column, effect_threshold,
                      rank_column, rank_method,
                      gene_column, min_pct){

  if(!class(res) %in% c("list", "SimpleList") | is.null(names(res))){
    stop("res must be a named list", call.=FALSE)
  }

  if(!all(grepl(delim, names(res)))) stop("delim was not found in all names of the res")

  invisible(match.arg(arg=class(signif_threshold), choices="numeric"))
  invisible(match.arg(arg=class(effect_threshold), choices="numeric"))

  check.gf <- sum(unlist(lapply(res, function(x)
    if(!gene_column %in% colnames(x)) return(1) else return(0))))
  check.sf <- sum(unlist(lapply(res, function(x)
    if(!signif_column %in% colnames(x)) return(1) else return(0))))
  check.ef <- sum(unlist(lapply(res, function(x)
    if(!effect_column %in% colnames(x)) return(1) else return(0))))
  check.rk <- sum(unlist(lapply(res, function(x)
    if(!rank_column %in% colnames(x)) return(1) else return(0))))

  if(check.gf>0) stop("gene_column does not exist in all entries of res.list!")
  if(check.sf>0) stop("signif_column does not exist in all entries of res.list!")
  if(check.ef>0) stop("effect_column does not exist in all entries of res.list!")
  if(check.rk>0) stop("rank_column does not exist in all entries of res.list!")

  rank_method <- match.arg(rank_method, c("decreasing", "increasing"))

  # Ranking
  nm <- data.frame(t(as.data.frame(strsplit(names(res), delim), check.names = FALSE)), check.names = FALSE)
  colnames(nm) <- c("first", "second")
  rownames(nm) <- NULL
  nm2 <- nm[,2:1]
  nm$direction <- "forward"
  nm2$direction <- "reverse"
  colnames(nm2) <- colnames(nm)
  nm <- rbind(nm, nm2)

  u <- unique(c(nm$first, nm$second))

  l <- sapply(u, function(f){

    y <- nm[nm$first %in% f,]

    y2 <- lapply(1:nrow(y), function(z){

      yy <- y[z,,drop=FALSE]

      if(yy$direction=="forward"){

        current_contrast <- paste0(yy$first, delim, yy$second)

        here <- res[[current_contrast]]
        here <- here[here$pct.1 >= min_pct,]

      } else {

        current_contrast <- paste0(yy$second, delim, yy$first)
        here <- res[[current_contrast]]
        here[[effect_column]] <- here[[effect_column]] * -1
        here <- here[here$pct.2 >= min_pct,]

      }

      h <- here[here[[signif_column]] < signif_threshold & here[[effect_column]] > effect_threshold,]
      k <- if(rank_method=="decreasing") TRUE else FALSE
      o <- order(h[[rank_column]], decreasing = k)
      h[o,][[gene_column]]

    })

    names(y2) <- y$second
    return(y2)


  }, simplify = FALSE)

  l <- as(l, "SimpleList")
  return(l)

}

#' Convert ranks to signatures
#'
#' @param ranked The output of \code{rank_degs}
#' @author Alexander Bender
#' @keywords internal
#' @inheritParams create_signatures
#'
ranks2signatures <- function(ranked, min_prop, n, exclude_groups){

  if(!is.null(exclude_groups)){
    if(!sum(exclude_groups %in% names(ranked))==length(exclude_groups))
      stop("Specified exclude_groups does not exist as a group name")
  }

  if(is.null(names(ranked))) stop("ranked has no names")
  if(min_prop > 1 | min_prop <= 0) stop("min_prop must be between >= 0 and <= 1")
  if(is.finite(n)){
    if(n <= 0 | (n %% 1 > 0)) stop("n must be an even number larger 0 or Inf")
  }

  # optionally remove exclude_groups from ranked
  sdf <- setdiff(names(ranked), exclude_groups)
  ranked <- ranked[sdf]
  ranked <- sapply(names(ranked), function(x){

    r <- ranked[[x]]
    k <- intersect(names(r), sdf)
    r[k]

  }, simplify = FALSE)

  nms <- names(ranked)
  s <- sapply(nms, function(x){

    genes_ranked <- ranked[[x]]
    genes_ranked[[x]] <- NULL

    # rank matrix, low values means high ranks and viceversa, 1 means not present in that pairwise list:
    rmat <- .rankmatrix(genes_ranked)

    # Get the genes that qualify as markers given the number of groups and the min_prop
    isna <- !is.na(rmat)
    from <- ncol(rmat)
    to   <- floor(min_prop*from)
    to <- if(to==0) 1 else to

    med  <- function(y) median(y,na.rm=TRUE)

    # take the rankmatrix (so each entry is the rank of the gene in the individual ranking),
    # and then calculate the median of the ranks -- that is the final ranking metric
    final <- lapply(from:to, function(n){

      idx <- base::rowSums(isna)==n
      if(sum(idx)==0) return(NULL)

      rmat_x <- rmat[idx,,drop=FALSE]

      # get the median of the ranks
      mat <- apply(rmat_x, 1, med)

      # rank by the median
      o <- order(mat)

      markers.df <- data.frame(rmat_x[o,,drop=FALSE], check.names=FALSE)

      return(markers.df)

    })

    final <- do.call(rbind, final)

    final[!is.na(final)] <- 1
    final[is.na(final)]  <- 0
    final <- head(final, n=n)
    if(length(final) == 0){
      final <- NULL
    }

    return(final)

  },simplify = FALSE)

  # Now, based on min_prop filter for groupwise signature genes
  sigs <- lapply(names(s), function(x){

    sx <- s[[x]]
    if(is.null(sx)) return(NULL)

    min_cols <- floor(min_prop * ncol(sx))

    keep <- apply(sx == 1, 1, sum) >= min_cols
    out <- rownames(sx[keep,,drop=FALSE])
    out

  })

  names(sigs) <- names(s)

  return(sigs)

}

#' Create group-specific gene signatures based on pairwise differential analysis results
#'
#' @param x List with DE results from \code{de_limma} run with \code{mode="pairwise"}
#' @param n The maximum number of top signature genes per group to return
#' @param min_prop A number greater zero and smaller or equal than 1. Defines the stringency of the signature. It is the proportion
#'   of groups that a gene needs to be overexpressed to. A value of 1 means a gene must be overexpressed versus all other groups. A value
#'   of 0.75 means overexpressed versus 75\% of other groups.
#' @param exclude_groups Vector with group names to exclude from signature creation. Can be helpful if a group is expected to be transitional between two other groups
#'   so one would not expect unique genes defining this group. Can help returning more genes for the remaining groups.
#' @param signif_column Name of the significance column in \code{x}.
#' @param signif_threshold Threshold for significance.
#' @param effect_column Name of the effect column in \code{x}.
#' @param effect_threshold Threshold for effect size \(below this is "underexpressed", above is "overexpressed"\)
#' @param gene_column Name of the gene column in \code{x}.
#' @param rank_column Name of the column in \code{x} to use for ranking.
#' @param rank_method Either "decreasing" or "increasing" ranking direction
#'
#' @author Alexander Bender
#'
#' @examples
#' set.seed(1)
#' sce <- scuttle::mockSCE(ncells=1000, ngenes=1000, nspikes=0)
#' sce$group <- rep(LETTERS[1:4], each=250)
#' sce$donor <- rep(LETTERS[1:4], 250)
#' res_pairwise <- de_limma(x = sce, aggregate_by = c("group", "donor"),
#'                          main_covariate = "group", other_covariates = c("donor"))
#' create_signatures <- create_signatures(x = res_pairwise, signif_threshold = .9)
#' create_signatures
#'
#' @export
create_signatures <- function(x, n=Inf, min_prop=1, exclude_groups=NULL,
                              signif_column="adj.P.Val", signif_threshold=0.05,
                              effect_column="logFC", effect_threshold=0,
                              rank_column="P.Value", rank_method="increasing",
                              gene_column="Gene"){

  mode_found <- x$params$mode
  delim_found <- x$params$delim
  pct_found <- x$params$min_pct
  results_found <- x$results

  if(is.null(x$params)){
    stop("This does not look like the (unmodified) output of de_limma!")
  }

  if(is.null(mode_found) | is.null(delim_found) | is.null(pct_found) | is.null(results_found)){
    stop("This does not look like the (unmodified) output of de_limma!")
  }

  ranked <- rank_degs(res=x$results, delim=delim_found, signif_column=signif_column, signif_threshold=signif_threshold,
                      effect_column=effect_column, effect_threshold=effect_threshold,
                      rank_column=rank_column, rank_method=rank_method,
                      gene_column=gene_column, min_pct=x$params$min_pct)

  signatures <- ranks2signatures(ranked=ranked, n=n, min_prop=min_prop, exclude_groups=exclude_groups)

  return(signatures)

}
