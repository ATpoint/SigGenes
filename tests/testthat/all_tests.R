library(testthat)
library(SingleCellExperiment)

make_test_sce <- function(ncells = 120, ngenes = 100) {
  set.seed(1)

  sce <- mockSCE(ncells = ncells, ngenes = ngenes, nspikes = 0)

  sce$group <- rep(LETTERS[1:3], each = ncells / 3)
  sce$donor <- rep(rep(paste0("D", 1:4), each = ncells / 12), 3)

  sce
}

test_that("pairwise pseudobulk analysis returns expected structure", {
  sce <- make_test_sce()

  res <- de_limma(
    x = sce,
    aggregate_by = c("group", "donor"),
    main_covariate = "group",
    other_covariates = "donor",
    mode = "pairwise",
    prefilter_method = "none",
    return_object = TRUE
  )

  expect_s4_class(res, "SimpleList")
  expect_s4_class(res$results, "SimpleList")
  expect_named(res$results, c("A_vs_B", "A_vs_C", "B_vs_C"))
  expect_true(all(vapply(res$results, is.data.frame, logical(1))))
  expect_s4_class(res$sce, "SingleCellExperiment")
  expect_type(res$params, "list")

  expect_true(all(c("Gene", "logFC", "P.Value", "adj.P.Val") %in% colnames(res$results$A_vs_B)))
  expect_identical(res$params$mode, "pairwise")
  expect_identical(res$params$limma_method, "trend")
  expect_identical(res$params$filtering_method, "none")
})

test_that("average pseudobulk analysis returns one result per group", {
  sce <- make_test_sce()

  res <- de_limma(
    x = sce,
    aggregate_by = c("group", "donor"),
    main_covariate = "group",
    other_covariates = "donor",
    mode = "average",
    prefilter_method = "none"
  )

  expect_s4_class(res, "SimpleList")
  expect_s4_class(res$results, "SimpleList")
  expect_named(res$results, c("A", "B", "C"))
  expect_true(all(vapply(res$results, is.data.frame, logical(1))))
  expect_null(res$sce)
})

test_that("ordinary limma moderation runs", {
  sce <- make_test_sce()

  res <- de_limma(
    x = sce,
    aggregate_by = c("group", "donor"),
    main_covariate = "group",
    other_covariates = "donor",
    limma_method = "ordinary",
    prefilter_method = "none"
  )

  expect_identical(res$params$limma_method, "ordinary")
  expect_true(length(res$results) > 0)
})

test_that("percent-expressed filtering adds percentage columns", {
  sce <- make_test_sce()

  res <- de_limma(
    x = sce,
    aggregate_by = c("group", "donor"),
    main_covariate = "group",
    other_covariates = "donor",
    prefilter_method = "percent_expressed",
    min_pct = 0
  )

  expect_true(all(c("pct.1", "pct.2") %in% colnames(res$results$A_vs_B)))
})

test_that("filterByExpr filtering runs", {
  sce <- make_test_sce()

  res <- de_limma(
    x = sce,
    aggregate_by = c("group", "donor"),
    main_covariate = "group",
    other_covariates = "donor",
    prefilter_method = "filterByExpr"
  )

  expect_s4_class(res, "SimpleList")
  expect_true(length(res$results) == 3)
})

test_that("use_directly skips filtering and normalization", {
  sce <- make_test_sce()

  logcounts <- edgeR::cpm(counts(sce), log = TRUE)
  assay(sce, "logcounts") <- logcounts

  res <- de_limma(
    x = sce,
    use_assay = "logcounts",
    use_directly = TRUE,
    main_covariate = "group",
    prefilter_method = "filterByExpr"
  )

  expect_identical(res$params$filtering_method, "none")
  expect_true(length(res$results) == 3)
})

test_that("sample weights run", {
  sce <- make_test_sce()

  res <- de_limma(
    x = sce,
    aggregate_by = c("group", "donor"),
    main_covariate = "group",
    other_covariates = "donor",
    prefilter_method = "none",
    use_weights = TRUE
  )

  expect_s4_class(res, "SimpleList")
  expect_true(length(res$results) == 3)
})

test_that("x must be a SingleCellExperiment", {
  expect_error(
    de_limma(matrix(1:20, nrow = 5), main_covariate = "group"),
    "x must be a SingleCellExperiment"
  )
})

test_that("use_assay must identify an existing assay", {
  sce <- make_test_sce()

  expect_error(
    de_limma(sce, use_assay = "missing", main_covariate = "group"),
    "use_assay must name an assay in x"
  )
})

test_that("aggregate_by must contain valid colData columns", {
  sce <- make_test_sce()

  expect_error(
    de_limma(sce, aggregate_by = "missing", main_covariate = "group"),
    "aggregate_by must contain column names from colData"
  )
})

test_that("aggregate_by must include all model covariates", {
  sce <- make_test_sce()

  expect_error(
    de_limma(
      sce,
      aggregate_by = "donor",
      main_covariate = "group",
      other_covariates = "donor"
    ),
    "Not all covariates are present in metadata"
  )

  expect_error(
    de_limma(
      sce,
      aggregate_by = "group",
      main_covariate = "group",
      other_covariates = "donor"
    ),
    "Not all covariates are present in metadata"
  )
})

test_that("main_covariate must be one valid colData column", {
  sce <- make_test_sce()

  expect_error(
    de_limma(sce, main_covariate = "missing"),
    "main_covariate must name one column in colData"
  )

  expect_error(
    de_limma(sce, main_covariate = c("group", "donor")),
    "main_covariate must name one column in colData"
  )

  expect_error(
    de_limma(sce, main_covariate = character()),
    "main_covariate must name one column in colData"
  )
})

test_that("other_covariates must be valid colData columns", {
  sce <- make_test_sce()

  expect_error(
    de_limma(sce, main_covariate = "group", other_covariates = "missing"),
    "other_covariates must contain column names from colData"
  )
})

test_that("numeric thresholds are validated", {
  sce <- make_test_sce()

  expect_error(
    de_limma(sce, main_covariate = "group", min_pct = -1),
    "min_pct must be a numeric scalar greater than or equal to zero"
  )

  expect_error(
    de_limma(sce, main_covariate = "group", min_pct = NA_real_),
    "min_pct must be a numeric scalar greater than or equal to zero"
  )

  expect_error(
    de_limma(sce, main_covariate = "group", min_pct = c(0, 10)),
    "min_pct must be a numeric scalar greater than or equal to zero"
  )

  expect_error(
    de_limma(sce, main_covariate = "group", min_expr = -1),
    "min_expr must be a numeric scalar greater than or equal to zero"
  )

  expect_error(
    de_limma(sce, main_covariate = "group", min_fc = 0.5),
    "min_fc must be a numeric scalar greater than or equal to one"
  )

  expect_error(
    de_limma(sce, main_covariate = "group", min_fc = NA_real_),
    "min_fc must be a numeric scalar greater than or equal to one"
  )
})

test_that("logical arguments are validated", {
  sce <- make_test_sce()

  expect_error(
    de_limma(sce, main_covariate = "group", use_directly = "yes"),
    "use_directly must be TRUE or FALSE"
  )

  expect_error(
    de_limma(sce, main_covariate = "group", use_weights = "yes"),
    "use_weights must be TRUE or FALSE"
  )

  expect_error(
    de_limma(sce, main_covariate = "group", return_object = 1),
    "return_object must be TRUE or FALSE"
  )

  expect_error(
    de_limma(sce, main_covariate = "group", verbose = NA),
    "verbose must be TRUE or FALSE"
  )
})

test_that("delim must be a non-empty character scalar", {
  sce <- make_test_sce()

  expect_error(
    de_limma(sce, main_covariate = "group", delim = ""),
    "delim must be a non-empty character scalar"
  )

  expect_error(
    de_limma(sce, main_covariate = "group", delim = 1),
    "delim must be a non-empty character scalar"
  )
})

test_that("enumerated arguments are validated", {
  sce <- make_test_sce()

  expect_error(
    de_limma(sce, main_covariate = "group", mode = "invalid")
  )

  expect_error(
    de_limma(sce, main_covariate = "group", limma_method = "voom")
  )

  expect_error(
    de_limma(sce, main_covariate = "group", prefilter_method = "invalid")
  )
})

test_that("pairwise delimiter is reflected in result names", {
  sce <- make_test_sce()

  res <- de_limma(
    x = sce,
    aggregate_by = c("group", "donor"),
    main_covariate = "group",
    other_covariates = "donor",
    prefilter_method = "none",
    delim = "__against__"
  )

  expect_named(
    res$results,
    c("A__against__B", "A__against__C", "B__against__C")
  )
})

test_that("create_signatures returns a list", {
  sce <- make_test_sce()

  res <- de_limma(
    x = sce,
    aggregate_by = c("group", "donor"),
    main_covariate = "group",
    other_covariates = "donor",
    prefilter_method = "none"
  )

  sig <- create_signatures(res, signif_threshold = 0.9)

  expect_type(sig, "list")
})

test_that("create_signatures requires expected result structure", {
  sce <- make_test_sce()

  res <- de_limma(
    x = sce,
    aggregate_by = c("group", "donor"),
    main_covariate = "group",
    other_covariates = "donor",
    prefilter_method = "none"
  )

  bad <- res
  bad$params <- NULL
  expect_error(create_signatures(bad))

  bad <- res
  bad$results <- NULL
  expect_error(create_signatures(bad))

  bad <- res
  bad$params$min_pct <- NULL
  expect_error(create_signatures(bad))
})

test_that("create_signatures requires expected result columns", {
  sce <- make_test_sce()

  res <- de_limma(
    x = sce,
    aggregate_by = c("group", "donor"),
    main_covariate = "group",
    other_covariates = "donor",
    prefilter_method = "none"
  )

  for (column in c("logFC", "adj.P.Val", "P.Value", "Gene")) {
    bad <- res
    names(bad$results[[1]])[names(bad$results[[1]]) == column] <- "missing_column"
    expect_error(create_signatures(bad))
  }
})

test_that("create_signatures validates arguments", {
  sce <- make_test_sce()

  res <- de_limma(
    x = sce,
    aggregate_by = c("group", "donor"),
    main_covariate = "group",
    other_covariates = "donor",
    prefilter_method = "none"
  )

  expect_error(create_signatures(res, min_prop = 2))
  expect_error(create_signatures(res, n = -1))
  expect_error(create_signatures(res, exclude_groups = "not_existing"))
})

