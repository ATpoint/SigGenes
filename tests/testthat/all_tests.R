#devtools::document()
#library(SingleCellExperiment)
#library(testthat)
#testthat::test_file(path = "tests/testthat/all_tests.R")

test_that("de_testing works", {
  set.seed(1)
  sce <- .mockSCE(ncells = 100, ngenes = 11, nspikes = 0)
  sce$group <- rep(LETTERS[1:4], each = 25)
  sce$donor <- rep(LETTERS[1:4], 25)

  res1 <- de_limma(x = sce, aggregate_by = c("group", "donor"), main_covariate = "group", other_covariates = c("donor"), return_object = TRUE)
  expect_true(is(res1, "SimpleList"))
  expect_true(is(res1$results, "SimpleList"))
  expect_true(is(res1$results$A_vs_B, "data.frame"))
  expect_true(is(res1$DGEList, "DGEList"))
  expect_true(is(res1$params, "list"))
  rm(res1)

  res2 <- de_limma(x = sce, aggregate_by = c("group", "donor"), main_covariate = "group", other_covariates = c("donor"), return_object = TRUE, mode = "average")
  expect_true(is(res2, "SimpleList"))
  expect_true(is(res2$results, "SimpleList"))
  expect_true(is(res2$results$A, "data.frame"))
  expect_true(is(res2$DGEList, "DGEList"))
  expect_true(is(res2$params, "list"))
  rm(res2)

  expect_error(de_limma(x = sce, main_covariate = "group", use_assay = "not_existing"))
  expect_error(de_limma(x = sce, main_covariate = "group", block = "not_existing"))
  expect_error(de_limma(x = sce, other_covariates = "g"))
  expect_error(de_limma(x = sce, aggregate_by = "g"))
  expect_error(de_limma(x = assay(sce), aggregate_by = "g"))
  expect_error(de_limma(x = assay(sce), min_pct = -1))
  expect_error(de_limma(x = assay(sce), min_fc = -1))
  expect_error(de_limma(x = assay(sce), mode = "3"))
})

# Lot of ChatGPT suggestions
test_that("Fails if input is not a SingleCellExperiment", {
  expect_error(de_limma(x = matrix(1:4, nrow = 2)), "x must be a SingleCellExperiment")
})

test_that("Fails if assay is missing", {
  sce <- SingleCellExperiment(list(counts = matrix(1:4, nrow = 2)))
  expect_error(de_limma(sce, use_assay = "logcounts", main_covariate = "group"), "x contains no assay named")
})

test_that("Fails if aggregate_by contains invalid columns", {
  sce <- SingleCellExperiment(list(counts = matrix(1:4, nrow = 2)))
  colData(sce)$group <- c("A", "B")
  expect_error(de_limma(sce, aggregate_by = "invalid", main_covariate = "group"), "Not all entries in aggregate_by are in colData")
})

test_that("Fails if use_directly is not logical", {
  sce <- SingleCellExperiment(list(counts = matrix(1:4, nrow = 2)))
  colData(sce)$group <- c("A", "B")
  expect_error(de_limma(sce, use_directly = "yes", main_covariate = "group"), "use_directly must be logical")
})

test_that("Fails if main_covariate is not in colData", {
  sce <- SingleCellExperiment(list(counts = matrix(1:4, nrow = 2)))
  expect_error(de_limma(sce, main_covariate = "group"), "main_covariate not in colData")
})

test_that("Fails if main_covariate is not length 1", {
  sce <- SingleCellExperiment(list(counts = matrix(1:4, nrow = 2)))
  colData(sce)$group <- c("A", "B")
  expect_error(de_limma(sce, main_covariate = c("group", "group2")), "main_covariate must be length 1")
})

test_that("Fails if other_covariates not in colData", {
  sce <- SingleCellExperiment(list(counts = matrix(1:4, nrow = 2)))
  colData(sce)$group <- c("A", "B")
  expect_error(de_limma(sce, main_covariate = "group", other_covariates = "not_in_cd"), "other_covariates not in colData")
})

test_that("Fails if block is not in colData", {
  sce <- SingleCellExperiment(list(counts = matrix(1:4, nrow = 2)))
  colData(sce)$group <- c("A", "B")
  expect_error(de_limma(sce, main_covariate = "group", block = "missing"), "Name of block not in colData")
})

test_that("Fails if min_pct is invalid", {
  sce <- SingleCellExperiment(list(counts = matrix(1:4, nrow = 2)))
  colData(sce)$group <- c("A", "B")
  expect_error(de_limma(sce, main_covariate = "group", min_pct = -1), "min_pct must be numeric and between 0 and 100")
})

test_that("Fails if min_fc is too low", {
  sce <- SingleCellExperiment(list(counts = matrix(1:4, nrow = 2)))
  colData(sce)$group <- c("A", "B")
  expect_error(de_limma(sce, main_covariate = "group", min_fc = 0.5), "min_fc must be >= 1")
})

test_that("Fails if use_weights is not logical", {
  sce <- SingleCellExperiment(list(counts = matrix(1:4, nrow = 2)))
  colData(sce)$group <- c("A", "B")
  expect_error(de_limma(sce, main_covariate = "group", use_weights = "TRUE"), "use_weights must be logical")
})

test_that("Fails if both aggregate_by and use_directly are set", {
  sce <- SingleCellExperiment(list(counts = matrix(1:4, nrow = 2)))
  colData(sce)$group <- c("A", "B")
  colData(sce)$sample <- c("S1", "S2")
  expect_error(de_limma(sce, aggregate_by = "sample", use_directly = TRUE, main_covariate = "group"), "aggregate_by and use_directly cannot be used both at the same time")
})

test_that("Fails if no replication in any group", {
  mat <- matrix(rpois(10, lambda = 1), ncol = 5)
  sce <- SingleCellExperiment(list(counts = mat))
  colData(sce)$group <- LETTERS[1:5]
  expect_error(de_limma(sce, main_covariate = "group"))
})

test_that("create_signatures works", {
  set.seed(1)
  sce <- scuttle::mockSCE(ncells = 100, ngenes = 11, nspikes = 0)
  sce$group <- rep(LETTERS[1:4], each = 25)
  sce$donor <- rep(LETTERS[1:4], 25)

  res1 <- de_limma(x = sce, aggregate_by = c("group", "donor"), main_covariate = "group", other_covariates = c("donor"), return_object = TRUE)
  sig <- create_signatures(res1, signif_threshold = .9)
  expect_true(is(sig, "list"))

  res_backup <- res1
  res1$params <- NULL
  expect_error(create_signatures(res1))

  res1 <- res_backup
  res1$results <- NULL
  expect_error(create_signatures(res1))

  res1 <- res_backup
  res1$params$min_pct <- NULL
  expect_error(create_signatures(res1))

  cn <- colnames(res1$results$A_vs_B)
  res1 <- res_backup
  colnames(res1$results$A_vs_B)[which(cn == "logFC")] <- "x"
  expect_error(create_signatures(res1))

  res1 <- res_backup
  colnames(res1$results$A_vs_B)[which(cn == "adj.P.Val")] <- "x"
  expect_error(create_signatures(res1))

  res1 <- res_backup
  colnames(res1$results$A_vs_B)[which(cn == "P.Value")] <- "x"
  expect_error(create_signatures(res1))

  res1 <- res_backup
  colnames(res1$results$A_vs_B)[which(cn == "Gene")] <- "x"
  expect_error(create_signatures(res1))

  expect_error(create_signatures(res_backup, min_prop = 2))
  expect_error(create_signatures(res_backup, n = -1))
  expect_error(create_signatures(res_backup, exclude_groups = "not_existing"))
})
