# testthat::test_file("tests/testthat/all_tests.R")

test_that("de_testing works", {

  set.seed(1)
  sce <- scuttle::mockSCE(ncells=100, ngenes=11, nspikes=0)
  sce$group <- rep(LETTERS[1:4], each=25)
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

test_that("create_signatures works", {

  set.seed(1)
  sce <- scuttle::mockSCE(ncells=100, ngenes=11, nspikes=0)
  sce$group <- rep(LETTERS[1:4], each=25)
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
  colnames(res1$results$A_vs_B)[which(cn=="logFC")] <- "x"
  expect_error(create_signatures(res1))

  res1 <- res_backup
  colnames(res1$results$A_vs_B)[which(cn=="adj.P.Val")] <- "x"
  expect_error(create_signatures(res1))

  res1 <- res_backup
  colnames(res1$results$A_vs_B)[which(cn=="P.Value")] <- "x"
  expect_error(create_signatures(res1))

  res1 <- res_backup
  colnames(res1$results$A_vs_B)[which(cn=="Gene")] <- "x"
  expect_error(create_signatures(res1))

  expect_error(create_signatures(res_backup, min_prop = 2))
  expect_error(create_signatures(res_backup, n = -1))
  expect_error(create_signatures(res_backup, exclude_groups = "not_existing"))

})


