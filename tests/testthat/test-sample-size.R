test_that("nt2r returns one r-score per repeat, and larger training sets score higher", {
  fx <- load_fixture()
  set.seed(1)
  r10 <- nt2r(fx$g, 10, n_iter = 3)
  r40 <- nt2r(fx$g, 40, n_iter = 3)
  expect_length(r10, 3)
  expect_true(all(is.finite(r10)))
  expect_gt(mean(r40), mean(r10))
})

test_that("FGCM fits a logistic curve and SSDFGS returns the v2 result structure", {
  skip_on_cran()
  fx <- load_fixture()
  set.seed(2)
  nt <- c(10, 15, 20, 30, 40, 50)
  fit <- FGCM(fx$g, nt = nt, n_iter = 2)
  expect_named(fit, c("alpha", "beta", "gamma", "sim"))
  expect_equal(nrow(fit$sim), length(nt) * 2)
  expect_true(fit$gamma > 0)
  out <- capture.output(res <- SSDFGS(fx$g, nt = nt, n_iter = 2))
  expect_named(res, c("OC.fig", "GC.fig", "parameter", "OC.fit"))
  expect_named(res$parameter, c("alpha", "beta", "gamma"))
  expect_s3_class(res$OC.fit, "data.frame")
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    expect_s3_class(res$OC.fig, "ggplot")
    expect_s3_class(res$GC.fig, "ggplot")
  } else {
    expect_null(res$OC.fig)
  }
})

test_that("v2 entry points and arguments are still available", {
  expect_true(all(c("optTrain", "r_score", "pev_score", "cd_score", "SSDFGS", "FGCM", "nt2r") %in% getNamespaceExports("TSDFGS")))
  expect_named(formals(optTrain), c("geno", "cand", "n.train", "subpop", "test", "method", "min.iter", "console", "seed"))
  expect_named(formals(SSDFGS), c("geno", "nt", "n_iter", "multi.threads"))
  expect_named(formals(FGCM), c("geno", "nt", "n_iter", "multi.threads"))
  expect_named(formals(nt2r), c("geno", "nt", "n_iter", "multi.threads"))
  expect_named(formals(r_score), c("X", "X0"))
})
