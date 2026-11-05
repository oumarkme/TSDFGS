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
  # v3 only appends arguments (lambda), so positional v2 calls keep working
  expect_equal(names(formals(optTrain))[1:9], c("geno", "cand", "n.train", "subpop", "test", "method", "min.iter", "console", "seed"))
  expect_equal(names(formals(SSDFGS))[1:4], c("geno", "nt", "n_iter", "multi.threads"))
  expect_equal(names(formals(FGCM))[1:4], c("geno", "nt", "n_iter", "multi.threads"))
  expect_equal(names(formals(nt2r))[1:4], c("geno", "nt", "n_iter", "multi.threads"))
  expect_equal(names(formals(r_score))[1:2], c("X", "X0"))
})

test_that("the parallel simulation does not depend on the number of threads", {
  fx <- load_fixture()
  set.seed(7); a <- nt2r(fx$g, 20, n_iter = 6, threads = 1)
  set.seed(7); b <- nt2r(fx$g, 20, n_iter = 6, threads = 3)
  expect_identical(a, b)
  expect_length(a, 6)
  expect_true(all(is.finite(a) & a > 0))
})

test_that("the parallel simulation gives r-scores like the serial run", {
  skip_on_cran()
  fx <- load_fixture()
  set.seed(3); serial <- nt2r(fx$g, 20, n_iter = 12)
  set.seed(3); par <- nt2r(fx$g, 20, n_iter = 12, threads = 2)
  expect_lt(abs(mean(serial) - mean(par)), 3 * sd(serial) / sqrt(12) + 0.01)
})

test_that("a bad threads argument is an error", {
  fx <- load_fixture()
  expect_error(nt2r(fx$g, 20, n_iter = 2, threads = 0), "threads")
})
