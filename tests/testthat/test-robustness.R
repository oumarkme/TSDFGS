# Input checks, degenerate cases, and the dual (Gram) path with more columns than rows.
set.seed(3)
wide <- matrix(rnorm(60 * 120), 60, 120)   # p > n for every training set below

test_that("min.iter must be a single finite number of at least 1", {
  g <- wide[, 1:5]
  for (bad in list(NA, Inf, -1, 0, c(10, 20), "a")) {
    expect_error(optTrain(g, 1:50, 10, min.iter = bad, console = FALSE), "min.iter")
  }
})

test_that("n.train, cand, and test are validated", {
  g <- wide[, 1:5]
  expect_error(optTrain(g, 1:50, 10.7, console = FALSE), "whole number")
  expect_error(optTrain(g, 1:50, NA, console = FALSE), "whole number")
  expect_error(optTrain(g, c(1:49, 1), 10, console = FALSE), "duplicated")
  expect_error(optTrain(g, 1:50, 10, test = 51, console = FALSE), "at least 2")
})

test_that("a score that is not finite stops the search with a clear error", {
  g <- wide[1:40, 1:10]
  expect_error(optTrain(g, 1:40, 39, min.iter = 5, console = FALSE), "not finite|at least 2")
  expect_error(optTrain(g, 1:39, 10, method = "CDmean", test = c(40, 40), min.iter = 5, console = FALSE),
               "not finite")
})

test_that("sub-populations with a test set never duplicate a training individual", {
  g <- wide[, 1:6]
  sp <- rep(c("A", "B"), 30)
  res <- optTrain(g, 1:40, 20, subpop = sp, test = 41:60, min.iter = 60, console = FALSE, seed = 5)
  expect_equal(anyDuplicated(res$OPTtrain), 0L)
  expect_error(optTrain(g, 1:40, 20, subpop = c(rep("A", 40), rep("B", 20)), test = 41:60,
                        min.iter = 20, console = FALSE),
               "fewer candidates")
})

test_that("score_many equals score_rows for r-score, PEV, and CD, with and without a Gram matrix", {
  sets <- vapply(1:6, function(i) sample(1:40, 15), integer(15))
  test <- 41:60
  for (id in 0:2) {
    ref <- apply(sets, 2, function(s) score_rows(id, wide, s, test, 0, FALSE))
    expect_equal(score_many(id, wide, sets, test, 0, FALSE, 2L, NULL), ref, tolerance = 1e-10)
    expect_equal(score_many(id, wide, sets, test, 0, FALSE, 1L, tcrossprod(wide)), ref, tolerance = 1e-10)
  }
})

test_that("the dual form (p > n) matches the explicit ridge formula", {
  tr <- 1:30; te <- 31:60
  x <- wide[tr, ]; x0 <- wide[te, ]
  b <- solve(crossprod(x) + diag(1, ncol(x)))
  # CD of the package: ridge on the columns, num = ||x0 (X'X + lam I)^-1 X'||^2 per row
  m <- x0 %*% b %*% t(x)
  expect_equal(cd_score(x, x0), mean(rowSums(m^2) / rowSums(x0^2)) , tolerance = 1e-8)
})

test_that("the v3 criteria run through optTrain with a test set", {
  g <- wide[, 1:8]
  for (m in c("CDmean", "CDmax", "MSPE", "PAM", "AvgGRM")) {
    res <- optTrain(g, 1:40, 12, test = 41:60, method = m, min.iter = 40, console = FALSE, seed = 2)
    expect_length(res$OPTtrain, 12)
    expect_true(all(is.finite(res$TOPscore)), info = m)
  }
})

test_that("optTrain improves over the random start", {
  g <- wide[, 1:8]
  res <- optTrain(g, 1:50, 15, min.iter = 150, console = FALSE, seed = 4)
  expect_gt(tail(res$TOPscore, 1), res$TOPscore[1] + 1e-12)
})

test_that("every criterion runs through all four search routes of optTrain", {
  # Routes: (1) no sub-populations, no test set; (2) test set; (3) sub-populations; (4) both.
  set.seed(21)
  g <- matrix(rnorm(80 * 40), 80, 40)   # 40 columns: GVoverall needs at least n.train columns
  cand <- 1:60; test <- 61:80; n <- 12
  sp_cand <- rep(c("A", "B"), 30)
  sp_all <- rep(c("A", "B"), 40)
  check <- function(fit, spec, label) {
    expect_length(fit$OPTtrain, n)
    expect_equal(anyDuplicated(fit$OPTtrain), 0L, info = label)
    expect_true(all(fit$OPTtrain %in% cand), info = label)
    expect_true(all(is.finite(fit$TOPscore)), info = label)
    d <- diff(fit$TOPscore)   # the best score never gets worse
    expect_true(all(if (spec$maximize) d >= -1e-12 else d <= 1e-12), info = label)
  }
  for (m in criteria()$method) {
    spec <- criteria()[criteria()$method == m, ]
    run <- function(...) optTrain(g, cand, n, method = m, min.iter = 30, console = FALSE, seed = 1, ...)
    check(run(), spec, paste(m, "plain"))
    check(run(subpop = sp_cand), spec, paste(m, "subpop"))
    if (spec$uses_test) {
      check(run(test = test), spec, paste(m, "test"))
      check(run(subpop = sp_all, test = test), spec, paste(m, "subpop+test"))
    } else {
      expect_warning(run(test = test), "test set is ignored")
      expect_warning(run(subpop = sp_all, test = test), "test set is ignored")
    }
  }
})
