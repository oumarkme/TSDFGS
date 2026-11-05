# Regression tests against v2.4.2. Scores are rounded to single precision (the option below)
# so that accept/reject decisions equal those of the float scorers of v2.4.2.
legacy_args <- function(fx, pop, target, method) {
  args <- list(geno = fx$g, cand = 1:100, n.train = 20, method = method,
               min.iter = 230, console = FALSE)
  if (pop) args$subpop <- if (target) fx$sp else fx$sp[1:100]
  if (target) args$test <- 101:120
  args
}

test_that("all twelve routes return the training sets recorded from v2.4.2", {
  skip_on_cran()
  fx <- load_fixture()
  withr_opt <- options(TSDFGS.legacy_float = TRUE)
  on.exit(options(withr_opt), add = TRUE)
  for (pop in c(FALSE, TRUE)) for (target in c(FALSE, TRUE)) for (m in c("rScore", "PEV", "CD")) {
    key <- paste(if (pop) "pop" else "npop", if (target) "target" else "untarget", m, sep = "_")
    set.seed(7)
    res <- do.call(optTrain, legacy_args(fx, pop, target, m))
    old <- fx$runs[[key]]
    expect_equal(res$OPTtrain, old$OPTtrain, info = key)
    expect_equal(res$TOPscore, old$TOPscore, info = key)
    expect_equal(res$ITERscore, old$ITERscore, info = key)
  }
})

test_that("seed gives identical results and a different seed does not", {
  fx <- load_fixture()
  a <- optTrain(fx$g, 1:100, 20, min.iter = 230, console = FALSE, seed = 11)
  b <- optTrain(fx$g, 1:100, 20, min.iter = 230, console = FALSE, seed = 11)
  d <- optTrain(fx$g, 1:100, 20, min.iter = 230, console = FALSE, seed = 12)
  expect_identical(a, b)
  expect_false(identical(a$OPTtrain, d$OPTtrain))
})

test_that("a small min.iter no longer errors (v2.4.2 failed with a negative look-back index)", {
  fx <- load_fixture()
  res <- optTrain(fx$g, 1:100, 20, min.iter = 5, console = FALSE, seed = 1)
  expect_length(res$OPTtrain, 20)
  expect_true(res$n.iter >= 5)
})

test_that("candidates that are not 1:nrow(geno) are respected with sub-populations", {
  fx <- load_fixture()
  set.seed(3)
  cand <- sort(sample(1:120, 80))
  sp <- fx$sp[cand]
  res <- optTrain(fx$g, cand, 16, subpop = sp, min.iter = 150, console = FALSE, seed = 4)
  expect_true(all(res$OPTtrain %in% cand))
  expect_equal(anyDuplicated(res$OPTtrain), 0)
  expected <- ceiling(16 * table(sp) / sum(table(sp)))
  while (sum(expected) > 16) { i <- which.max(expected); expected[i] <- expected[i] - 1 }
  expect_equal(as.vector(table(factor(fx$sp[res$OPTtrain], levels = names(expected)))), as.vector(expected))
})

test_that("a sub-population with one candidate does not break the search", {
  fx <- load_fixture()
  sp <- fx$sp[1:100]
  sp[5] <- "SINGLE"
  res <- optTrain(fx$g, 1:100, 20, subpop = sp, min.iter = 230, console = FALSE, seed = 5)
  expect_length(res$OPTtrain, 20)
  expect_true(all(res$OPTtrain %in% 1:100))
})

test_that("untargeted search improves the criterion over a random start", {
  fx <- load_fixture()
  for (m in c("rScore", "PEV", "CD")) {
    res <- optTrain(fx$g, 1:100, 20, method = m, min.iter = 230, console = FALSE, seed = 6)
    top <- res$TOPscore
    if (m == "PEV") expect_lt(tail(top, 1), top[1]) else expect_gt(tail(top, 1), top[1])
  }
})

test_that("input errors are clear", {
  fx <- load_fixture()
  expect_error(optTrain(fx$g, 1:100, 20, method = "XX", console = FALSE), "Method not found")
  expect_error(optTrain(fx$g, 1:100, 200, console = FALSE), "larger than the number of candidates")
  expect_error(optTrain(fx$g, c(1:100, 500), 20, console = FALSE), "row numbers of geno")
  expect_error(optTrain(fx$g, 1:100, 20, subpop = fx$sp[1:50], console = FALSE), "Input data is not correct")
  g <- fx$g; g[1, 1] <- NA
  expect_error(optTrain(g, 1:100, 20, console = FALSE), "missing or infinite values")
})
