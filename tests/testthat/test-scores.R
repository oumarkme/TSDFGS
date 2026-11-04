test_that("scores match the v2.4.2 reference formulas for p < n, p = n, and p > n", {
  set.seed(1)
  for (p in c(5, 25, 40)) {  # n = 25 training rows
    X <- matrix(rnorm(25 * p), 25, p)
    X0 <- matrix(rnorm(40 * p), 40, p)
    expect_equal(r_score(X, X0), ref_r_score(X, X0), tolerance = 1e-9)
    expect_equal(pev_score(X, X0), ref_pev_score(X, X0), tolerance = 1e-9)
    expect_equal(cd_score(X, X0), ref_cd_score(X, X0), tolerance = 1e-9)
  }
})

test_that("scores reproduce the values recorded from v2.4.2", {
  fx <- load_fixture()
  g <- fx$g
  # v2.4.2 returned single-precision floats, so agreement is to about 7 digits.
  tol <- 1e-6
  test_idx <- 101:120
  expect_equal(sapply(fx$sets, function(s) r_score(g[s, ], g[test_idx, ])), fx$scores$r, tolerance = tol)
  expect_equal(sapply(fx$sets, function(s) pev_score(g[s, ], g[test_idx, ])), fx$scores$pev, tolerance = tol)
  expect_equal(sapply(fx$sets, function(s) cd_score(g[s, ], g[test_idx, ])), fx$scores$cd, tolerance = tol)
  expect_equal(sapply(fx$sets, function(s) r_score(g[s, ], g[setdiff(1:100, s), ])), fx$scores$r_un, tolerance = tol)
  expect_equal(sapply(fx$sets, function(s) pev_score(g[s, ], g[setdiff(1:100, s), ])), fx$scores$pev_un, tolerance = tol)
  expect_equal(sapply(fx$sets, function(s) cd_score(g[s, ], g[setdiff(1:100, s), ])), fx$scores$cd_un, tolerance = tol)
})

test_that("scores are returned in double precision and accept integer matrices", {
  set.seed(2)
  X <- matrix(rnorm(60), 12, 5); X0 <- matrix(rnorm(30), 6, 5)
  v <- r_score(X, X0)
  expect_false(v == as_float(v))  # has digits beyond single precision
  Xi <- matrix(sample(-1:1, 60, TRUE), 12, 5); X0i <- matrix(sample(-1:1, 30, TRUE), 6, 5)
  storage.mode(Xi) <- "integer"; storage.mode(X0i) <- "integer"
  expect_equal(pev_score(Xi, X0i), ref_pev_score(Xi + 0, X0i + 0), tolerance = 1e-9)
})

test_that("scores reject inputs with different numbers of columns", {
  expect_error(r_score(matrix(1, 5, 3), matrix(1, 4, 2)), "same number of columns")
  expect_error(cd_score(matrix(1, 5, 3), matrix(1, 4, 2)), "same number of columns")
  expect_error(pev_score(matrix(1, 5, 3), matrix(1, 4, 2)), "same number of columns")
})
