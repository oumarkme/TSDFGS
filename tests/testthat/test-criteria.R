# Dense R references for the criteria added in v3, written from the published definitions.

set.seed(11)
pcs <- matrix(rnorm(60 * 6), 60, 6)
tr <- pcs[1:25, ]
te <- pcs[26:60, ]
wide_tr <- matrix(rnorm(20 * 40), 20, 40)   # more columns than rows
wide_te <- matrix(rnorm(15 * 40), 15, 40)

ridge_inv <- function(x, lambda) solve(crossprod(x) + diag(lambda, ncol(x)))
with_int <- function(x) cbind(1, x)

ref_aopt <- function(x, lambda = 1 / ncol(x)) sum(diag(ridge_inv(x, lambda)))
ref_dopt <- function(x, lambda = 1 / ncol(x)) as.numeric(determinant(crossprod(x) + diag(lambda, ncol(x)))$modulus)
ref_eopt <- function(x, lambda = 1 / ncol(x)) max(eigen(ridge_inv(x, lambda), symmetric = TRUE)$values)
ref_pev_matrix <- function(x, x0, lambda = 1 / ncol(x)) {
  xe <- with_int(x); x0e <- with_int(x0)
  x0e %*% solve(crossprod(xe) + diag(lambda, ncol(xe))) %*% t(x0e)
}
ref_pevmax <- function(x, x0, lambda = 1 / ncol(x)) max(diag(ref_pev_matrix(x, x0, lambda)))
ref_gopt <- function(x, x0, lambda = 1 / ncol(x)) max(eigen(ref_pev_matrix(x, x0, lambda), symmetric = TRUE)$values)
ref_mspe <- function(x, x0, lambda = 1) {
  a <- t(x) %*% solve(x %*% t(x) + diag(lambda, nrow(x)))
  1 + sum(diag(x0 %*% a %*% t(a) %*% t(x0))) / nrow(x0) + sum((x0 - x0 %*% a %*% x)^2)
}
ref_gvo <- function(x) as.numeric(determinant(tcrossprod(x) / ncol(x))$modulus)
ref_gva <- function(x) sum(diag(tcrossprod(x) / ncol(x)))
ref_grm <- function(x, y) mean(x %*% t(y)) / ncol(x)
ref_cd <- function(x, x0, lambda = 1, v2 = FALSE) {
  all <- rbind(x, x0); nt <- nrow(x); np <- nrow(x0); p <- ncol(x)
  k <- tcrossprod(all) / p
  tt <- seq_len(nt); pp <- nt + seq_len(np)
  m <- diag(nt) - matrix(1 / nt, nt, nt)
  w <- solve(m %*% k[tt, tt] + diag(lambda, nt)) %*% m
  cov_all <- k[, tt] %*% w %*% k[tt, ]
  if (v2) return(mean(diag(cov_all)[pp] / diag(k)[pp]))
  cd <- vapply(seq_len(np), function(i) {
    cc <- numeric(nt + np); cc[pp] <- -1 / np; cc[nt + i] <- cc[nt + i] + 1
    drop(t(cc) %*% cov_all %*% cc / (t(cc) %*% k %*% cc))
  }, numeric(1))
  mean(cd)
}
ref_pam <- function(x, x0) {
  d <- sqrt(pmax(outer(rowSums(x0^2), rowSums(x^2), "+") - 2 * x0 %*% t(x), 0))
  sum(apply(d, 1, min))
}

test_that("optimal-design criteria equal the dense formulas for p < n and p > n", {
  for (m in list(tr, wide_tr)) {
    expect_equal(aopt_score(m), ref_aopt(m), tolerance = 1e-9)
    expect_equal(dopt_score(m), ref_dopt(m), tolerance = 1e-9)
  }
  expect_equal(eopt_score(tr), ref_eopt(tr), tolerance = 1e-9)
  expect_equal(eopt_score(wide_tr), 1 / (1 / ncol(wide_tr)), tolerance = 1e-12)
  expect_equal(aopt_score(tr, lambda = 0.5), ref_aopt(tr, 0.5), tolerance = 1e-9)
})

test_that("PEVmax and G-optimality equal the dense formulas for p < n and p > n", {
  expect_equal(pevmax_score(tr, te), ref_pevmax(tr, te), tolerance = 1e-9)
  expect_equal(gopt_score(tr, te), ref_gopt(tr, te), tolerance = 1e-9)
  expect_equal(pevmax_score(wide_tr, wide_te), ref_pevmax(wide_tr, wide_te), tolerance = 1e-9)
  expect_equal(gopt_score(wide_tr, wide_te), ref_gopt(wide_tr, wide_te), tolerance = 1e-9)
  expect_equal(pevmax_score(tr, te, lambda = 2), ref_pevmax(tr, te, 2), tolerance = 1e-9)
})

test_that("MSPE of ridge regression equals the dense formula for p < n and p > n", {
  expect_equal(mspe_score(tr, te), ref_mspe(tr, te), tolerance = 1e-9)
  expect_equal(mspe_score(wide_tr, wide_te), ref_mspe(wide_tr, wide_te), tolerance = 1e-9)
})

test_that("genomic variation and relationship criteria equal the dense formulas", {
  expect_equal(gv_overall_score(wide_tr), ref_gvo(wide_tr), tolerance = 1e-9)
  expect_equal(gv_average_score(tr), ref_gva(tr), tolerance = 1e-12)
  expect_equal(avg_grm_self_score(tr), ref_grm(tr, tr), tolerance = 1e-12)
  expect_equal(avg_grm_score(tr, te), ref_grm(tr, te), tolerance = 1e-12)
  expect_equal(avg_grm_minmax_score(tr, te), ref_grm(tr, te) - ref_grm(tr, tr), tolerance = 1e-12)
  expect_error(gv_overall_score(tr), "at least as many columns")  # fewer columns than rows: singular
})

test_that("CDmean and CDmean(v2) equal the dense GBLUP formulas", {
  expect_equal(cdmean_score(wide_tr, wide_te), ref_cd(wide_tr, wide_te), tolerance = 1e-8)
  expect_equal(cdmean2_score(wide_tr, wide_te), ref_cd(wide_tr, wide_te, v2 = TRUE), tolerance = 1e-8)
  expect_equal(cdmean_score(tr, te, lambda = 0.3), ref_cd(tr, te, 0.3), tolerance = 1e-8)
  cd <- cdmean_score(wide_tr, wide_te)
  expect_true(cd > 0 && cd < 1)
})

test_that("the PAM cost equals the dense k-medoids cost", {
  expect_equal(pam_score(tr, te), ref_pam(tr, te), tolerance = 1e-9)
})

test_that("the CDmax score is the largest per-individual CD of cd_score", {
  b <- solve(crossprod(tr) + diag(ncol(tr)))
  cd_i <- rowSums((te %*% (b %*% crossprod(tr) %*% b)) * te) / rowSums(te^2)
  expect_equal(cdmax_score(tr, te), max(cd_i), tolerance = 1e-9)
  expect_equal(mean(cd_i), cd_score(tr, te), tolerance = 1e-6)
})

test_that("criteria() lists every method that optTrain accepts", {
  cr <- criteria()
  expect_equal(nrow(cr), 18L)
  expect_true(all(vapply(cr$score_fun, function(f) is.function(getExportedValue("TSDFGS", f)), logical(1))))
  expect_error(optTrain(pcs, 1:60, 10, method = "nope"), "Method not found")
})

test_that("optTrain improves every criterion over random training sets", {
  wide <- matrix(rnorm(60 * 40), 60, 40)
  cand <- 1:60
  for (m in setdiff(criteria()$method, c("rScore", "PEV", "CD"))) {
    g <- if (m == "GVoverall") wide else pcs
    n <- 12
    fit <- optTrain(g, cand, n, method = m, min.iter = 150, console = FALSE, seed = 3)
    expect_length(unique(fit$OPTtrain), n)
    spec <- criteria()[criteria()$method == m, ]
    f <- getExportedValue("TSDFGS", spec$score_fun)
    sc <- function(s) {
      rest <- setdiff(cand, s)
      if (!spec$uses_test) f(g[s, ]) else f(g[s, ], if (m %in% c("AvgGRM", "MSPE", "CDmean", "CDmeanV2")) g[cand, ] else g[rest, ])
    }
    set.seed(99)
    rnd <- vapply(1:30, function(i) sc(sample(cand, n)), numeric(1))
    best <- sc(fit$OPTtrain)
    if (spec$maximize) expect_gte(best, max(rnd)) else expect_lte(best, min(rnd))
  }
})

test_that("a set-only criterion warns when a test set is given", {
  expect_warning(optTrain(pcs, 1:40, 10, test = 41:60, method = "AOpt", min.iter = 30, console = FALSE),
                 "test set is ignored")
})

test_that("GVoverall asks for enough columns", {
  expect_error(optTrain(pcs, 1:40, 10, method = "GVoverall", min.iter = 30, console = FALSE), "at least n.train columns")
})
