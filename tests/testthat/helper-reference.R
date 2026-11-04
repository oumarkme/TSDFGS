# Direct R transcription of the v2.4.2 C++ scorers (src/rscore.cpp, pevScore.cpp, cdScore.cpp),
# in double precision. They form the n0 x n0 matrices on purpose: they are the reference
# that the faster v3 scorers must match.
ref_r_score <- function(X, X0) {
  nr <- nrow(X); nc <- ncol(X); n0r <- nrow(X0)
  A <- t(X) %*% solve(X %*% t(X) + diag(nr) / nc)
  IJ <- diag(n0r) - diag(n0r) / n0r  # v2.4.2: scaled identity, not the centering matrix
  q1 <- n0r - 1 + sum((IJ %*% X0)^2)
  IJX0A <- IJ %*% X0 %*% A
  IJX0AX <- IJX0A %*% X
  q2 <- sum(IJX0A^2) + sum(IJX0AX^2)
  q12 <- sum(diag(t(X0) %*% IJX0AX))
  q12 / sqrt(q1 * q2)
}

ref_pev_score <- function(X, X0) {
  p <- ncol(X)
  x <- cbind(1, X); x0 <- cbind(1, X0)
  sum(diag(x0 %*% solve(t(x) %*% x + diag(p + 1) / p) %*% t(x0)))
}

ref_cd_score <- function(X, X0) {
  p <- ncol(X)
  cd <- X0 %*% solve(t(X) %*% X + diag(p)) %*% t(X)
  mean(diag(cd %*% t(cd)) / diag(X0 %*% t(X0)))
}

load_fixture <- function() {
  readRDS(testthat::test_path("fixtures", "baseline_v242.rds"))
}

# Rounds to single precision (R has no float type).
as_float <- function(x) {
  readBin(writeBin(as.double(x), raw(), size = 4L), "numeric", size = 4L, n = length(x))
}
