# Criterion values for a given training set and test set. The computation is in C++
# (src/scores.cpp); these wrappers check and convert the input.

.as_score_matrix = function(x, name){
  x = as.matrix(x)
  if(!is.numeric(x)) stop(name, " must be a numeric matrix.")
  if(any(!is.finite(x))) stop(name, " must not contain missing or non-finite values.")
  storage.mode(x) = "double"
  x
}

.score_pair = function(id, X, X0){
  score_matrices(id, .as_score_matrix(X, "X"), .as_score_matrix(X0, "X0"), FALSE)
}

#' r-score
#'
#' This function calculate r-score <doi:10.1007/s00122-019-03387-0> by given training set and test set.
#'
#' @details The value equals the r-score of v2.4.2 (to its single-precision accuracy of about 1e-7).
#' Since v3 the computation uses p x p matrices when the number of columns p is at most the
#' training set size n, and n x n matrices otherwise, and the result is returned in double
#' precision. The term \code{1 - 1/n0} (n0: test set size) of v2.4.2 is kept as published
#' in the software; see NEWS.md.
#'
#' @author Jen-Hsiang Ou
#'
#' @param X A numeric matrix. The training set genotypic information matrix can be given as genotype matrix (coded as -1, 0, 1) or principal component matrix (row: sample; column: marker or PC).
#' @param X0 A numeric matrix. The test set genotypic information matrix, with the same columns as X.
#' @return A number, the r-score.
#'
#' @references Ou, J.-H., & Liao, C.-T. (2019). Training set determination for genomic selection. Theoretical and Applied Genetics, 132, 2781-2792. <doi:10.1007/s00122-019-03387-0>
#'
#' @export
#' @examples
#' data(geno)
#' \dontrun{r_score(geno[1:50, ], geno[51:100, ])}
r_score = function(X, X0) .score_pair(0L, X, X0)

#' PEV score
#'
#' This function calculate prediction error variance (PEV) score <doi:10.1186/s12711-015-0116-6> by given training set and test set.
#'
#' @author Jen-Hsiang Ou
#'
#' @inheritParams r_score
#' @return A number, the PEV score (smaller is better).
#'
#' @export
#' @examples
#' data(geno)
#' \dontrun{pev_score(geno[1:50, ], geno[51:100, ])}
pev_score = function(X, X0) .score_pair(1L, X, X0)

#' CD-score
#'
#' This function calculate CD-score <doi:10.1186/1297-9686-28-4-359> by given training set and test set.
#'
#' @author Jen-Hsiang Ou
#'
#' @inheritParams r_score
#' @return A number, the CD score (larger is better).
#'
#' @export
#' @examples
#' data(geno)
#' \dontrun{cd_score(geno[1:50, ], geno[51:100, ])}
cd_score = function(X, X0) .score_pair(2L, X, X0)
