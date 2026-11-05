# Criterion values for a given training set and test set. The computation is in C++
# (src/scores.cpp); these wrappers check and convert the input.

.as_score_matrix = function(x, name){
  x = as.matrix(x)
  if(!is.numeric(x)) stop(name, " must be a numeric matrix.")
  if(any(!is.finite(x))) stop(name, " must not contain missing or non-finite values.")
  storage.mode(x) = "double"
  x
}

.lambda_value = function(lambda){
  if(is.null(lambda)) return(0)  # 0 selects the default penalty of the criterion
  if(!is.numeric(lambda) || length(lambda) != 1 || !is.finite(lambda) || lambda <= 0){
    stop("lambda must be a single positive number or NULL.")
  }
  as.numeric(lambda)
}

.score_pair = function(id, X, X0, lambda = NULL){
  score_matrices(id, .as_score_matrix(X, "X"), .as_score_matrix(X0, "X0"), .lambda_value(lambda), FALSE)
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
#' @param lambda A positive number, the ridge penalty, or NULL (default) for the default of the criterion (r-score and PEV: 1 / number of columns; CD: 1).
#' @return A number, the r-score.
#'
#' @references Ou, J.-H., & Liao, C.-T. (2019). Training set determination for genomic selection. Theoretical and Applied Genetics, 132, 2781-2792. <doi:10.1007/s00122-019-03387-0>
#'
#' @export
#' @examples
#' data(geno)
#' \dontrun{r_score(geno[1:50, ], geno[51:100, ])}
r_score = function(X, X0, lambda = NULL) .score_pair(0L, X, X0, lambda)

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
pev_score = function(X, X0, lambda = NULL) .score_pair(1L, X, X0, lambda)

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
cd_score = function(X, X0, lambda = NULL) .score_pair(2L, X, X0, lambda)

.score_set = function(id, X, lambda = NULL){
  X = .as_score_matrix(X, "X")
  score_matrices(id, X, X[1L, , drop = FALSE], .lambda_value(lambda), FALSE)  # X0 is not used
}

#' A-optimality score
#'
#' \eqn{tr[(X'X + \lambda I)^{-1}]}: the trace of the inverse ridge information matrix of the training set. Smaller is better. It uses the training set only.
#'
#' @details See \code{\link{criteria}} for the list of criteria and the package website for one page per criterion.
#'
#' @param X A numeric matrix of the training set (row: individual; column: marker or PC).
#' @param lambda A positive number, the ridge penalty, or NULL (default) for the default of the criterion.
#' @return A number.
#' @references Akdemir, D., & Isidro-Sanchez, J. (2019). Design of training populations for selective phenotyping in genomic prediction. Scientific Reports, 9, 1446. <doi:10.1038/s41598-018-38081-6>
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' data(geno)
#' aopt_score(geno[1:50, ])
aopt_score = function(X, lambda = NULL) .score_set(3L, X, lambda)

#' D-optimality score
#'
#' \eqn{\log\det(X'X + \lambda I)} of the training set. Larger is better (a larger determinant of the information matrix). It uses the training set only.
#'
#' @details See \code{\link{criteria}} for the list of criteria and the package website for one page per criterion.
#'
#' @param X A numeric matrix of the training set (row: individual; column: marker or PC).
#' @param lambda A positive number, the ridge penalty, or NULL (default) for the default of the criterion.
#' @return A number.
#' @references Akdemir, D., & Isidro-Sanchez, J. (2019). Design of training populations for selective phenotyping in genomic prediction. Scientific Reports, 9, 1446. <doi:10.1038/s41598-018-38081-6>
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' data(geno)
#' dopt_score(geno[1:50, ])
dopt_score = function(X, lambda = NULL) .score_set(4L, X, lambda)

#' E-optimality score
#'
#' Largest eigenvalue of \eqn{(X'X + \lambda I)^{-1}}, that is \eqn{1 / (\lambda_{min}(X'X) + \lambda)}. Smaller is better. It uses the training set only. When X has more columns than rows, \eqn{X'X} is singular and the score is \eqn{1/\lambda} for every training set.
#'
#' @details See \code{\link{criteria}} for the list of criteria and the package website for one page per criterion.
#'
#' @param X A numeric matrix of the training set (row: individual; column: marker or PC).
#' @param lambda A positive number, the ridge penalty, or NULL (default) for the default of the criterion.
#' @return A number.
#' @references Akdemir, D., & Isidro-Sanchez, J. (2019). Design of training populations for selective phenotyping in genomic prediction. Scientific Reports, 9, 1446. <doi:10.1038/s41598-018-38081-6>
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' data(geno)
#' eopt_score(geno[1:50, ])
eopt_score = function(X, lambda = NULL) .score_set(5L, X, lambda)

#' PEVmax score
#'
#' Largest prediction error variance among the target rows, with an intercept column added: \eqn{\max_i x_{0i}' (X'X + \lambda I)^{-1} x_{0i}}. Smaller is better.
#'
#' @details See \code{\link{criteria}} for the list of criteria and the package website for one page per criterion.
#'
#' @param X A numeric matrix of the training set (row: individual; column: marker or PC).
#' @param X0 A numeric matrix of the target rows (the test set), with the same columns as X.
#' @param lambda A positive number, the ridge penalty, or NULL (default) for the default of the criterion.
#' @return A number.
#' @references Akdemir, D., & Isidro-Sanchez, J. (2019). Design of training populations for selective phenotyping in genomic prediction. Scientific Reports, 9, 1446. <doi:10.1038/s41598-018-38081-6>
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' data(geno)
#' pevmax_score(geno[1:50, ], geno[51:100, ])
pevmax_score = function(X, X0, lambda = NULL) .score_pair(6L, X, X0, lambda)

#' G-optimality score (GOPTPEV)
#'
#' Largest eigenvalue of the prediction error variance matrix of the target rows, \eqn{X_0 (X'X + \lambda I)^{-1} X_0'}, with an intercept column added. Smaller is better.
#'
#' @details See \code{\link{criteria}} for the list of criteria and the package website for one page per criterion.
#'
#' @param X A numeric matrix of the training set (row: individual; column: marker or PC).
#' @param X0 A numeric matrix of the target rows (the test set), with the same columns as X.
#' @param lambda A positive number, the ridge penalty, or NULL (default) for the default of the criterion.
#' @return A number.
#' @references Akdemir, D., & Isidro-Sanchez, J. (2019). Design of training populations for selective phenotyping in genomic prediction. Scientific Reports, 9, 1446. <doi:10.1038/s41598-018-38081-6>
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' data(geno)
#' gopt_score(geno[1:50, ], geno[51:100, ])
gopt_score = function(X, X0, lambda = NULL) .score_pair(7L, X, X0, lambda)

#' MSPE score of ridge regression
#'
#' The approximation of the expected mean squared prediction error of ridge regression, \eqn{1 + tr[X_0 A A' X_0']/n_0 + tr[(X_0 - X_0 A X)'(X_0 - X_0 A X)]} with \eqn{A = X'(XX' + \lambda I)^{-1}}. Smaller is better. The default penalty is \eqn{\lambda = 1}.
#'
#' @details See \code{\link{criteria}} for the list of criteria and the package website for one page per criterion.
#'
#' @param X A numeric matrix of the training set (row: individual; column: marker or PC).
#' @param X0 A numeric matrix of the target rows (the test set), with the same columns as X.
#' @param lambda A positive number, the ridge penalty, or NULL (default) for the default of the criterion.
#' @return A number.
#' @references Chen, S.-P., Sung, W.-H., & Liao, C.-T. (2024). Constructing training sets for genomic selection to identify superior genotypes in candidate populations. Theoretical and Applied Genetics, 137, 270. <doi:10.1007/s00122-024-04766-y>
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' data(geno)
#' mspe_score(geno[1:50, ], geno[51:100, ])
mspe_score = function(X, X0, lambda = NULL) .score_pair(8L, X, X0, lambda)

#' GVoverall score
#'
#' \eqn{\log\det K_T} with \eqn{K_T = XX'/p}, the relationship matrix of the training set (overall genomic variation). Larger is better. It needs at least as many columns as rows, otherwise \eqn{K_T} is singular and the score is minus infinity.
#'
#' @details See \code{\link{criteria}} for the list of criteria and the package website for one page per criterion.
#'
#' @param X A numeric matrix of the training set (row: individual; column: marker or PC).
#' @return A number.
#' @references Chen, S.-P., Sung, W.-H., & Liao, C.-T. (2024). Constructing training sets for genomic selection to identify superior genotypes in candidate populations. Theoretical and Applied Genetics, 137, 270. <doi:10.1007/s00122-024-04766-y>
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' data(geno)
#' gv_overall_score(geno[1:50, ])
gv_overall_score = function(X) .score_set(9L, X)

#' GVaverage score
#'
#' \eqn{tr(K_T)} with \eqn{K_T = XX'/p} (average genomic variation). Larger is better. It uses the training set only.
#'
#' @details See \code{\link{criteria}} for the list of criteria and the package website for one page per criterion.
#'
#' @param X A numeric matrix of the training set (row: individual; column: marker or PC).
#' @return A number.
#' @references Chen, S.-P., Sung, W.-H., & Liao, C.-T. (2024). Constructing training sets for genomic selection to identify superior genotypes in candidate populations. Theoretical and Applied Genetics, 137, 270. <doi:10.1007/s00122-024-04766-y>
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' data(geno)
#' gv_average_score(geno[1:50, ])
gv_average_score = function(X) .score_set(10L, X)

#' Avg_GRM score
#'
#' Mean of the relationship matrix \eqn{G = XX'/p} between the training set and the target rows (a = 1, b = 0 of the generalized Avg_GRM). Larger is better.
#'
#' @details See \code{\link{criteria}} for the list of criteria and the package website for one page per criterion.
#'
#' @param X A numeric matrix of the training set (row: individual; column: marker or PC).
#' @param X0 A numeric matrix of the target rows (the test set), with the same columns as X.
#' @return A number.
#' @references Fernandez-Gonzalez, J., Akdemir, D., & Isidro y Sanchez, J. (2023). A comparison of methods for training population optimization in genomic selection. Theoretical and Applied Genetics, 136, 30. <doi:10.1007/s00122-023-04265-6>
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' data(geno)
#' avg_grm_score(geno[1:50, ], geno[51:100, ])
avg_grm_score = function(X, X0) .score_pair(11L, X, X0)

#' Avg_GRM_self score
#'
#' Mean of \eqn{G = XX'/p} within the training set. Smaller is better. It uses the training set only.
#'
#' @details See \code{\link{criteria}} for the list of criteria and the package website for one page per criterion.
#'
#' @param X A numeric matrix of the training set (row: individual; column: marker or PC).
#' @return A number.
#' @references Fernandez-Gonzalez, J., Akdemir, D., & Isidro y Sanchez, J. (2023). A comparison of methods for training population optimization in genomic selection. Theoretical and Applied Genetics, 136, 30. <doi:10.1007/s00122-023-04265-6>
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' data(geno)
#' avg_grm_self_score(geno[1:50, ])
avg_grm_self_score = function(X) .score_set(12L, X)

#' Avg_GRM_MinMax score
#'
#' Mean relationship between the training set and the target rows minus the mean relationship within the training set (a = 1, b = 1). Larger is better.
#'
#' @details See \code{\link{criteria}} for the list of criteria and the package website for one page per criterion.
#'
#' @param X A numeric matrix of the training set (row: individual; column: marker or PC).
#' @param X0 A numeric matrix of the target rows (the test set), with the same columns as X.
#' @return A number.
#' @references Fernandez-Gonzalez, J., Akdemir, D., & Isidro y Sanchez, J. (2023). A comparison of methods for training population optimization in genomic selection. Theoretical and Applied Genetics, 136, 30. <doi:10.1007/s00122-023-04265-6>
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' data(geno)
#' avg_grm_minmax_score(geno[1:50, ], geno[51:100, ])
avg_grm_minmax_score = function(X, X0) .score_pair(13L, X, X0)

#' CDmean score (GBLUP form)
#'
#' Mean coefficient of determination of the contrasts between each target row and the target mean, for a GBLUP model with \eqn{K = XX'/p}, a general mean as the only fixed effect, and penalty \eqn{\lambda = \sigma_e^2/\sigma_A^2} (default 1). Larger is better.
#'
#' @details See \code{\link{criteria}} for the list of criteria and the package website for one page per criterion.
#'
#' @param X A numeric matrix of the training set (row: individual; column: marker or PC).
#' @param X0 A numeric matrix of the target rows (the test set), with the same columns as X.
#' @param lambda A positive number, the ridge penalty, or NULL (default) for the default of the criterion.
#' @return A number.
#' @references Rincent, R., et al. (2012). Maximizing the reliability of genomic selection by optimizing the calibration set of reference individuals. Genetics, 192, 715-728. <doi:10.1534/genetics.112.141473>
#'   Chen, S.-P., Sung, W.-H., & Liao, C.-T. (2024). Constructing training sets for genomic selection to identify superior genotypes in candidate populations. Theoretical and Applied Genetics, 137, 270. <doi:10.1007/s00122-024-04766-y>
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' data(geno)
#' cdmean_score(geno[1:50, ], geno[51:100, ])
cdmean_score = function(X, X0, lambda = NULL) .score_pair(14L, X, X0, lambda)

#' CDmean(v2) score
#'
#' Mean of \eqn{diag(K_{PT} W K_{TP}) / diag(K_{PP})}, the heuristic form of CDmean (GBLUP, \eqn{K = XX'/p}, default \eqn{\lambda = 1}). Larger is better.
#'
#' @details See \code{\link{criteria}} for the list of criteria and the package website for one page per criterion.
#'
#' @param X A numeric matrix of the training set (row: individual; column: marker or PC).
#' @param X0 A numeric matrix of the target rows (the test set), with the same columns as X.
#' @param lambda A positive number, the ridge penalty, or NULL (default) for the default of the criterion.
#' @return A number.
#' @references Chen, S.-P., Sung, W.-H., & Liao, C.-T. (2024). Constructing training sets for genomic selection to identify superior genotypes in candidate populations. Theoretical and Applied Genetics, 137, 270. <doi:10.1007/s00122-024-04766-y>
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' data(geno)
#' cdmean2_score(geno[1:50, ], geno[51:100, ])
cdmean2_score = function(X, X0, lambda = NULL) .score_pair(15L, X, X0, lambda)

#' PAM (k-medoids) cost
#'
#' Sum, over the target rows, of the Euclidean distance to the nearest training row: the cost that partitioning around medoids minimizes when the training rows are the medoids. Smaller is better.
#'
#' @details See \code{\link{criteria}} for the list of criteria and the package website for one page per criterion.
#'
#' @param X A numeric matrix of the training set (row: individual; column: marker or PC).
#' @param X0 A numeric matrix of the target rows (the test set), with the same columns as X.
#' @return A number.
#' @references Kaufman, L., & Rousseeuw, P. J. (1990). Finding groups in data: an introduction to cluster analysis. Wiley. <doi:10.1002/9780470316801>
#'   Guo, T., et al. (2019). Optimal designs for genomic selection in hybrid crops. Molecular Plant, 12, 390-401. <doi:10.1016/j.molp.2018.12.022>
#'   Rio, S., et al. (2021). Genomic prediction and training set optimization in a structured Mediterranean oat population. Theoretical and Applied Genetics, 134, 3595-3609. <doi:10.1007/s00122-021-03916-w>
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' data(geno)
#' pam_score(geno[1:50, ], geno[51:100, ])
pam_score = function(X, X0) .score_pair(16L, X, X0)

#' CDmax score
#'
#' Largest, over the target rows, of the coefficient of determination of the ridge model on the columns of \code{X} (the per-individual value averaged by \code{\link{cd_score}}; default \eqn{\lambda = 1}). It protects the best-determined individual, as the CDMAX of STPGA does. Larger is better.
#'
#' @details See \code{\link{criteria}} for the list of criteria and the package website for one page per criterion.
#'
#' @param X A numeric matrix of the training set (row: individual; column: marker or PC).
#' @param X0 A numeric matrix of the target rows (the test set), with the same columns as X.
#' @param lambda A positive number, the ridge penalty, or NULL (default) for the default of the criterion.
#' @return A number.
#' @references Akdemir, D., Sanchez, J. I., & Jannink, J.-L. (2015). Optimization of genomic selection training populations with a genetic algorithm. Genetics Selection Evolution, 47, 38. <doi:10.1186/s12711-015-0116-6>
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' data(geno)
#' cdmax_score(geno[1:50, ], geno[51:100, ])
cdmax_score = function(X, X0, lambda = NULL) .score_pair(17L, X, X0, lambda)
