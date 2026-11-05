## usethis namespace: start
#' @description We propose an optimality criterion to determine the required training set, r-score, which is derived directly from Pearson's correlation between the genomic estimated breeding values and phenotypic values of the test set <doi:10.1007/s00122-019-03387-0>. This package provides two main functions to determine a good training set and its size.
#' @importFrom Rcpp evalCpp
#' @importFrom stats nls coef SSlogis
#' @useDynLib TSDFGS, .registration = TRUE
"_PACKAGE"

# Keep Eigen single-threaded inside one score; parallel work runs independent searches (OpenMP).
.onLoad = function(libname, pkgname){
  set_eigen_threads(1L)
}
## usethis namespace: end
NULL
