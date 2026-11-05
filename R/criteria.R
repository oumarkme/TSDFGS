# Registry of the training set criteria. The C++ code uses the integer codes (src/scores.h).
# One row per criterion; the website has one page per row (vignettes/articles/criterion-*.Rmd).

.criteria_registry = data.frame(
  method   = c("rScore", "PEV", "CD", "AOpt", "DOpt", "EOpt", "PEVmax", "GOpt", "MSPE", "GVoverall",
               "GVaverage", "AvgGRM", "AvgGRMself", "AvgGRMminmax", "CDmean", "CDmeanV2", "PAM", "CDmax"),
  id       = 0:17,
  maximize = c(TRUE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, TRUE,
               TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE),
  set_only = c(FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE,
               TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE),
  score_fun = c("r_score", "pev_score", "cd_score", "aopt_score", "dopt_score", "eopt_score",
                "pevmax_score", "gopt_score", "mspe_score", "gv_overall_score", "gv_average_score",
                "avg_grm_score", "avg_grm_self_score", "avg_grm_minmax_score", "cdmean_score",
                "cdmean2_score", "pam_score", "cdmax_score"),
  type     = c("model-based", "model-based", "model-based", "optimal design", "optimal design",
               "optimal design", "model-based", "optimal design", "model-based", "genomic variation",
               "genomic variation", "relationship", "relationship", "relationship", "model-based",
               "model-based", "clustering", "model-based"),
  name     = c("r-score", "PEVmean (ridge, PCs)", "CD (ridge, PCs)", "A-optimality", "D-optimality",
               "E-optimality", "PEVmax", "G-optimality (GOPTPEV)", "MSPE of ridge regression",
               "GVoverall", "GVaverage", "Avg_GRM", "Avg_GRM_self", "Avg_GRM_MinMax", "CDmean",
               "CDmean(v2)", "Partitioning around medoids (PAM)", "CDmax (ridge, PCs)"),
  stringsAsFactors = FALSE
)

#' List of training set criteria
#'
#' Returns the criteria that \code{\link{optTrain}} can optimize, with the name to give to its
#' \code{method} argument, the function that scores a given training set, the direction of
#' optimization, and whether the criterion uses a test set.
#'
#' @details Criteria marked \code{uses_test = FALSE} depend on the training set only. When a
#' test set is given to \code{optTrain}, such a criterion ignores it (with a warning).
#' One page per criterion on the package website gives the definition, the original
#' reference, and the license review.
#'
#' @return A data frame with the columns \code{method}, \code{name}, \code{type},
#' \code{score_fun}, \code{maximize}, and \code{uses_test}.
#' @author Jen-Hsiang Ou
#' @export
#' @examples
#' criteria()
criteria = function(){
  r = .criteria_registry
  data.frame(method = r$method, name = r$name, type = r$type, score_fun = r$score_fun,
             maximize = r$maximize, uses_test = !r$set_only, stringsAsFactors = FALSE)
}
