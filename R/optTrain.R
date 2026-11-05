####################################################
#                             \-- r-score (RPT)    #
#              \-- target ----\-- pev-score (PPT)  #
#              \              \-- CD-score (CPT)   #
#   \-- pop ---\                                   #
#   \          \              \-- r-score(RPN)     #
#   \          \-- untarget --\-- pev-score(PPN)   #
#   \                         \-- CD-score(CPN)    #
#   \                                              #
#   \                         \-- r-score (RNT)    #
#   \          \-- target ----\-- pev-score (PNT)  #
#   \          \              \-- CD-score (CNT)   #
#   \-- npop --\                                   #
#              \              \-- r-score (RNN)    #
#              \-- untarget --\-- pev-score(PNN)   #
#                             \-- CD-score (CNN)   #
####################################################

#' Optimal training set determination
#'
#' This function is designed for determining optimal training set.
#'
#' All twelve search routes of v2.4.2 (with or without sub-populations, targeted or
#' un-targeted, and the three criteria) are kept. The search moves, their random-number
#' order, the default numbers of iterations, and the stopping rule are the same as in
#' v2.4.2, so that a given seed returns the training set of v2.4.2 (apart from the
#' float-to-double change of the scores, see Details of \code{\link{r_score}}). The
#' r-score is the criterion of Ou and Liao (2019) <doi:10.1007/s00122-019-03387-0>.
#'
#' Only the r-score, PEV, and CD are tested carefully: they reproduce the output recorded from
#' v2.4.2 on all twelve routes. The 15 criteria added in v3 follow their papers and run on all
#' routes, but they are not carefully tested; evaluate them on your own data before relying on them.
#'
#' @author Jen-Hsiang Ou
#'
#' @param geno A numeric matrix of principal components (rows: individuals; columns: PCs).
#' @param cand An integer vector of which rows of individuals are candidates of the training set in the geno matrix.
#' @param n.train The size of the target training set. This could be determined with the help of the ssdfgp function provided in this package.
#' @param subpop A character vector of sub-population's group name. The algorithm will ignore the population structure if it remains NULL. Without \code{test}, it has one entry per candidate (same order as \code{cand}); with \code{test}, it has one entry per row of \code{geno}.
#' @param test An integer vector of which rows of individuals are in the test set in the geno matrix. The algorithm will use an un-target method if it remains NULL.
#' @param method A criterion name; see \code{\link{criteria}} for the list. "rScore" (default), "PEV" and "CD" are the criteria of v2.4.2.
#' @param min.iter Minimum iteration of all methods can be appointed. One should always check if the algorithm is converged or not. A minimum iteration will set by considering the candidate and test set size if it remains NULL.
#' @param console Default: TRUE. Set it to FALSE if you don't want the function printing out the number count of each iteration.
#' @param lambda Default: NULL. A positive ridge penalty; NULL uses the default of each criterion (see \code{\link{criteria}}).
#' @param threads Default: NULL (one thread). The number of threads that Eigen may use inside the matrix products of one search (OpenMP). For the exchange search it only helps when geno has many columns (several hundred; about 1.3 to 1.5 times faster at 404 columns on 2 threads) and the training set found is the same apart from rounding in the last digits. For the targeted r-score and PEV search without sub-populations (a genetic algorithm), the 30 candidate sets of each generation are scored in parallel, and the result does not depend on \code{threads}. To run many independent searches in parallel use \code{\link{SSDFGS}}.
#' @param seed Default: NULL. An integer passed to \code{set.seed} before the search, which changes the global random-number state. Leave it NULL to use the current state.
#'
#' @return This function will return a list with OPTtrain (a vector of chosen optimal training set), TOPscore (highest scores of before iteration), ITERscore (criteria scores of each iteration), n.iter (number of iterations run), and converged (TRUE if the search stopped by the convergence rule, FALSE if it stopped at twice \code{min.iter}).
#'
#' @references Ou, J.-H., & Liao, C.-T. (2019). Training set determination for genomic selection. Theoretical and Applied Genetics, 132, 2781-2792. <doi:10.1007/s00122-019-03387-0>
#'
#' @export
#' @examples
#' data(geno)
#' \dontrun{optTrain(geno, cand = 1:404, n.train = 100)}
#'
optTrain = function(geno, cand, n.train, subpop=NULL, test=NULL, method="rScore", min.iter=NULL, console=TRUE, seed=NULL, lambda=NULL, threads=NULL)
{
  if(length(method) != 1 || !method %in% .criteria_registry$method){
    stop("Method not found. Please choose one from (", paste(.criteria_registry$method, collapse = ", "), "); see criteria().")
  }
  if(!is.null(threads)){
    if(length(threads) != 1 || is.na(threads) || threads < 1) stop("threads must be a single positive integer.")
    threads = .limit_check_threads(as.integer(threads))
    set_eigen_threads(threads)
    on.exit(set_eigen_threads(1L), add = TRUE)
  }
  geno = as.matrix(geno)
  if(is.numeric(geno)) storage.mode(geno) = "double"
  .check_optTrain_input(geno, cand, n.train, test)
  if(!is.null(seed)) set.seed(seed)

  spec = .criterion_spec(method, lambda)
  if(method == "GVoverall" && ncol(geno) < n.train){
    stop("GVoverall needs at least n.train columns in geno (for example all PCs, or the markers).")
  }
  pop = !is.null(subpop)
  target = !is.null(test)
  Nc = length(cand)
  if(is.null(min.iter)){
    min.iter = round(sqrt(Nc * n.train) * .default_iter_factor(pop, target, method))
  }
  .check_min_iter(min.iter)
  tol = .default_tol(pop, target, method)

  if(spec$set_only && target){
    warning("The criterion ", method, " uses the training set only; the test set is ignored.")
  }
  if(!pop && target && method %in% c("rScore", "PEV")){
    fit = .search_ga(geno, cand, n.train, test, spec, min.iter, tol, console, if(is.null(threads)) 1L else as.integer(threads))
  }else if(pop && target){
    fit = .search_pop_target(geno, cand, n.train, subpop, test, spec, min.iter, tol, console)
  }else{
    fit = .search_exchange(geno, cand, n.train, subpop, test, spec, min.iter, tol, console)
  }
  list(
    OPTtrain = as.numeric(fit$sol),
    TOPscore = as.numeric(fit$top),
    ITERscore = as.numeric(fit$iter),
    n.iter = fit$n.iter,
    converged = fit$converged
  )
}
