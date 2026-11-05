# roxygen2::roxygenise()

##### Parallel simulation #####
# Number of OpenMP threads for n_tasks independent searches. `threads` overrides everything;
# otherwise multi.threads = TRUE uses 75% of the available processors (half when there are 4 or fewer);
# n_procs() respects OMP_NUM_THREADS and OMP_THREAD_LIMIT.
.n_threads = function(multi.threads, threads, n_tasks){
  if(!is.null(threads)){
    if(length(threads) != 1 || is.na(threads) || threads < 1) stop("threads must be a single positive integer.")
    return(max(1L, min(.limit_check_threads(as.integer(threads)), as.integer(n_tasks))))
  }
  if(!isTRUE(multi.threads)) return(1L)
  cores = n_procs()
  n = if(cores <= 4) max(round(cores / 2), 1) else round(cores * 3 / 4)
  max(1L, min(.limit_check_threads(as.integer(n)), as.integer(n_tasks)))
}

# During R CMD check (_R_CHECK_LIMIT_CORES_ is set) at most two threads are used, as the CRAN
# policy requires.
.limit_check_threads = function(n){
  lim = Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
  if(nzchar(lim) && !identical(tolower(lim), "false")) min(n, 2L) else n
}

# Best r-score of each of length(nt) untargeted searches, run in parallel by C++ (OpenMP).
# Every search has its own random generator, seeded from R's generator, so the values depend
# on set.seed() but not on the number of threads.
.simulate_rscore = function(geno, nt, threads){
  geno = as.matrix(geno)
  if(is.numeric(geno)) storage.mode(geno) = "double"
  cand = seq_len(nrow(geno))
  .check_optTrain_input(geno, cand, max(nt), NULL)
  spec = .criterion_spec("rScore")
  n_cand = length(cand)
  size = sqrt(n_cand * nt)
  min.iter = round(size * .default_iter_factor(FALSE, FALSE, "rScore"))
  sols = lapply(nt, function(k) as.integer(.sampleK(cand, k)))
  seeds = as.numeric(sample.int(.Machine$integer.max, length(nt)))
  # Longest searches first (larger training sets run longer), so that no thread starts a long
  # search last; the seeds belong to the searches, so the order does not change the values.
  ord = order(nt, decreasing = TRUE)
  fit = exchange_search_multi(geno, as.integer(cand), sols[ord], NULL, spec$id, spec$maximize,
                              as.numeric(min.iter)[ord], .default_tol(FALSE, FALSE, "rScore"),
                              (size * 5)[ord], spec$lambda, spec$legacy, seeds[ord],
                              as.integer(threads))
  score = numeric(length(nt))
  score[ord] = fit$score
  score
}

##### nt2r #####
#' @name nt2r
#' @title Simulate r-scores of each training set size
#' @description Calculate r-scores (un-target) by in parallel.
#' @param geno A numeric dataframe of genotype, column represent sites (genotype coding as 1, 0, -1)
#' @param nt Numeric. Number of training set size
#' @param n_iter Times of iteration. (default = 30)
#' @param multi.threads Default: FALSE. If TRUE, the iterations run in parallel (OpenMP, any system with OpenMP support) on 75\% of the available processors (limited by \code{OMP_NUM_THREADS}), or half of them when there are 4 or fewer.
#' @param threads Default: NULL. The number of threads; it overrides \code{multi.threads}. The parallel iterations use their own random generators, seeded from R's generator. The values then depend on \code{set.seed()} but not on \code{threads}, and they are not those of the serial run with the same seed.
#' @return A vector of r-scores of each iteration
#' @export
#' @examples
#' data(geno)
#' \dontrun{nt2r(geno, 50)}
nt2r = function(geno, nt, n_iter = 30, multi.threads = FALSE, threads = NULL){

  # Calculate r-scores
  if(isTRUE(multi.threads) || !is.null(threads)){
    r = .simulate_rscore(geno, rep(nt, n_iter), .n_threads(multi.threads, threads, n_iter))
  }else{
    r = sapply(1:n_iter, function(i) return(max(TSDFGS::optTrain(geno, seq(nrow(geno)), nt, console=F)$TOPscore)))
  }

  return(r)
}
##### END #####



##### FGCM #####
# Fit growth curve model
#' @name FGCM
#' @title Fit logistic growth curve model
#' @description A function for fitting logisti growth model
#' @param geno Genotype information saved as a dataframe. Columns represent variants (SNPs or PCs).
#' @param nt A numerical vector of training set sample size for estimating logistic growth curve parameters
#' @param n_iter Number of simulation of each training set size. Automatically gave a suitable number by default.
#' @param multi.threads Default: FALSE. If TRUE, the simulations run in parallel (OpenMP) on 75\% of the available processors (limited by \code{OMP_NUM_THREADS}), or half of them when there are 4 or fewer.
#' @param threads Default: NULL. The number of threads; it overrides \code{multi.threads}. See \code{\link{nt2r}} for the random numbers of the parallel run.
#' @return Estimation of parameters.
#' @export
#' @examples
#' data(geno)
#' \dontrun{FGCM(geno)}
FGCM = function(geno, nt = NULL, n_iter = NULL, multi.threads = FALSE, threads = NULL){

  # Basic setting
  if(is.null(nt)){
    nt = round(seq(nrow(geno)/16, nrow(geno)/2, length.out = 10))
  }

  if(is.null(n_iter)) n_iter = length(nt)

  # Start simulation
  result = data.frame(n = rep(nt, each = n_iter))
  use_parallel = isTRUE(multi.threads) || !is.null(threads)

  if(use_parallel){
    result$r = .simulate_rscore(geno, result$n, .n_threads(multi.threads, threads, nrow(result)))
  }else{
    for(i in 1:nrow(result)){
      result$r[i] = max(TSDFGS::optTrain(geno, seq(nrow(geno)), result$n[i], console=F)$TOPscore)
    }
  }


  # Fit growth model
  fit = stats::coef(stats::nls(r ~ SSlogis(n, asym, xmid, scal), data = result))
  alpha = as.numeric(fit[1])
  gamma = as.numeric(1 / fit[3])
  beta = as.numeric(fit[2] * gamma)

  # Return parameters
  return(list(alpha=alpha, beta=beta, gamma=gamma, sim = result))
}



##### SSDFGS #####
#' Sample size determination for genomic selection
#'
#' This function is designed to generate an operating curve for sample size determination
#'
#' @author Jen-Hsiang Ou & Po-Ya Wu
#'
#' @param geno A numeric data frame carried genotype information (column: PCs, row: sample)
#' @param nt A numeric vector carried training set sizes for r-score simulation.
#' @param n_iter Number of iterations for estimating parameters.
#' @param multi.threads Default: FALSE. If TRUE, the simulations run in parallel (OpenMP) on 75\% of the available processors (limited by \code{OMP_NUM_THREADS}), or half of them when there are 4 or fewer.
#' @param threads Default: NULL. The number of threads; it overrides \code{multi.threads}. See \code{\link{nt2r}} for the random numbers of the parallel run.
#' @return A list with the operating curve and growth curve figures (OC.fig and GC.fig; ggplot objects, or NULL when the ggplot2 package is not installed), the logistic growth curve parameters, and the operating curve data (OC.fit).
#' @export
#' @examples
#' data(geno)
#' \dontrun{SSDFGS(geno)}
#'
SSDFGS = function(geno, nt = NULL, n_iter = NULL, multi.threads = FALSE, threads = NULL){

  cat('Fitting logistic growth curve model...')
  par = FGCM(geno, nt, n_iter, multi.threads, threads)

  # Operating curve functions
  GC = function(n){
    par$alpha / (1+exp(par$beta - par$gamma * n))
  }
  GC.fit = data.frame(nt = seq(min(par$sim$n), max(par$sim$n)), r.score = sapply(seq(min(par$sim$n), max(par$sim$n)), GC))

  RErs = function(n){
    (1 + exp(par$beta - par$gamma * nrow(geno)))/(1 + exp(par$beta - par$gamma * n))
  }
  OC.fit = data.frame(nt = seq(min(par$sim$n), max(par$sim$n)), r.score = sapply(seq(min(par$sim$n), max(par$sim$n)), RErs))

  # Points of the operating curve where RErs reaches 0.95 and 0.99
  anno = data.frame(nt = c(), r.score = c())
  if(min(OC.fit$r.score) <= 0.95 & max(OC.fit$r.score) >= 0.95){
    anno = rbind(anno, OC.fit[which.min(abs(OC.fit$r.score - 0.95)),])
  }
  if(min(OC.fit$r.score) <= 0.99 & max(OC.fit$r.score) >= 0.99){
    anno = rbind(anno, OC.fit[which.min(abs(OC.fit$r.score - 0.99)),])
  }

  # Figures need ggplot2, which is an optional package; all numbers are returned without it.
  cat('Plotting...')
  plotOC = plotGC = NULL
  if(requireNamespace("ggplot2", quietly = TRUE)){
    plotOC = .plot_operating_curve(OC.fit, anno)
    plotGC = .plot_growth_curve(par, GC.fit)
  }else{
    message("Package ggplot2 is not installed: the figures are NULL; the operating curve data are in OC.fit.")
  }

  return(list(OC.fig = plotOC, GC.fig = plotGC,
              parameter = list(alpha = par$alpha, beta = par$beta, gamma = par$gamma),
              OC.fit = OC.fit))
}
