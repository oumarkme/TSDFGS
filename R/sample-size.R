# roxygen2::roxygenise()

##### nt2r #####
#' @name nt2r
#' @title Simulate r-scores of each training set size
#' @description Calculate r-scores (un-target) by in parallel.
#' @import parallel
#' @param geno A numeric dataframe of genotype, column represent sites (genotype coding as 1, 0, -1)
#' @param nt Numeric. Number of training set size
#' @param n_iter Times of iteration. (default = 30)
#' @param multi.threads Default: FALSE. Multi-thread function is only avalyble for mac or linux systems.
#' @return A vector of r-scores of each iteration
#' @export
#' @examples
#' data(geno)
#' \dontrun{nt2r(geno, 50)}
nt2r = function(geno, nt, n_iter = 30, multi.threads = FALSE){

  # Calculate r-scores
  if(multi.threads & .Platform$OS.type == 'unix'){
    cores = parallel::detectCores()
    if(cores <= 4){
      n_core = max(round(cores/2), 1)
    }else{
      n_core = round(cores * 3 / 4)
      if(n_iter < cores) n_core = n_iter
    }
    r = unlist(parallel::mclapply(1:n_iter, function(i) return(max(TSDFGS::optTrain(geno, seq(nrow(geno)), nt, console=F)$TOPscore)), mc.cores = n_core))
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
#' @import parallel
#' @param geno Genotype information saved as a dataframe. Columns represent variants (SNPs or PCs).
#' @param nt A numerical vector of training set sample size for estimating logistic growth curve parameters
#' @param n_iter Number of simulation of each training set size. Automatically gave a suitable number by default.
#' @param multi.threads Default: FALSE. Multi-thread function is only avalyble for mac or linux systems.
#' @return Estimation of parameters.
#' @export
#' @examples
#' data(geno)
#' \dontrun{FGCM(geno)}
FGCM = function(geno, nt = NULL, n_iter = NULL, multi.threads = FALSE){

  # Basic setting
  if(is.null(nt)){
    nt = round(seq(nrow(geno)/16, nrow(geno)/2, length.out = 10))
  }

  if(is.null(n_iter)) n_iter = length(nt)

  n_core = detectCores()
  if(multi.threads & n_core > 4 & .Platform$OS.type == 'unix'){
    n_core = round(n_core * 4 / 5)
  }else{
    n_core = 1
  }

  # Start simulation
  result = data.frame(n = rep(nt, each = n_iter))

  if(n_core > 1){
      result$r = unlist(parallel::mclapply(1:nrow(result), function(i){
      r.score = max(TSDFGS::optTrain(geno, seq(nrow(geno)), result$n[i], console=F)$TOPscore)
      return(r.score)
    }, mc.cores = n_core))
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
#' @import parallel
#' @param geno A numeric data frame carried genotype information (column: PCs, row: sample)
#' @param nt A numeric vector carried training set sizes for r-score simulation.
#' @param n_iter Number of iterations for estimating parameters.
#' @param multi.threads Default: False. If TRUE, this function will use 75\% of threads if the computer has more than 4 threads. Multi-thread computing is only avalible in mac and linux environments.
#' @return A list with the operating curve and growth curve figures (OC.fig and GC.fig; ggplot objects, or NULL when the ggplot2 package is not installed), the logistic growth curve parameters, and the operating curve data (OC.fit).
#' @export
#' @examples
#' data(geno)
#' \dontrun{SSDFGS(geno)}
#'
SSDFGS = function(geno, nt = NULL, n_iter = NULL, multi.threads = FALSE){

  if(.Platform$OS.type != 'unix') multi.threads = FALSE

  cat('Fitting logistic growth curve model...')
  par = FGCM(geno, nt, n_iter, multi.threads)

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
