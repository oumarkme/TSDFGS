# roxygen2::roxygenise()

##### nt2r #####
#' @name nt2r
#' @title Simulate r-scores of each training set size
#' @description Calculate r-scores (un-target) by in parallel.
#' @import parallel
#' @param geno A numeric dataframe of genotype, column represent sites (genotype coding as 1, 0, -1)
#' @param nt Numeric. Number of training set size
#' @param n_iter Times of iteration. (default = 30)
#' @param multi.threads Default: TRUE
#' @return A vector of r-scores of each iteration
#' @examples
#' data(geno)
#' \dontrun{nt2r(geno, 50)}
nt2r = function(geno, nt, n_iter = 30, multi.threads = TRUE){
  
  # Calculate r-scores
  if(multi.threads){
    cores = parallel::detectCores()
    if(cores <= 4){
      n_core = max(round(cores/2), 1)
    }else{
      n_core = round(cores * 3 / 4)
      if(n_iter < cores) n_core = n_iter
    }
    r = unlist(parallel::mclapply(1:n_iter, function(i) return(max(TSDFGS::optTrain(geno, seq(nrow(geno)), nt)$TOPscore)), mc.cores = n_core))
  }else{
    r = sapply(1:n_iter, function(i) return(max(TSDFGS::optTrain(geno, seq(nrow(geno)), nt)$TOPscore)))
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
#' @import dplyr
#' @param geno Genotype information saved as a dataframe. Columns represent variants (SNPs or PCs).
#' @param nt A numerical vector of training set sample size for estimating logistic growth curve parameters
#' @param n_iter Number of simulation of each training set size. Automatically gave a suitable number by default.
#' @param multi.threads Default: TRUE. Set as FALSE if you just want to run it by single thread.
#' @return Estimation of parameters.
#' @examples
#' data(geno)
#' \dontrun{FGCM(geno)}
FGCM = function(geno, nt = NULL, n_iter = NULL, multi.threads = TRUE){
  
  # Basic setting
  if(is.null(nt)){
    nt = round(seq(nrow(geno)/16, nrow(geno)/2, length.out = 10))
  }
  
  if(is.null(n_iter)) n_iter = length(nt)
  
  n_core = detectCores()
  if(multi.threads | n_core > 4){
    n_core = round(n_core * 4 / 5)
  }else{
    n_core = 1
  }
  
  # Start simulation
  result = data.frame(n = rep(nt, each = n_iter))
  
  result$r = unlist(parallel::mclapply(1:nrow(result), function(i){
    r.score = max(TSDFGS::optTrain(geno, seq(nrow(geno)), result$n[i])$TOPscore)
    return(r.score)
  }, mc.cores = n_core))
  
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
#' @import dplyr
#' @import ggplot2
#' @import latex2exp
#' @param geno A numeric data frame carried genotype information (column: PCs, row: sample)
#' @param nt A numeric vector carried training set sizes for r-score simulation.
#' @param n_iter Number of iterations for estimating parameters.
#' @param multi.threads Default (multi.threads = TRUE) use 75% of threads if the computer has more than 4 threads.
#' @return An operating curve and its information.
#' @export
#' @examples
#' data(geno)
#' \dontrun{SSDFGS(geno)}
#' 
SSDFGS = function(geno, nt = NULL, n_iter = NULL, multi.threads = TRUE){
  
  par = FGCM(geno, nt, n_iter, multi.threads)
  r.score = r = NULL
  
  # Operating curve functions
  GC = function(n){
    par$alpha / (1+exp(par$beta - par$gamma * n))
  }
  GC.fit = data.frame(nt = seq(min(par$sim$n), max(par$sim$n)), r.score = sapply(seq(min(par$sim$n), max(par$sim$n)), GC))
  
  RErs = function(n){
    (1 + exp(par$beta - par$gamma * nrow(geno)))/(1 + exp(par$beta - par$gamma * n))
  }
  OC.fit = data.frame(nt = seq(min(par$sim$n), max(par$sim$n)), r.score = sapply(seq(min(par$sim$n), max(par$sim$n)), RErs))
  
  # Plot (Operating curve)
  anno = data.frame(nt = c(), r.score = c())
  if(min(OC.fit$r.score) <= 0.95 & max(OC.fit$r.score) >= 0.95){
    anno = rbind(anno, OC.fit[which.min(abs(OC.fit$r.score - 0.95)),])
  }
  if(min(OC.fit$r.score) <= 0.99 & max(OC.fit$r.score) >= 0.99){
    anno = rbind(anno, OC.fit[which.min(abs(OC.fit$r.score - 0.99)),])
  }
  
  plotOC = ggplot() +
    geom_line(aes(x = nt, y = r.score), OC.fit) +
    geom_point(aes(x = nt, y = r.score), anno, color = "red") +
    annotate("segment", x = -Inf, xend = anno$nt, y = anno$r.score, yend = anno$r.score, color = "red", linetype = "dashed") +
    annotate("segment", x = anno$nt, xend = anno$nt, y = anno$r.score, yend = -Inf, color = "red", linetype = "dashed") +
    geom_text(aes(x = nt, y = r.score, label = paste0(" (", nt, ", ", round(r.score, digits = 2), ")")), anno, color = "red", hjust = "left", vjust = "top") +
    labs(x = expression('Training set size (n'[t]*')'), y = expression('RErs(n'[t]*')'), title = "Operating curve")
  
  # Plot(Growth curve)
  plotGC = ggplot() +
    geom_point(aes(x = n, y = r), data = par$sim, size = 0.8, alpha = 0.5) +
    geom_line(aes(x = nt, y = r.score), data = GC.fit, color = "blue") +
    labs(title = TeX(paste0("Logistic growth curve: $y=\\frac{", round(par$alpha, 3), "}{1+\\exp{(", round(par$beta, 3), "-", round(par$gamma, 3), "x)}}")),
         x = expression('Training set size (n'[t]*')'), y = "r-score")
  
  return(list(OC.fig = plotOC, GC.fig = plotGC, 
              parameter = list(alpha = par$alpha, beta = par$beta, gamma = par$gamma), 
              OC.fit = OC.fit))
}
