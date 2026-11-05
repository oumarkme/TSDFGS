# Internal code of optTrain(): defaults of v2.4.2, the stopping rule, and the search routines.
# The exchange search (hill climbing) runs in C++ (src/engine.cpp); the genetic algorithm
# and the sub-population search with a test set run in R and score in C++.

##### Internal helpers #####

# Input checks that v2.4.2 did not make; they only turn later cryptic errors into clear ones.
.check_optTrain_input = function(geno, cand, n.train, test){
  if(!is.numeric(geno)) stop("geno must be a numeric matrix.")
  if(anyNA(geno) || any(!is.finite(geno))) stop("geno must not contain missing or infinite values.")
  if(length(n.train) != 1 || !is.numeric(n.train) || is.na(n.train) || n.train < 2 || n.train != round(n.train)){
    stop("n.train must be a single whole number of at least 2.")
  }
  if(length(cand) < 1 || !is.numeric(cand) || anyNA(cand) || any(!cand %in% seq_len(nrow(geno)))){
    stop("cand must be row numbers of geno.")
  }
  if(anyDuplicated(cand)) stop("cand must not contain duplicated row numbers.")
  if(n.train > length(cand)) stop("n.train is larger than the number of candidates.")
  if(!is.null(test)){
    if(!is.numeric(test) || anyNA(test) || any(!test %in% seq_len(nrow(geno)))) stop("test must be row numbers of geno.")
    if(length(test) < 2) stop("test must contain at least 2 individuals.")
  }
  invisible(TRUE)
}

# min.iter must be a single finite number of at least 1 (NA or Inf would never stop the search).
.check_min_iter = function(min.iter){
  if(length(min.iter) != 1 || !is.numeric(min.iter) || !is.finite(min.iter) || min.iter < 1){
    stop("min.iter must be a single finite number of at least 1.")
  }
  invisible(TRUE)
}

# Criterion code used by the C++ code (0: r-score, 1: PEV, 2: CD) and optimization direction.
# Option TSDFGS.legacy_float = TRUE gives the single-precision scores of v2.4.2, which
# makes a seed return the v2.4.2 training set (used by the regression tests).
.criterion_spec = function(method, lambda = NULL){
  row = .criteria_registry[.criteria_registry$method == method, ]
  list(id = row$id, maximize = row$maximize, set_only = row$set_only, lambda = .lambda_value(lambda),
       legacy = isTRUE(getOption("TSDFGS.legacy_float", FALSE)))
}

.score_rows = function(spec, geno, train, test){
  value = score_rows(spec$id, geno, as.integer(train), as.integer(test), spec$lambda, spec$legacy)
  if(!is.finite(value)){
    stop("The score is not finite; check that the test set has at least 2 individuals and that geno has no constant rows or columns.")
  }
  value
}

# Default minimum iterations, in units of sqrt(candidates x training size). Same values as v2.4.2.
.default_iter_factor = function(pop, target, method){
  if(!method %in% c("rScore", "PEV", "CD")) return(50)  # criteria added in v3
  tab = list(
    npop_untarget = c(PEV = 33, CD = 65, rScore = 80),
    npop_target   = c(PEV = 65, CD = 40, rScore = 40),
    pop_untarget  = c(PEV = 50, CD = 50, rScore = 70),
    pop_target    = c(PEV = 50, CD = 60, rScore = 60)
  )
  key = paste(if(pop) "pop" else "npop", if(target) "target" else "untarget", sep = "_")
  unname(tab[[key]][method])
}

# Convergence tolerance; v2.4.2 used 1e-10 for two CD routes and 1e-6 elsewhere.
.default_tol = function(pop, target, method){
  if(method == "CD" && ((!pop && target) || (pop && !target))) 1e-10 else 1e-6
}

# sample() on a length-1 numeric vector draws from 1:x; these helpers do not.
# For vectors longer than 1 they use the same random draws as sample().
.sample1 = function(x) x[sample.int(length(x), 1L)]
.sampleK = function(x, k) x[sample.int(length(x), k)]

# Stopping rule of v2.4.2: after min.iter iterations stop when the best score has not
# changed by more than tol over the last 5*sqrt(Nc*n) iterations; always stop at 2*min.iter.
# v2.4.2 failed with an error when the look-back index was below 1; here it is skipped.
# (The exchange search in C++ applies the same rule.)
.is_converged = function(top, iter, min.iter, window, tol){
  if(iter > min.iter){
    back = trunc(iter - window)
    if(back >= 1 && abs(top[iter] - top[back]) < tol) return(1L)
  }
  if(iter > min.iter * 2) return(2L)
  0L
}

# Number of training individuals per sub-population (proportional, rounded up, then trimmed to n).
.pop_ratio = function(subpop, n){
  tab = table(subpop)
  ratio = ceiling(n * (as.numeric(tab) / sum(as.numeric(tab))))
  while(sum(ratio) > n){
    i = which(ratio == max(ratio))[1]
    ratio[i] = ratio[i] - 1
  }
  list(pops = names(tab), ratio = ratio)
}

# Exchange (hill-climbing) search in C++: no sub-populations, or sub-populations without a
# test set, or no sub-populations with CD and a test set. Initial set and group sizes are
# drawn here so that the random numbers are those of v2.4.2.
.search_exchange = function(geno, cand, n, subpop, test, spec, min.iter, tol, console){
  pop = !is.null(subpop)
  target = !is.null(test)
  window = sqrt(length(cand) * n) * 5
  group = NULL
  n.groups = 0L
  if(pop){
    subpop = as.character(subpop)
    if(length(cand) != length(subpop)) stop("Input data is not correct.")
    pr = .pop_ratio(subpop, n)
    sol = c()
    for(i in seq_along(pr$pops)){
      sol = c(sol, .sampleK(cand[subpop == pr$pops[i]], pr$ratio[i]))
    }
    group = match(subpop, pr$pops)
    n.groups = length(pr$pops)
  }else{
    sol = .sampleK(cand, n)
  }
  fit = exchange_search(geno, as.integer(cand), as.integer(sol),
                        if(target) as.integer(test) else NULL,
                        if(pop) as.integer(group) else NULL, as.integer(n.groups),
                        spec$id, spec$maximize, as.numeric(min.iter), tol, window,
                        spec$lambda, spec$legacy, console)
  list(sol = fit$sol, top = fit$top, iter = fit$iter, n.iter = fit$n.iter, converged = fit$converged)
}

# Exchange search with sub-populations and a test set (a new member is drawn with
# probability proportional to the group size, which uses sample() with prob).
.search_pop_target = function(geno, cand, n, subpop, test, spec, min.iter, tol, console){
  better = if(spec$maximize) function(a, b) a > b else function(a, b) a < b
  window = sqrt(length(cand) * n) * 5
  subpop = as.character(subpop)
  if(length(subpop) != nrow(geno)) stop("With test, subpop must have one entry per row of geno.")
  pr = .pop_ratio(subpop, n)
  pops = pr$pops
  n.avail = vapply(pops, function(g) sum(subpop[cand] == g), numeric(1))
  if(any(n.avail < pr$ratio)){
    bad = pops[n.avail < pr$ratio]
    stop("Sub-population(s) ", paste(bad, collapse = ", "), " have fewer candidates than the training individuals they need; ",
         "subpop must give the group of every row of geno, and each group needs enough candidates.")
  }
  sol = c()
  for(i in seq_along(pops)){
    sol = c(sol, .sampleK(cand[subpop[cand] == pops[i]], pr$ratio[i]))
  }
  eval_set = function(s) .score_rows(spec, geno, s, test)

  maxlen = ceiling(min.iter * 2) + 3
  top = numeric(maxlen)
  iter.score = numeric(maxlen)
  score = eval_set(sol)
  top[1] = score
  iter.score[1] = score
  iter = 1
  status = 0L
  while(status == 0L){
    p = sample.int(length(pops), 1L, prob = pr$ratio)
    pool = which(subpop[cand] == pops[p])
    # A new member must not already be in the training set. (v2.4.2 could draw a current member,
    # which duplicated it; the legacy option keeps that draw so that the v2.4.2 results are kept.)
    if(!spec$legacy) pool = pool[!cand[pool] %in% sol]
    slots = which(subpop[sol] == pops[p])
    new.sol = sol
    if(length(pool) > 0 && length(slots) > 0){
      value = cand[.sample1(pool)]
      new.sol[.sample1(slots)] = value
    }
    new.score = eval_set(new.sol)
    if(better(new.score, score)){
      sol = new.sol
      score = new.score
    }
    iter.score[iter + 1] = new.score
    top[iter + 1] = score
    if(console) cat(iter, "..", sep = "")
    iter = iter + 1
    status = .is_converged(top, iter, min.iter, window, tol)
  }
  list(sol = sol, top = top[2:iter], iter = iter.score[2:iter], n.iter = iter - 1L,
       converged = status == 1L)
}

# Genetic algorithm of the un-structured targeted routes for PEV and r-score: 30 solutions,
# 3 replaced per generation, 3% of each solution mutated.
.search_ga = function(geno, cand, n, test, spec, min.iter, tol, console, threads = 1L){
  maximize = spec$maximize
  window = sqrt(length(cand) * n) * 5
  npop = 30L
  k.mut = ceiling(n * .03)
  # More columns than training rows: the dual form needs only the Gram matrix of all rows,
  # built once here.
  gram = if(spec$id <= 2L && ncol(geno) > n && nrow(geno) <= 5000) tcrossprod(geno) else NULL
  sol = matrix(NA_integer_, n, npop)
  for(i in seq_len(npop)) sol[, i] = .sampleK(cand, n)
  storage.mode(sol) = "integer"
  score = score_many(spec$id, geno, sol, as.integer(test), spec$lambda, spec$legacy, threads, gram)

  best = if(maximize) max else min
  maxlen = ceiling(min.iter * 2) + 3
  top = numeric(maxlen)
  iter.score = numeric(npop + maxlen)
  top[1] = best(score)
  iter.score[seq_len(npop)] = score
  iter = 1
  status = 0L
  while(status == 0L){
    if(maximize){
      elite = which(rank(score) > 27)
      del = sample(seq_len(npop)[-elite], 3, prob = ((1 / score[-elite]) / sum(1 / score[-elite])))
    }else{
      elite = which(rank(score) < 4)
      del = sample(seq_len(npop)[-elite], 3, prob = (score[-elite] / sum(score[-elite])))
    }
    for(i in seq_along(del)){
      par = sample(seq_len(npop)[-del], 2)
      sol[, del[i]] = .sampleK(unique(c(sol[, par[1]], sol[, par[2]])), n)
    }
    # The candidate sets are drawn first, in the order of v2.4.2, and scored together.
    trial = sol
    for(i in seq_len(npop)){
      if(i %in% del){
        value = .sampleK(cand[!cand %in% sol[, i]], k.mut)
        sol[sample.int(n, k.mut), i] = value
        trial[, i] = sol[, i]
      }else{
        new.sol = sol[, i]
        value = .sampleK(cand[!cand %in% new.sol], k.mut)
        new.sol[sample.int(n, k.mut)] = value
        trial[, i] = new.sol
      }
    }
    trial.score = score_many(spec$id, geno, trial, as.integer(test), spec$lambda, spec$legacy, threads, gram)
    for(i in seq_len(npop)){
      if(i %in% del){
        score[i] = trial.score[i]
      }else if((maximize && trial.score[i] > score[i]) || (!maximize && trial.score[i] < score[i])){
        sol[, i] = trial[, i]
        score[i] = trial.score[i]
      }
    }
    iter.score[npop + iter] = mean(score)
    top[iter + 1] = best(score)
    if(console) cat(iter, "..", sep = "")
    iter = iter + 1
    status = .is_converged(top, iter, min.iter, window, tol)
  }
  winner = which(score == best(score))[1]
  list(sol = sol[, winner], top = top[2:iter], iter = iter.score[2:(npop + iter - 1)],
       n.iter = iter - 1L, converged = status == 1L)
}
