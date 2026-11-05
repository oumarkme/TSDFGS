# TSDFGS 3.0.0

R stays the interface; the computation moved to C++ (Rcpp and RcppEigen). All exported
functions, their names, and their arguments are those of v2.4.2, so v2 scripts run unchanged.

## New criteria and license

* 15 training set criteria were added (A-, D-, E-, G-optimality, PEVmax, MSPE, GVoverall,
  GVaverage, Avg_GRM, Avg_GRM_self, Avg_GRM_MinMax, CDmean, CDmean(v2), PAM, CDmax). Select one
  with `optTrain(method = )`; score a given set with its `*_score()` function; list them
  with `criteria()`. Each has a page on the website with the original reference, the
  authors' repository, and a license review.
* `optTrain()` has a `lambda` argument (ridge penalty).
* The license is GPL (>= 3), as in v2.4.2. STPGA (GPL-3) and TSODBG (GPL >= 3) were used to check results;
  no code was copied.

## Multithreading

* `SSDFGS()`, `FGCM()`, and `nt2r()` run independent searches in parallel with OpenMP
  (`multi.threads = TRUE` or `threads = n`) instead of `parallel::mclapply`. Each search has
  its own random generator, so values do not depend on the number of threads. They differ
  from the serial run with the same seed. Benchmark: `bench/bench_openmp.R`.
* The package no longer imports `parallel`; the default thread count comes from OpenMP
  (processors available to the process, limited by `OMP_NUM_THREADS`).
* `optTrain(threads = )` sets the threads Eigen may use inside one search (useful for several
  hundred columns).

## Published method unchanged

* `r_score()`, `optTrain()`, and `SSDFGS()` keep the method of Ou & Liao (2019,
  <doi:10.1007/s00122-019-03387-0>). The criterion values equal the v2.4.2 values to about
  1e-7 (v2.4.2 returned single-precision numbers).
* `optTrain()` keeps the search moves, the order of the random numbers, the default number
  of iterations, and the stopping rule. With `options(TSDFGS.legacy_float = TRUE)` and the
  same seed it returns the training sets recorded from v2.4.2 for all twelve routes
  (with or without sub-populations, targeted or not, three criteria). This is tested.
* `r_score()` still uses the scalar `1 - 1/n0` where a centering matrix may be intended. It
  is kept on purpose, so that results do not change.

## Faster

* The scores no longer form n0 x n0 matrices. They use matrices of size min(n, p) (primal or
  dual form), a Cholesky solve instead of an explicit inverse, and a precomputed test-set
  term when the test set is fixed. See the article "How the scores are computed".
* With more columns than training rows (for example markers), the scores use the Gram matrix of all rows, built once per search. On the rice PCs this makes 404 PCs 13 to 106 times faster than v2.4.2 on one core.
* The search (hill climbing) runs in C++: no row copies in R, no growing vectors.
* The score vectors of the search are allocated once.
* Benchmark: `bench/bench_speed.R`; results in the article "Moving from v2 to v3".

## Fixes

* `optTrain()` no longer fails with "the condition has length > 1" when `min.iter` is smaller
  than `5 * sqrt(length(cand) * n.train)` (negative look-back index in the stopping rule).
* Drawing from a group or pool of size 1 no longer samples from `1:x`.
* With sub-populations and no test set, `cand` need not be `1:nrow(geno)`.
* Clear errors for missing values, `n.train` larger than `cand`, and row numbers outside `geno`.
* `min.iter` must be a single finite number of at least 1; `NA` or `Inf` used to run forever.
* A score that is not finite (for example a test set of one individual) stops the search with an
  error instead of returning `NaN`. A test set needs at least 2 individuals.
* `n.train` must be a whole number, and `cand` must not contain duplicates.
* With sub-populations and a test set, a new training individual is never one that is already in
  the training set (v2.4.2 could duplicate a member; `TSDFGS.legacy_float` keeps that draw).
  A sub-population with too few candidates gives a clear error.
* Parallel searches share one Gram matrix instead of copying it for each search.
* Scores are returned in double precision.
* Removed compiled files (`*.o`, `*.so`) from the repository.

## New

* `optTrain(seed = )`; the result also has `n.iter` and `converged`.
* Plots of `SSDFGS()` use the Okabe-Ito colors; ggplot2 is optional (without it the figures are
  `NULL` and the numbers are still returned).
* Dependencies: `dplyr`, `ggplot2`, and `latex2exp` are no longer imported (`dplyr` was never
  used). Imports are now `Rcpp` and `stats`; `RcppEigen` is link-time only.
* Tests (`tests/testthat`) with regression tests against v2.4.2 output. Only the r-score, PEV, and CD (the criteria of v2.4.2) are tested carefully; the 15 added criteria are provided but not carefully tested.
* `env/`: Apptainer recipe (`TSDFGS.def`) and pinned conda list for building and testing.
* License: GPL (>= 3), unchanged from v2.4.2.
* `environment.yml` and `env/Dockerfile`: the tested software versions as a conda environment and as a Docker image (GitHub Container Registry, built by `.github/workflows/docker.yaml`).
* New website (pkgdown) with a getting-started guide and articles on criteria, methods, and
  migration.

