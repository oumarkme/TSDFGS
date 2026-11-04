# TSDFGS 3.0.0.9000 (development version, branch v3.0-alpha)

R stays the interface; the computation moved to C++ (Rcpp and RcppEigen). All exported
functions, their names, and their arguments are those of v2.4.2, so v2 scripts run unchanged.

## Published method unchanged

* `r_score()`, `optTrain()`, and `SSDFGS()` keep the method of Ou & Liao (2019,
  <doi:10.1007/s00122-019-03387-0>). The criterion values equal the v2.4.2 values to about
  1e-7 (v2.4.2 returned single-precision numbers).
* `optTrain()` keeps the search moves, the order of the random numbers, the default number
  of iterations, and the stopping rule. With `options(TSDFGS.legacy_float = TRUE)` and the
  same seed it returns the training sets recorded from v2.4.2 for all twelve routes
  (with or without sub-populations, targeted or not, three criteria). This is tested.
* `r_score()` still uses the scalar `1 - 1/n0` where a centering matrix may be intended. It
  is kept on purpose until checked against the paper (`update_plan.md`, finding A2).

## Faster

* The scores no longer form n0 x n0 matrices. They use matrices of size min(n, p) (primal or
  dual form), a Cholesky solve instead of an explicit inverse, and a precomputed test-set
  term when the test set is fixed. See the article "How the scores are computed".
* The search (hill climbing) runs in C++: no row copies in R, no growing vectors.
* The score vectors of the search are allocated once.
* Benchmark: `bench/bench_speed.R`; results in the article "Moving from v2 to v3".

## Fixes

* `optTrain()` no longer fails with "the condition has length > 1" when `min.iter` is smaller
  than `5 * sqrt(length(cand) * n.train)` (negative look-back index in the stopping rule).
* Drawing from a group or pool of size 1 no longer samples from `1:x`.
* With sub-populations and no test set, `cand` need not be `1:nrow(geno)`.
* Clear errors for missing values, `n.train` larger than `cand`, and row numbers outside `geno`.
* Scores are returned in double precision.
* Removed compiled files (`*.o`, `*.so`) from the repository.

## New

* `optTrain(seed = )`; the result also has `n.iter` and `converged`.
* Plots of `SSDFGS()` use the Okabe-Ito colors; ggplot2 is optional (without it the figures are
  `NULL` and the numbers are still returned).
* Dependencies: `dplyr`, `ggplot2`, and `latex2exp` are no longer imported (`dplyr` was never
  used). Imports are now `Rcpp` and `parallel`; `RcppEigen` is link-time only.
* Tests (`tests/testthat`) with regression tests against v2.4.2 output.
* `env/`: Apptainer recipe (`TSDFGS.def`) and pinned conda list for building and testing.
* `environment.yml` and `env/Dockerfile`: the tested software versions as a conda environment and as a Docker image (GitHub Container Registry, built by `.github/workflows/docker.yaml`).
* New website (pkgdown) with a getting-started guide and articles on criteria, methods, and
  migration.
* `update_plan.md`: literature review of training set criteria and the plan for v3.0.

## Known issues kept from v2.4.2

* In the targeted sub-population route, a new individual may duplicate one already in the
  training set. Planned for a later alpha.
