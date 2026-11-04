# TSDFGS update plan (v2.4.2 to v3.0)

Date: 2026-11-04. Status: the plan was written first; the `v3.0-alpha` branch now implements the parts marked below.

## 0. Status on the v3.0-alpha branch

| Plan item | State |
|---|---|
| Phase 0: tests, baseline, remove compiled files | Done. `tests/testthat` has regression tests against output recorded from v2.4.2 (all 12 `optTrain` routes reproduce it exactly with `options(TSDFGS.legacy_float = TRUE)`). |
| Phase 1: exact faster scores (4.1) | Done in C++ with RcppEigen: primal or dual form, Cholesky solves, double precision. A1, A3, A4, A5 fixed. |
| Phase 2: C++ search engine (4.3) | Done for the exchange search (five of the six non-genetic routes). The genetic algorithm and the sub-population search with a test set still run in R and score in C++. A6, A7, A8, A9, A10 fixed. Rank-1 updates (4.2) not done; the untargeted r-score and PEV avoid copying the test rows by subtracting from a precomputed total. |
| Phase 3: OpenMP, seeds, faster `SSDFGS` (4.4, 4.5) | `seed` done. OpenMP, warm starts, and a faster `SSDFGS` are not done. |
| Phase 4-5: new criteria | Not done (kept out of this alpha so that the published method stays unchanged). |
| A2: centering matrix in the r-score | Kept as published; open. |
| Phase 6: docs, website, release | README with v2 to v3 changelog, NEWS, pkgdown site with articles: done. |

Speed measured against v2.4.2 is in `bench/bench_speed_results.tsv` and in the article "Moving from v2 to v3".

## 1. Summary

TSDFGS (Training Set Determination For Genomic Selection) ships three criteria (r-score, prediction error variance [PEV], and coefficient of determination [CD]) and one search routine (`optTrain`). Since the package was written, the field has added model-based criteria for multi-environment and top-genotype goals, relationship- and diversity-based criteria, and general-purpose optimizers (STPGA, an R package using a genetic algorithm [GA], and TrainSel). A 2023 benchmark found that the best criteria were CDmean for targeted sets and minimum average relationship for untargeted sets (Fernández-González et al., 2023). The package has neither the second one nor a CDmean that matches the literature definition.

The package is slow for three reasons that can be fixed without changing results:

1. Every score is computed from scratch, with n x n or n0 x n0 matrices, although the r-score and PEV can be written with p x p matrices (p = number of principal components [PCs]).
2. The search loop runs in R and copies sub-matrices on each iteration.
3. `SSDFGS` repeats the full search 100 times (default) to fit one curve.

I checked the first point numerically (section 4.1). For p = 10 PCs the p x p form returned the same r-score as the current C++ code (0.9396054) and ran 140 times faster (0.028 s vs 0.0002 s per call, n = 300, n0 = 400). For p = 300 it was not faster, so rank-1 updates are also needed (section 4.2).

## 2. Audit of the current code

| # | Finding | Where | Effect |
|---|---|---|---|
| A1 | Scores are returned as `float` | `src/rscore.cpp`, `src/pevScore.cpp`, `src/cdScore.cpp` (return type) | About 7 significant digits. The `1e-10` convergence tolerance (`R/optTrain.R:176`, `:288`) cannot be resolved. |
| A2 | `IJ <- I - I/n0` is a scaled identity, not the centering matrix `I - J/n0` | `src/rscore.cpp` (`IJ`) | The r-score is not computed on centered test genotypes. Check against Ou and Liao (2019) before changing; if the paper centers, this is a bug that changes results. |
| A3 | `cd_score` forms an n0 x n0 matrix to read its diagonal | `src/cdScore.cpp` | O(n0^2) memory and time; row sums of squares give the same number. |
| A4 | `pev_score` forms an n0 x n0 matrix to take its trace | `src/pevScore.cpp` | Same; `trace(B x0'x0)` is p x p. |
| A5 | All Eigen inputs are passed by value, and `.inverse()` is used | all three `.cpp` files | A copy per call and a slower, less stable inverse than a Cholesky solve. |
| A6 | Twelve near-identical branches (pop x target x method) | `R/optTrain.R:47-447` | About 400 duplicated lines; bugs must be fixed twelve times. |
| A7 | `geno[sol, ]` sub-setting and `c(top.score, ...)` growth inside the loop | `R/optTrain.R` (every branch) | R-level copies each iteration; vector growth is O(iterations^2). |
| A8 | Convergence index `top.score[iter - sqrt(Nc*n)*5]` is a non-integer | `R/optTrain.R:70` and equivalents | Silently truncated; the stopping rule is hard to reason about. Hard cap at `2*min.iter`. |
| A9 | In the sub-population untargeted branches, `which(subpop == pops[p])` indexes positions in `cand`, but `sol` holds row numbers of `geno` | `R/optTrain.R:238`, `:276`, `:313` | Correct only when `cand` is `1:nrow(geno)`. Needs a test with a shuffled `cand`. |
| A10 | `sample(x, 1)` where `x` has length 1 samples from `1:x` | many lines, for example `:58` | Wrong draw if only one candidate remains in a group. |
| A11 | Stale build products `src/*.o` and `src/TSDFGS.so` are in the tree | `src/` | Remove and add to `.gitignore`. |
| A12 | `parallel::mclapply` only; no Windows path | `R/SSDFGS.R` | Parallel runs are Linux and macOS only. |
| A13 | No tests, no benchmark, no seed argument | package root | Nothing guards a refactor. |
| A14 | `geno` documentation says 404 PCs for 404 samples | `R/data.R` | With p = n, the p x p form has no advantage, so benchmarks on this data hide the gain. |

## 3. Criteria available now

"Targeted" means the test set is known when the training set is chosen. Status refers to TSDFGS v2.4.2. Verification labels follow the `trusted-sources` rules (section 8).

### 3.1 Model-based (expected accuracy from the mixed model or ridge model)

| Criterion | Idea | Status in TSDFGS | Source |
|---|---|---|---|
| r-score | Closed form of Pearson r between GEBV and phenotype in the test set (ridge) | Implemented | Ou & Liao, 2019 |
| PEVmean | Mean prediction error variance of contrasts between candidates and the population mean | Implemented as PCA-ridge variant | Rincent et al., 2012; Akdemir et al., 2015 |
| CDmean | Mean expected reliability (coefficient of determination) of the same contrasts | Implemented as a PCA-ridge variant (`lambda = 1`); not the G-matrix form with `lambda = (1 - h2)/h2` | Rincent et al., 2012 |
| CDmean(v2) | Algebraically equivalent to CDmean, much cheaper to compute | Missing | Chen et al., 2024 |
| MSPE-ridge | Ridge-regression expected mean squared prediction error (MSPE); recommended when the candidate set has no strong structure and the goal is finding top genotypes | Missing | Chen et al., 2024 |
| GV_overall | D-optimality-like criterion that maximizes total genomic variation; preferred with strong sub-population structure | Missing | Chen et al., 2024 |
| GV_average | Ranking rule used to shrink a large candidate pool before a heuristic search | Missing | Chen et al., 2024 |
| CDmean.MET | CDmean for multi-environment trials (MET) with a genotype-by-environment (G x E) term | Missing | Rio et al., 2022 (cited in Liu & Liao, 2026; not read directly) |
| A-, D-, E-optimality and related design criteria | Classical optimal-design functions of the training information matrix | Missing | Akdemir, 2017 (preprint) |
| G x E x trait kernel PEV | PEV on a genomic-by-environment-by-trait (and enviromic) kernel | Missing | Gevartosky et al., 2021 |
| Hybrid design (choose crosses to phenotype) | PEV-type criteria over all possible single crosses | Missing | Kadam et al., 2021; Heslot & Feoktistov, 2017 |

### 3.2 Relationship- and diversity-based (no phenotype model)

| Criterion | Idea | Status | Source |
|---|---|---|---|
| Avg_GRM_self | Minimize the mean genomic relationship within the training set | Missing; best untargeted criterion in a seven-data-set benchmark | Fernández-González et al., 2023 |
| Uniform coverage | Sample the genetic (PC) space uniformly, optionally within strata | Missing | Bustos-Korts et al., 2016 |
| Stratified sampling and stratified CDmean | Allocate training size by sub-population, then optimize | Partly (proportional allocation inside `optTrain`) | Isidro et al., 2015 |
| PAM (partitioning around medoids), FURS, MaxCD | Cluster-, graph-, and connectedness-based selection | Missing | Rio et al., 2021; Tibbs-Cortes et al., 2022 |
| Core Hunter metrics (A-NE, E-NE, allelic richness) | Core-collection diversity metrics with multi-objective weights | Missing | De Beukelaer et al., 2018 |

### 3.3 Evaluation metrics (not criteria, but needed to compare them)

Normalized discounted cumulative gain (NDCG), Spearman correlation, and rank sum ratio measure how well a training set finds the top genotypes, which Pearson r does not (Chen et al., 2024; Liu & Liao, 2026). Add them to a benchmark helper, not to the search.

### 3.4 What the recent evidence says about choosing

* Targeted optimization beat untargeted, with the largest gap at low heritability (no single gap value is reported in the abstracts I read). CDmean was best but slow (Fernández-González et al., 2023; Akdemir et al., 2019).
* A training set of 50-55% of the candidates (targeted) or 65-85% (untargeted) reached 95% of the maximum accuracy (Fernández-González et al., 2023). This bears on `SSDFGS` defaults.
* PEVmean and CDmean beat partitioning around medoids (PAM) under strong structure in oat (Rio et al., 2021), but PEVmean was worst in a multi-parental wheat population (Zhang et al., 2021). No criterion wins everywhere, so the package should expose several and a comparison helper.
* The objective matters: criteria that maximize mean accuracy are not the best for picking top genotypes (Chen et al., 2024).
* Gap: I found no peer-reviewed deep-learning or active-learning method for training set selection in genomic prediction. One review only proposes it as a future direction (Alemu et al., 2024). I did not search livestock-specific criteria in depth because the PubMed server hit its rate limit; run that search before release.

## 4. Speed and efficiency plan

### 4.1 Cheaper scores (exact, same value)

Write each score with p x p matrices, where `B = (X'X + lambda I)^-1`:

* r-score: `X'(XX' + lambda I)^-1 = B X'` (push-through identity). With `S0 = c^2 X0'X0` (`c = 1 - 1/n0` reproduces the current code), `q1 = n0 - 1 + tr(S0)`, `q2 = tr(B X'X B S0) + tr(M' S0 M)`, `q12 = c tr(X0'X0 M)`, `M = B X'X`. Verified numerically against `src/rscore.cpp` logic (difference 0 to 7 digits).
* PEV: `tr(B x0'x0)`; `x0'x0` is computed once for a fixed test set.
* CD: per-row quadratic forms `x_i' Q x_i / x_i'x_i` with `Q = B X'X B`; no n0 x n0 matrix.

Use the dual (n x n) form when n < p, and the primal form when p < n. The bundled `geno` has p = 404 = n, so also add an `n.pc` option (keep the leading PCs) and report that it changes the criterion.

### 4.2 Exchange moves with rank-1 updates

One search step swaps one individual. Then `X'X` changes by two rank-1 terms, and `B` can be updated with Sherman-Morrison-Woodbury in O(p^2) instead of an O(p^3) inverse. For untargeted search the test set is "candidates minus training set", so keep totals of `sum x` and `sum xx'` and subtract the training-set part; no re-subsetting. Expected cost per step: O(p^2) to O(p^3), independent of n and n0.

### 4.3 Move the loop to C++

* One C++ engine: inputs are `geno`, candidate indices, test indices or "complement", strata, criterion id, n.train, seed, and stopping rules. Return best set, best-score trace, and iteration trace.
* Pass `Eigen::Map<const MatrixXd>`; use `LDLT` solves; return `double`.
* Replace the twelve R branches with one dispatcher in R (A6). Keep the current argument names.
* Fix A8 with an explicit rule: stop when the best score has not improved by more than `tol` for `patience` iterations, with `max.iter` as a hard cap.
* Keep `top.score` and `iter.score` in preallocated vectors (A7).

### 4.4 Parallelism and reproducibility

* OpenMP (`SHLIB_OPENMP_CXXFLAGS`) over GA population members and independent restarts; `threads` argument defaults to 1, never to all cores. Works on Windows, replacing `mclapply` (A12).
* Per-thread random streams seeded from one user `seed`, so results are identical across thread counts.
* Add `-O2` and optional `-march=native` through `Makevars` only for local builds (CRAN forbids it by default).

### 4.5 Faster `FGCM` and `SSDFGS`

* Default: 10 sizes x 10 repeats = 100 full searches. Use fewer repeats with early stopping, run them in parallel (4.4), and cap `max.iter` using a pilot run.
* Warm start: seed the search for size n+1 from the best set of size n.
* Use `nls` with a fallback start (`SSlogis` can fail on flat curves) and return the fit status.
* Use the Fernández-González et al. (2023) size fractions as the plotted reference range.

### 4.6 Large candidate pools

* Pre-filter with a cheap ranking (GV_average, Chen et al., 2024) before the heuristic search.
* Optional low-rank `geno` through truncated SVD (`RSpectra`) when n is above a few thousand.
* Optional GA and tabu memory as in STPGA (Akdemir, 2017) as a second optimizer; the exchange algorithm stays the default.

## 5. API design for v3.0

```r
optTrain(geno, cand, n.train, test = NULL, subpop = NULL,
         criterion = c("rScore", "CDmean", "CDmean2", "PEV", "AvgGRM",
                       "MSPE", "GVoverall", "Uniform"),
         optimizer = c("exchange", "ga"), lambda = NULL, n.pc = NULL,
         max.iter = NULL, tol = 1e-8, patience = NULL,
         threads = 1L, seed = NULL, verbose = FALSE)
```

* `method` stays as a deprecated alias of `criterion` for two releases.
* Return value keeps `OPTtrain`, `TOPscore`, `ITERscore`; add `criterion`, `seed`, `converged`, `n.iter`.
* New: `compareCriteria()` (run several criteria on one split and report selection-focused metrics), `scoreSet()` (score a user-supplied set), a user-supplied criterion hook.
* Export `r_score`, `pev_score`, `cd_score` unchanged in signature.
* Move `Imports` that are used only for plotting (`dplyr`, `ggplot2`, `latex2exp`) to `SSDFGS` internals and drop `dplyr` if unused.

## 6. Validation and tests

1. Equality tests: new scorers versus the v2.4.2 scorers on 20 random sets, tolerance 1e-6 (float precision of the old code); then versus an R reference with 1e-10.
2. Edge cases: p > n, p < n, n.train = 1, one individual left in a stratum (A10), shuffled `cand` (A9), `test` overlapping `cand`, missing values in `geno` (error).
3. Update-rule test: after 1000 exchange moves, the rank-1-updated `B` equals a fresh inverse to 1e-8.
4. Reproducibility: same `seed` gives the same set for `threads` = 1, 4.
5. Statistical check: on the bundled rice PCs (and a public wheat set, if available), compare selected sets from each criterion with random sets by cross-validated prediction accuracy and by NDCG; one table. State the data source and version; do not claim a ranking beyond what this run shows.
6. A2 decision: reproduce the paper's r-score on its example before and after any centering change; document the choice.

## 7. Roadmap

| Phase | Work | Output | Risk |
|---|---|---|---|
| 0 | Add tests and baseline benchmark for v2.4.2; remove `.o`/`.so` (A11); decide A2 | `tests/`, `bench/`, baseline timings | A2 may change published numbers |
| 1 | New double-precision scorers (4.1) with equality tests; fix A1, A3, A4, A5 | Same results, faster per call | Low |
| 2 | C++ engine with exchange rank-1 updates (4.2, 4.3); R dispatcher replaces 12 branches; fix A8-A10 | `optTrain` v3 core | Medium; stopping rule changes run length |
| 3 | OpenMP, seeds, Windows support (4.4); faster `FGCM`/`SSDFGS` (4.5) | Threaded runs | Medium; reproducibility across threads |
| 4 | New criteria: CDmean (G-matrix form), CDmean(v2), AvgGRM, MSPE-ridge, GV_overall, uniform coverage | Criteria table in section 3 | Needs the papers read in full to match definitions |
| 5 | Multi-environment and hybrid criteria (CDmean.MET, crosses) | Optional, after phase 4 | Needs Rio et al. (2022) read directly |
| 6 | Pkgdown vignette, README, NEWS, version 3.0.0; submit | Release | Low |

Benchmark protocol for phases 0-3: run on gs4 (not on this WSL machine), at most 50 cores, `threads` set explicitly, scratch in a private `mktemp -d /raid/tmp/TSDFGS.XXXXXX` directory, results copied to `bench/` in the project. Data sets: bundled rice PCs (404 x 404 and truncated to 10, 50 PCs), plus simulated 5,000 x 50 and 20,000 x 100 PCs. Report wall time at equal final criterion value, median of 5 seeds.

## 8. References and verification status

Labels: `fetched` = record read through PubMed or Consensus with title, year, and journal matching; `snippet` = abstract only; `secondary` = found only as a citation in another paper; `preprint` = Tier 2. Queried 2026-11-04 through PubMed and Consensus (both open access to abstracts). DOIs are given where the record returned one.

* Ou, J.-H., & Liao, C.-T. (2019). Training set determination for genomic selection. *Theoretical and Applied Genetics, 132*, 2781-2792. PMID 31267147, doi:10.1007/s00122-019-03387-0. fetched (abstract).
* Rincent, R., et al. (2012). Maximizing the reliability of genomic selection by optimizing the calibration set of reference individuals. *Genetics*. PMID 22865733. snippet.
* Akdemir, D., Sanchez, J. I., & Jannink, J.-L. (2015). Optimization of genomic selection training populations with a genetic algorithm. *Genetics Selection Evolution, 47*, 38. PMID 25943105, doi:10.1186/s12711-015-0116-6. fetched (abstract).
* Isidro, J., et al. (2015). Training set optimization under population structure in genomic selection. *Theoretical and Applied Genetics, 128*, 145-158. PMID 25367380, doi:10.1007/s00122-014-2418-4. fetched (abstract).
* Bustos-Korts, D., et al. (2016). Improvement of predictive ability by uniform coverage of the target genetic space. *G3, 6*, 3733-3747. PMID 27672112, doi:10.1534/g3.116.035410. fetched (abstract).
* Akdemir, D. (2017). STPGA: Selection of training populations with a genetic algorithm. bioRxiv. preprint; used only for the optimizer idea.
* Heslot, N., & Feoktistov, V. (2017). Optimization of selective phenotyping and population design for genomic prediction. *Journal of Agricultural, Biological and Environmental Statistics*. snippet (Consensus); PMID not found by citation lookup.
* De Beukelaer, H., et al. (2018). Core Hunter 3: flexible core subset selection. *BMC Bioinformatics*. PMID 29855322. snippet.
* Akdemir, D., et al. (2019). Design of training populations for selective phenotyping in genomic prediction. *Scientific Reports*. snippet (Consensus); citation lookup was ambiguous between two PMIDs, so the PMID is unconfirmed.
* Akdemir, D., Rio, S., & Isidro y Sánchez, J. (2021). TrainSel: an R package for selection of training populations. *Frontiers in Genetics, 12*, 655287. PMID 34025720, doi:10.3389/fgene.2021.655287. fetched (abstract).
* Gevartosky, R., et al. (2021). Enviromic-based kernels may optimize resource allocation with multi-trait multi-environment genomic prediction for tropical maize. *BMC Plant Biology*. snippet.
* Kadam, D. C., Rodriguez, O. R., & Lorenz, A. J. (2021). Optimization of training sets for genomic prediction of early-stage single crosses in maize. *Theoretical and Applied Genetics, 134*, 687-699. PMID 33398385, doi:10.1007/s00122-020-03722-w. fetched (abstract).
* Rio, S., et al. (2021). Genomic prediction and training set optimization in a structured Mediterranean oat population. *Theoretical and Applied Genetics, 134*, 3595-3609. PMID 34341832, doi:10.1007/s00122-021-03916-w. fetched (abstract).
* Zhang, W., et al. (2021). Evaluation of genomic prediction for Fusarium head blight resistance with a multi-parental population. *Biology, 10*, 756. PMID 34439988, doi:10.3390/biology10080756. fetched (abstract).
* Rio, S., et al. (2022). CDmean.MET (multi-environment CDmean). secondary: known only from Liu & Liao (2026); title, journal, and DOI not verified. Read directly before phase 5.
* Tibbs-Cortes, L. E., et al. (2022). Genomic prediction of tocochromanols in exotic-derived maize. *The Plant Genome, 16*, e20286. PMID 36575809, doi:10.1002/tpg2.20286. fetched (abstract).
* Fernández-González, J., et al. (2023). A comparison of methods for training population optimization in genomic selection. *Theoretical and Applied Genetics*. PMID 36892603. snippet.
* Chen, S.-P., et al. (2024). Constructing training sets for genomic selection to identify superior genotypes in candidate populations. *Theoretical and Applied Genetics*. snippet (Consensus); citation lookup returned many matches, so the PMID and DOI are unconfirmed.
* Alemu, A., et al. (2024). Genomic selection in plant breeding: key factors shaping two decades of progress. *Molecular Plant*. snippet; review.
* Liu, Z.-J., & Liao, C.-T. (2026). Optimizing training sets for genomic selection to identify superior genotypes across multiple environments. *G3, 16*(4). PMID 41668481, doi:10.1093/g3journal/jkag031. fetched (abstract).

Not verified or not found: Rio et al. (2022) bibliographic record; PMIDs for Heslot and Feoktistov (2017), Akdemir et al. (2019), Chen et al. (2024); livestock-specific criteria; any deep-learning training set selection method. Abstracts support only what they state. Criterion definitions in section 3 that come from abstracts (CDmean(v2), MSPE-ridge, GV_overall, GV_average, Avg_GRM_self) must be checked against the full text before coding.

## 9. Software facts checked locally (2026-11-04)

* Files read: `R/optTrain.R`, `R/SSDFGS.R`, `src/rscore.cpp`, `src/pevScore.cpp`, `src/cdScore.cpp`, `DESCRIPTION`, `NAMESPACE`, `R/data.R`.
* Numerical check, R on this machine, seed 1, n = 30, n0 = 50, p = 8: current-code r-score 0.9396054; p x p form with `c = 1 - 1/n0` 0.9396054; p x p form with true centering 0.9224974 (this difference is finding A2).
* Timing, n = 300, n0 = 400, mean of 5 calls: p = 10, 0.028 s vs 0.0002 s; p = 300, 0.172 s vs 0.212 s (no gain).
* Check script: `/tmp/claude-1000/-lab-nas-home-user-Mark-TSDFGS/43f525c3-6331-4088-98e4-19690bb96ade/scratchpad/chk.R`.
* PEV and CD forms in section 4.1 are derived but not yet tested numerically; phase 1 tests cover them.
