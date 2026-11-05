## Submission

This is an update of TSDFGS (CRAN version 2.0) to version 3.0.0.

Main changes (details in NEWS.md):

* The search and the scores run in C++ (Rcpp, RcppEigen). All exported functions keep their
  names and arguments. With `options(TSDFGS.legacy_float = TRUE)` and the same seed, `optTrain()`
  returns the training sets of version 2.4.2 (tested).
* 15 new training set criteria (see `criteria()`).
* The `parallel` package is no longer used. `SSDFGS()`, `FGCM()`, and `nt2r()` can run independent
  searches with OpenMP (`multi.threads` or `threads`); the default is one thread. During
  `R CMD check` (`_R_CHECK_LIMIT_CORES_` set) at most two threads are used.
* Imports are now Rcpp and stats (dplyr, latex2exp, and parallel were removed).

The license stays GPL (>= 3). The package contains no code from other packages; the criteria
are implemented from their published definitions.

## Test environments

* Linux (Apptainer image, R 4.4.3, GCC 16.2.0): `R CMD check --as-cran`, 0 errors, 0 warnings.
* Windows and macOS: not tested before submission; the CRAN checks run after submission.
* Without OpenMP (the default on macOS): all tests pass (compiled with the OpenMP flags removed).

## R CMD check results

0 errors | 0 warnings | 3 notes, all from the check container:

* "unable to verify current time" (the container has no clock service).
* "Compilation used the following non-portable flag(s): -march=nocona" (set in the container's R).
* "checking compiled code": "requires 'nm' to be on the PATH" (not in the container). The native
  routines are registered with `useDynLib(TSDFGS, .registration = TRUE)`.

## Maintainer

The maintainer address changes from jen-hsiang.ou@imbim.uu.se (CRAN, version 2.0) to
oumark.me@outlook.com. Please send the confirmation request to both addresses if needed.

## Reverse dependencies

None (checked with `tools::package_dependencies("TSDFGS", reverse = TRUE)` on the CRAN package database, 2026-11-05).
