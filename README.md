# TSDFGS

<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/TSDFGS)](https://CRAN.R-project.org/package=TSDFGS) [![R-CMD-check](https://github.com/oumarkme/TSDFGS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/oumarkme/TSDFGS/actions/workflows/R-CMD-check.yaml) [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

> **This is the `v3.0-alpha` branch.** It keeps every function of v2 and the published method, and it is much faster. See [What changed from v2 to v3](#what-changed-from-v2-to-v3).

TSDFGS helps you decide which individuals to phenotype, and how many, before building a genomic selection or prediction model. To predict the genomic estimated breeding values (GEBV) well, the training set should be as strongly related to the test set as possible. The package provides these criteria:

-   **r-score**, the Pearson correlation between GEBVs and phenotypic values of a test set ([Ou & Liao, 2019](https://doi.org/10.1007/s00122-019-03387-0))
-   Prediction error variance (PEV; [Akdemir et al., 2015](https://doi.org/10.1186/s12711-015-0116-6))
-   Generalized coefficient of determination (CD; [Laloë et al., 1996](https://doi.org/10.1186/1297-9686-28-4-359))

and two tools built on them: `optTrain()` finds a good training set, and `SSDFGS()` finds a reasonable training set size.

Method paper: Training set determination for genomic selection ([Ou & Liao, 2019](https://doi.org/10.1007/s00122-019-03387-0)).

## Installation

Development version (this branch):

``` r
# install.packages("remotes")
remotes::install_github("oumarkme/TSDFGS@v3.0-alpha")
```

Stable version from CRAN (v2):

``` r
install.packages("TSDFGS")
```

R \>= 4.0 and a C++17 compiler are needed. The C++ code uses [Rcpp](https://CRAN.R-project.org/package=Rcpp) and [RcppEigen](https://CRAN.R-project.org/package=RcppEigen); [ggplot2](https://CRAN.R-project.org/package=ggplot2) is optional and only draws the figures of `SSDFGS()`. Old versions are on the [releases page](https://github.com/oumarkme/TSDFGS/releases).

## Quick start

``` r
library(TSDFGS)
data(geno)                                   # 404 rice accessions, principal components
g <- as.matrix(geno)[, 1:15]

# 30 of the first 100 accessions to phenotype, to predict accessions 101-150
fit <- optTrain(g, cand = 1:100, n.train = 30, test = 101:150, seed = 1)
fit$OPTtrain                                 # rows to phenotype
r_score(g[fit$OPTtrain, ], g[101:150, ])     # criterion of the chosen set

# How many individuals are enough? (slow; use multi.threads = TRUE on Linux or macOS)
res <- SSDFGS(geno)
res$OC.fig
```

More in the [Getting started guide](https://github.markou.xyz/TSDFGS/articles/getting-started.html).

## Main functions

| Function | Purpose |
|----|----|
| `optTrain()` | Find the optimal training set ([more](https://github.markou.xyz/TSDFGS/reference/optTrain.html)) |
| `SSDFGS()` | Find a reasonable training set size ([more](https://github.markou.xyz/TSDFGS/reference/SSDFGS.html)) |
| `r_score()` | r-score of a given training set ([more](https://github.markou.xyz/TSDFGS/reference/r_score.html)) |
| `pev_score()`, `cd_score()` | PEV and CD of a given training set |
| `FGCM()`, `nt2r()` | Building blocks of `SSDFGS()` |

`cd_score()` and `pev_score()` are also available in the [STPGA](https://CRAN.R-project.org/package=STPGA) package by Deniz Akdemir, with the mixed-model forms of CDmean and PEVmean; [TrainSel](https://CRAN.R-project.org/package=TrainSel) offers a general optimizer. See the article [Choosing a criterion](https://github.markou.xyz/TSDFGS/articles/criteria.html) for other criteria.

## What changed from v2 to v3

**Nothing you call has changed.** All exported functions keep their names and arguments, and the method is the one published in 2019.

| | v2.4.2 | v3.0-alpha |
|----|----|----|
| **Structure** | Search loop in R (12 copies of the same code), scores in C++ | R is the interface; scores and the search run in C++ (Rcpp and RcppEigen) |
| **Speed** | Reference | 2 to 138 times faster per iteration in our benchmark (see below) |
| **Scores** | Single precision (`float`) | Double precision; equal to v2 to about 1e-7 |
| **Reproducing v2** | n/a | `options(TSDFGS.legacy_float = TRUE)` with the same seed returns the v2.4.2 training set for all 12 routes (tested) |
| **Random search** | `set.seed()` only | Also `optTrain(seed = )` |
| **Result of `optTrain()`** | `OPTtrain`, `TOPscore`, `ITERscore` | Same, plus `n.iter` and `converged` |
| **Small `min.iter`** | Error "the condition has length \> 1" below `5 * sqrt(length(cand) * n.train)` | Works |
| **Sub-populations without a test set** | `cand` had to be `1:nrow(geno)` | Any `cand` |
| **Input checks** | Cryptic errors | Clear errors for missing values, too large `n.train`, wrong row numbers |
| **Dependencies** | dplyr, ggplot2, latex2exp, parallel, Rcpp, RcppEigen | Rcpp, parallel (ggplot2 optional) |
| **Figures** | ggplot2, red and blue | ggplot2 if installed, colorblind-safe Okabe-Ito colors |
| **Quality** | No tests | Tests, including regression tests against recorded v2.4.2 output |
| **Build and test** | n/a | Apptainer recipe in `env/` |
| **Website** | One page per function | New site with a guide and articles |

Speed on the rice principal components (one core, same seed and iteration count; median of one run per case):

Milliseconds per iteration of `optTrain()` (rice principal components, n.train = 50 to 100):

| PCs (p) | Case | v2.4.2 | v3.0-alpha | Speed-up |
|---:|----|---:|---:|---:|
| 10 | CD, targeted, n=50, test=300 | 0.075 | 0.015 | 5.1x |
| 10 | PEV, untargeted, n=100 | 0.060 | 0.007 | 8.9x |
| 10 | r-score, untargeted, n=100 | 0.995 | 0.007 | 138.1x |
| 50 | CD, targeted, n=50, test=300 | 0.409 | 0.237 | 1.7x |
| 50 | PEV, untargeted, n=100 | 0.359 | 0.085 | 4.2x |
| 50 | r-score, untargeted, n=100 | 2.661 | 0.136 | 19.6x |
| 404 | CD, targeted, n=50, test=300 | 26.667 | 1.487 | 17.9x |
| 404 | PEV, untargeted, n=100 | 27.599 | 3.639 | 7.6x |
| 404 | r-score, untargeted, n=100 | 17.703 | 4.038 | 4.4x |

The full list is in [NEWS.md](NEWS.md), the article [Moving from v2 to v3](https://github.markou.xyz/TSDFGS/articles/migrating-from-v2.html) has the details, and [update_plan.md](update_plan.md) holds the review of newer training set criteria and the plan for v3.0.

**Kept on purpose.** The r-score code multiplies the test set by `1 - 1/n0` where a centering matrix may be intended. It stays as published until it is checked against the paper.

## Example dataset

Zhao et al. (2011) published the rice genome data used as an example. The raw dataset is on the [Rice Diversity website](http://www.ricediversity.org/data/). The pre-arranged data are in this repository and load with `data(geno)` and `data(subpop)`, or:

``` r
download.file("https://github.com/oumarkme/TSDFGS/raw/main/data/geno.rda", "geno.rda")
load("geno.rda")
```

## Reproducible environment

We tested v3 with R 4.4.3, Rcpp 1.1.2, RcppEigen 0.3.4.0.2, and GCC 16.2.0 (all versions are pinned in `environment.yml`). You can use the same software in either of two ways.

``` bash
# Docker image with R and TSDFGS installed (built by GitHub Actions)
docker run --rm -it -v "$PWD":/work ghcr.io/oumarkme/tsdfgs:v3.0-alpha R

# or a conda / mamba environment
mamba env create -f environment.yml
mamba activate TSDFGS
```

Build the image yourself with `docker build -f env/Dockerfile -t tsdfgs .` from the project root. The image tag is the branch name or release tag. The same packages can also be built as an Apptainer image (see below).

## Development

Build and test with the same software versions in an Apptainer image:

``` bash
env/build_sif.sh <scratch_dir>      # on gs4; writes env/TSDFGS.sif (git-ignored)
apptainer exec --bind $PWD:/work --pwd /work env/TSDFGS.sif Rscript -e 'devtools::test()'
```

## Authors

-   Jen-Hsiang Ou, Ph.D.
    -   Author, maintainer
    -   E-mail: [oumark.me\@outlook.com](mailto:oumark.me@outlook.com)
    -   Department of Medical Genetics, National Taiwan University Hospital, Taipei, Taiwan
-   Po-Ya Wu, Ph.D.
    -   Author
    -   E-mail: [Po-Ya.Wu\@hhu.de](mailto:Po-Ya.Wu@hhu.de)
    -   Institute for Quantitative Genetics and Genomics of Plants, Heinrich Heine University, Düsseldorf, Germany
-   Chen-Tuo Liao, Ph.D.
    -   Author, thesis advisor
    -   E-mail: [ctliao\@ntu.edu.tw](mailto:ctliao@ntu.edu.tw)
    -   Department of Agronomy, National Taiwan University, Taipei, Taiwan

## Citation

Ou, J.-H., & Liao, C.-T. (2019). Training set determination for genomic selection. *Theoretical and Applied Genetics, 132*, 2781-2792. <https://doi.org/10.1007/s00122-019-03387-0>
