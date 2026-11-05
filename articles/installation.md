# Installation and environments

This page lists what TSDFGS needs, the versions we tested, and four ways
to get a working setup. Pick the one that fits your computer; the
results are the same in all of them.

## Requirements

| Item | Needed | Notes |
|----|----|----|
| R | 4.0 or later | Tested with R-devel, release, and oldrel on every push. |
| C++ compiler | C++17 | GCC or Clang on Linux and macOS, Rtools on Windows. The package is compiled when you install it. |
| Rcpp, RcppEigen | Installed before TSDFGS | R installs them for you from CRAN. They are used at build time. |
| stats | Part of R |  |
| ggplot2 | Optional | Only for the figures of [`SSDFGS()`](https://www.oumark.net/TSDFGS/reference/SSDFGS.md). |
| OpenMP | Optional | Used by the parallel options `multi.threads` and `threads`. Linux and Windows have it. The default macOS compiler does not, so there the package runs on one thread and returns the same values. |

TSDFGS uses one thread unless you ask for more. The environment variable
`OMP_NUM_THREADS` limits the number of threads that
`multi.threads = TRUE` may use.

## Versions we tested

| Software  | Version   |
|-----------|-----------|
| R         | 4.4.3     |
| Rcpp      | 1.1.2     |
| RcppEigen | 0.3.4.0.2 |
| GCC       | 16.2.0    |

These versions are pinned in
[`environment.yml`](https://github.com/oumarkme/TSDFGS/blob/main/environment.yml)
and are the ones inside the Docker image. GitHub Actions also checks
every push on Linux (R-devel, release, oldrel), macOS, and Windows.

Check what you have installed:

``` r

packageVersion("TSDFGS")
sessionInfo()   # R, platform, and loaded package versions; add to bug reports
```

## 1. Install in R

``` r

# install.packages("remotes")
remotes::install_github("oumarkme/TSDFGS")          # current main branch
remotes::install_github("oumarkme/TSDFGS@v3.0.0")   # a given release

# Source package attached to each GitHub release
install.packages("https://github.com/oumarkme/TSDFGS/releases/download/v3.0.0/TSDFGS_3.0.0.tar.gz",
                 repos = NULL)

# CRAN, once version 3 is accepted there
install.packages("TSDFGS")
```

On Windows install
[Rtools](https://cran.r-project.org/bin/windows/Rtools/) first. On macOS
install the command line tools with `xcode-select --install`.

## 2. Docker

The image has R, the compilers, and TSDFGS installed. It is built from
this repository by GitHub Actions and tested before it is published.

``` bash
docker pull ghcr.io/oumarkme/tsdfgs:latest      # current main branch
docker pull ghcr.io/oumarkme/tsdfgs:3.0.0       # a release

# Interactive R in the current directory
docker run --rm -it -v "$PWD":/work ghcr.io/oumarkme/tsdfgs:latest R

# Run a script
docker run --rm -v "$PWD":/work ghcr.io/oumarkme/tsdfgs:latest Rscript my_analysis.R
```

The files of the current directory appear in `/work` inside the
container. Use a release tag in your own analyses, because `latest`
moves. To build the image yourself, run
`docker build -f env/Dockerfile -t tsdfgs .` in a clone of the
repository.

## 3. conda or mamba

A small environment with R, the compiler, and the packages TSDFGS needs,
then TSDFGS from GitHub:

``` bash
mamba create -n tsdfgs -c conda-forge r-base=4.4 r-rcpp r-rcppeigen r-ggplot2 r-remotes gxx_linux-64 make
mamba activate tsdfgs
Rscript -e 'remotes::install_github("oumarkme/TSDFGS@v3.0.0")'
```

The `gxx_linux-64` compiler is for Linux. On macOS use `clangxx_osx-64`
(Intel) or `clangxx_osx-arm64` (Apple silicon) instead.

To get exactly the versions in the table above, clone the repository and
use the pinned file. It also contains the tools to develop and test the
package (devtools, roxygen2, pkgdown, testthat):

``` bash
git clone https://github.com/oumarkme/TSDFGS.git
cd TSDFGS
mamba env create -f environment.yml
mamba activate TSDFGS
R CMD INSTALL .
```

## 4. Apptainer or Singularity

On a shared server without Docker, use the same image:

``` bash
apptainer exec --bind "$PWD":/work --pwd /work docker://ghcr.io/oumarkme/tsdfgs:3.0.0 Rscript my_analysis.R
```

The repository also has a recipe, `env/TSDFGS.def`, and a script,
`env/build_sif.sh`, that build the image from the pinned conda file.

## Older versions

| Version | Where | Use it for |
|----|----|----|
| 1.0 | [CRAN archive](https://cran.r-project.org/src/contrib/Archive/TSDFGS/TSDFGS_1.0.tar.gz) | The package used in the manuscript (Ou & Liao, 2019). |
| 2.4.2 | `remotes::install_github("oumarkme/TSDFGS@v2.4.2")` | The last version before the C++ rewrite. |
| 3.0.0 | This page | Current version. |

To get the training sets of version 2.4.2 with version 3, set
`options(TSDFGS.legacy_float = TRUE)` and use the same seed. [Moving
from v2 to
v3](https://www.oumark.net/TSDFGS/articles/migrating-from-v2.md)
explains the differences.
