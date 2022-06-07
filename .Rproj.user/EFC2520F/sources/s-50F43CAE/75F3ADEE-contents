
# TSDFGS

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/TSDFGS)](https://CRAN.R-project.org/package=TSDFGS)
[![R-CMD-check](https://github.com/oumarkme/TSDFGS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/oumarkme/TSDFGS/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
<!-- badges: end -->

In this package, we provide two useful tools to help determine a reasonable training set for genomic selection or prediction. In order to have a better prediction to the genomic estimated breeding values (GEBV), the training set should be optimized as highly genomic correlated with the test set as possible. Several criteria published previously, including:

- Prediction error variance (PEV; [Akdemir et al., 2015](https://doi.org/10.1186/s12711-015-0116-6))
- Generalized coefficient of determination (CD; [Laloë et al., 1993](https://doi.org/10.1186/1297-9686-28-4-359))

In our research, we provide an alternative criteria, **r-score**, which is derived from Pearson's correlation between GEBVs and phenotypic values of a test set. With the criteria, we could determine both reasonable training set size and an optimal training set for building prediction model. Both functions are provided in our package.

For more information of the method, please check our published article:

- Training set determination for genomic selection ([Ou et al., 2019](https://doi.org/10.1007/s00122-019-03387-0))


## Installation

The development version of TSDFGS can be installed from github (recommend):

``` r
# library(devtools)
install_github("oumark/TSDFGS")
```

You may also install the stable version from CRAN, which the most resent function may not included.

``` r
install.packages("TSDFGS")
```

## Main functions
