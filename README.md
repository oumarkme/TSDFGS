# TSDFGS

<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/TSDFGS)](https://CRAN.R-project.org/package=TSDFGS) [![R-CMD-check](https://github.com/oumarkme/TSDFGS/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/oumarkme/TSDFGS/actions/workflows/R-CMD-check.yaml) [![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)

<!-- badges: end -->

This package provides two useful tools for determining a proper train genomic selection or prediction training set. To better predict of the genomic estimated breeding values (GEBV), the training set should be optimized as highly genomic correlated with the test set as possible. Several criteria have been published previously, including:

-   Prediction error variance (PEV; [Akdemir et al., 2015](https://doi.org/10.1186/s12711-015-0116-6))
-   Generalized coefficient of determination (CD; [Laloë et al., 1993](https://doi.org/10.1186/1297-9686-28-4-359))

Our research provides an alternative criterion, **r-score**, which is Pearson's correlation between GEBVs and phenotypic values of a test set. We could determine both a reasonable training set size and an optimal training set for building a prediction model with the criteria. Both functions are provided in our package.

For more information on the method, please check our published article:

-   Training set determination for genomic selection ([Ou et al., 2019](https://doi.org/10.1007/s00122-019-03387-0))

## Installation

The development version of TSDFGS can be installed from GitHub (recommend):

``` r
# library(devtools)
install_github("oumarkme/TSDFGS", dependencies = TRUE, force = TRUE)
```

You may also install the stable version from CRAN, which the most recent function may not include.

``` r
install.packages("TSDFGS")
```

-   All functions were developed under r version 4.3.1 and tested in version 4.1.1. However, we recommend you use this package with R version \> 4.0.
-   Rcpp and RcppEigen were used in the package. In addition, the core C++ scripts were published on GitHub for those who want a better performance.

Old versions are available on [GitHub](https://github.com/oumarkme/TSDFGS/releases). You may download the source file and install it locally.

## Main functions

-   `r_score`: available for calculating r-score ([more](https://www.oumark.me/TSDFGS/reference/r_score.html)).
-   `pev_score`: Function for calculating PEV score ([more](https://www.oumark.me/TSDFGS/reference/pev_score.html)).
-   `cd_score`: Function for Calculating CD score ([more](https://www.oumark.me/TSDFGS/reference/cd_score.html)).
-   `optTrain`: Function for determining optimal training set ([more](https://www.oumark.me/TSDFGS/reference/optTrain.html)).
-   `SSDFGS`: Function for determining reasonable training set size ([more](https://www.oumark.me/TSDFGS/reference/SSDFGS.html)).

Note that `cd_score()` and `pev_score()` functions are also available from the [STPGA](https://CRAN.R-project.org/package=STPGA) package by Deniz Akdemir. Try their package for advanced usage.

If you want to install the recent version of the `r_score()` function independently, you may download the `rscore.cpp` script from my GitHub repo and install it by:

``` r
library(Rcpp, RcppEigen)
download.file("https://raw.githubusercontent.com/oumarkme/TSDFGS/main/src/rscore.cpp", "rscore.cpp")
Rcpp::sourceCpp("rscore.cpp")
```

## Example dataset

An example of data provided for testing this package. Zhao et al. (2011) published the rice genome data in their research. The raw dataset is available on the [Rice Diversity website](http://www.ricediversity.org/data/). The pre-arranged dataset is available in this GitHub repository, and you may be loaded in R by

``` r
download.file("https://github.com/oumarkme/TSDFGS/raw/main/data/geno.rda", "geno.rda")
load("geno.rda")
```

## Authors

-   Jen-Hsiang Ou
    -   Author, maintainer
    -   E-mail: [jen-hsiang.ou\@imbim.uu.se](mailto:jen-hsiang.ou@imbim.uu.se){.email}
    -   Department of Medical Biochemistry and Microbiology, Uppsala University, Uppsala, Sweden
-   Po-Ya Wu
    -   Author
    -   E-mail: [Po-Ya.Wu\@hhu.de](mailto:Po-Ya.Wu@hhu.de){.email}
    -   Institute for Quantitative Genetics and Genomics of Plants, Heinrich Heine University, Düsseldorf, Germany
-   Chen-Tuo Liao
    -   Author, thesis advisor
    -   E-mail: [ctliao\@ntu.edu.tw](mailto:ctliao@ntu.edu.tw){.email}
    -   Department of Agronomy, National Taiwan University, Taipei, Taiwan

## Citing this package

If you make use of the TSDFGS package in your research, we would appreciate a citation of the following papers:

-   Ou, J.-H., and Liao, C.-T. Training set determination for genomic selection. Theoretical and Applied Genetics 132, 2781--2792 (2019). <https://doi.org/10.1007/s00122-019-03387-0>
-   Wu, P.-Y., Ou, J.-H., and Liao, C.-T. Sample size determination for training set optimization in genomic prediction. Theoretical and Applied Genetics 136, 57 (2023). <https://doi.org/10.1007/s00122-023-04254-9>
