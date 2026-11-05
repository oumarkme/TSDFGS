# Sample size determination for genomic selection

This function is designed to generate an operating curve for sample size
determination

## Usage

``` r
SSDFGS(geno, nt = NULL, n_iter = NULL, multi.threads = FALSE, threads = NULL)
```

## Arguments

- geno:

  A numeric data frame carried genotype information (column: PCs, row:
  sample)

- nt:

  A numeric vector carried training set sizes for r-score simulation.

- n_iter:

  Number of iterations for estimating parameters.

- multi.threads:

  Default: FALSE. If TRUE, the simulations run in parallel (OpenMP) on
  75% of the available processors (limited by `OMP_NUM_THREADS`), or
  half of them when there are 4 or fewer.

- threads:

  Default: NULL. The number of threads; it overrides `multi.threads`.
  See [`nt2r`](https://www.oumark.net/TSDFGS/reference/nt2r.md) for the
  random numbers of the parallel run.

## Value

A list with the operating curve and growth curve figures (OC.fig and
GC.fig; ggplot objects, or NULL when the ggplot2 package is not
installed), the logistic growth curve parameters, and the operating
curve data (OC.fit).

## Author

Jen-Hsiang Ou & Po-Ya Wu

## Examples

``` r
data(geno)
if (FALSE) SSDFGS(geno) # \dontrun{}
```
