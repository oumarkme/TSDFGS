# Fit logistic growth curve model

A function for fitting logisti growth model

## Usage

``` r
FGCM(geno, nt = NULL, n_iter = NULL, multi.threads = FALSE, threads = NULL)
```

## Arguments

- geno:

  Genotype information saved as a dataframe. Columns represent variants
  (SNPs or PCs).

- nt:

  A numerical vector of training set sample size for estimating logistic
  growth curve parameters

- n_iter:

  Number of simulation of each training set size. Automatically gave a
  suitable number by default.

- multi.threads:

  Default: FALSE. If TRUE, the simulations run in parallel (OpenMP) on
  75% of the available processors (limited by `OMP_NUM_THREADS`), or
  half of them when there are 4 or fewer.

- threads:

  Default: NULL. The number of threads; it overrides `multi.threads`.
  See [`nt2r`](https://www.oumark.net/TSDFGS/reference/nt2r.md) for the
  random numbers of the parallel run.

## Value

Estimation of parameters.

## Examples

``` r
data(geno)
if (FALSE) FGCM(geno) # \dontrun{}
```
