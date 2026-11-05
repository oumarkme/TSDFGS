# Simulate r-scores of each training set size

Calculate r-scores (un-target) by in parallel.

## Usage

``` r
nt2r(geno, nt, n_iter = 30, multi.threads = FALSE, threads = NULL)
```

## Arguments

- geno:

  A numeric dataframe of genotype, column represent sites (genotype
  coding as 1, 0, -1)

- nt:

  Numeric. Number of training set size

- n_iter:

  Times of iteration. (default = 30)

- multi.threads:

  Default: FALSE. If TRUE, the iterations run in parallel (OpenMP, any
  system with OpenMP support) on 75% of the available processors
  (limited by `OMP_NUM_THREADS`), or half of them when there are 4 or
  fewer.

- threads:

  Default: NULL. The number of threads; it overrides `multi.threads`.
  The parallel iterations use their own random generators, seeded from
  R's generator. The values then depend on
  [`set.seed()`](https://rdrr.io/r/base/Random.html) but not on
  `threads`, and they are not those of the serial run with the same
  seed.

## Value

A vector of r-scores of each iteration

## Examples

``` r
data(geno)
if (FALSE) nt2r(geno, 50) # \dontrun{}
```
