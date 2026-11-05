# Optimal training set determination

This function is designed for determining optimal training set.

## Usage

``` r
optTrain(
  geno,
  cand,
  n.train,
  subpop = NULL,
  test = NULL,
  method = "rScore",
  min.iter = NULL,
  console = TRUE,
  seed = NULL,
  lambda = NULL,
  threads = NULL
)
```

## Arguments

- geno:

  A numeric matrix of principal components (rows: individuals; columns:
  PCs).

- cand:

  An integer vector of which rows of individuals are candidates of the
  training set in the geno matrix.

- n.train:

  The size of the target training set. This could be determined with the
  help of the ssdfgp function provided in this package.

- subpop:

  A character vector of sub-population's group name. The algorithm will
  ignore the population structure if it remains NULL. Without `test`, it
  has one entry per candidate (same order as `cand`); with `test`, it
  has one entry per row of `geno`.

- test:

  An integer vector of which rows of individuals are in the test set in
  the geno matrix. The algorithm will use an un-target method if it
  remains NULL.

- method:

  A criterion name; see
  [`criteria`](https://www.oumark.net/TSDFGS/reference/criteria.md) for
  the list. "rScore" (default), "PEV" and "CD" are the criteria of
  v2.4.2.

- min.iter:

  Minimum iteration of all methods can be appointed. One should always
  check if the algorithm is converged or not. A minimum iteration will
  set by considering the candidate and test set size if it remains NULL.

- console:

  Default: TRUE. Set it to FALSE if you don't want the function printing
  out the number count of each iteration.

- seed:

  Default: NULL. An integer passed to `set.seed` before the search,
  which changes the global random-number state. Leave it NULL to use the
  current state.

- lambda:

  Default: NULL. A positive ridge penalty; NULL uses the default of each
  criterion (see
  [`criteria`](https://www.oumark.net/TSDFGS/reference/criteria.md)).

- threads:

  Default: NULL (one thread). The number of threads that Eigen may use
  inside the matrix products of one search (OpenMP). For the exchange
  search it only helps when geno has many columns (several hundred;
  about 1.3 to 1.5 times faster at 404 columns on 2 threads) and the
  training set found is the same apart from rounding in the last digits.
  For the targeted r-score and PEV search without sub-populations (a
  genetic algorithm), the 30 candidate sets of each generation are
  scored in parallel, and the result does not depend on `threads`. To
  run many independent searches in parallel use
  [`SSDFGS`](https://www.oumark.net/TSDFGS/reference/SSDFGS.md).

## Value

This function will return a list with OPTtrain (a vector of chosen
optimal training set), TOPscore (highest scores of before iteration),
ITERscore (criteria scores of each iteration), n.iter (number of
iterations run), and converged (TRUE if the search stopped by the
convergence rule, FALSE if it stopped at twice `min.iter`).

## Details

All twelve search routes of v2.4.2 (with or without sub-populations,
targeted or un-targeted, and the three criteria) are kept. The search
moves, their random-number order, the default numbers of iterations, and
the stopping rule are the same as in v2.4.2, so that a given seed
returns the training set of v2.4.2 (apart from the float-to-double
change of the scores, see Details of
[`r_score`](https://www.oumark.net/TSDFGS/reference/r_score.md)). The
r-score is the criterion of Ou and Liao (2019)
\<doi:10.1007/s00122-019-03387-0\>.

Only the r-score, PEV, and CD are tested carefully: they reproduce the
output recorded from v2.4.2 on all twelve routes. The 15 criteria added
in v3 follow their papers and run on all routes, but they are not
carefully tested; evaluate them on your own data before relying on them.

## References

Ou, J.-H., & Liao, C.-T. (2019). Training set determination for genomic
selection. Theoretical and Applied Genetics, 132, 2781-2792.
\<doi:10.1007/s00122-019-03387-0\>

## Author

Jen-Hsiang Ou

## Examples

``` r
data(geno)
if (FALSE) optTrain(geno, cand = 1:404, n.train = 100) # \dontrun{}
```
