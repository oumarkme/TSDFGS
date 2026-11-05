# CDmax score

Largest, over the target rows, of the coefficient of determination of
the ridge model on the columns of `X` (the per-individual value averaged
by [`cd_score`](https://www.oumark.net/TSDFGS/reference/cd_score.md);
default \\\lambda = 1\\). It protects the best-determined individual, as
the CDMAX of STPGA does. Larger is better.

## Usage

``` r
cdmax_score(X, X0, lambda = NULL)
```

## Arguments

- X:

  A numeric matrix of the training set (row: individual; column: marker
  or PC).

- X0:

  A numeric matrix of the target rows (the test set), with the same
  columns as X.

- lambda:

  A positive number, the ridge penalty, or NULL (default) for the
  default of the criterion.

## Value

A number.

## Details

See [`criteria`](https://www.oumark.net/TSDFGS/reference/criteria.md)
for the list of criteria and the package website for one page per
criterion.

## References

Akdemir, D., Sanchez, J. I., & Jannink, J.-L. (2015). Optimization of
genomic selection training populations with a genetic algorithm.
Genetics Selection Evolution, 47, 38. \<doi:10.1186/s12711-015-0116-6\>

## Author

Jen-Hsiang Ou

## Examples

``` r
data(geno)
cdmax_score(geno[1:50, ], geno[51:100, ])
#> [1] 4.920842e-05
```
