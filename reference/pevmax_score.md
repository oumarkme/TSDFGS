# PEVmax score

Largest prediction error variance among the target rows, with an
intercept column added: \\\max_i x\_{0i}' (X'X + \lambda I)^{-1}
x\_{0i}\\. Smaller is better.

## Usage

``` r
pevmax_score(X, X0, lambda = NULL)
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

Akdemir, D., & Isidro-Sanchez, J. (2019). Design of training populations
for selective phenotyping in genomic prediction. Scientific Reports, 9,
1446. \<doi:10.1038/s41598-018-38081-6\>

## Author

Jen-Hsiang Ou

## Examples

``` r
data(geno)
pevmax_score(geno[1:50, ], geno[51:100, ])
#> [1] 12487101
```
