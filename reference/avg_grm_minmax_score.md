# Avg_GRM_MinMax score

Mean relationship between the training set and the target rows minus the
mean relationship within the training set (a = 1, b = 1). Larger is
better.

## Usage

``` r
avg_grm_minmax_score(X, X0)
```

## Arguments

- X:

  A numeric matrix of the training set (row: individual; column: marker
  or PC).

- X0:

  A numeric matrix of the target rows (the test set), with the same
  columns as X.

## Value

A number.

## Details

See [`criteria`](https://www.oumark.net/TSDFGS/reference/criteria.md)
for the list of criteria and the package website for one page per
criterion.

## References

Fernandez-Gonzalez, J., Akdemir, D., & Isidro y Sanchez, J. (2023). A
comparison of methods for training population optimization in genomic
selection. Theoretical and Applied Genetics, 136, 30.
\<doi:10.1007/s00122-023-04265-6\>

## Author

Jen-Hsiang Ou

## Examples

``` r
data(geno)
avg_grm_minmax_score(geno[1:50, ], geno[51:100, ])
#> [1] -1.570958
```
