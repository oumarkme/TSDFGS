# Avg_GRM score

Mean of the relationship matrix \\G = XX'/p\\ between the training set
and the target rows (a = 1, b = 0 of the generalized Avg_GRM). Larger is
better.

## Usage

``` r
avg_grm_score(X, X0)
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
avg_grm_score(geno[1:50, ], geno[51:100, ])
#> [1] -0.158039
```
