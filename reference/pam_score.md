# PAM (k-medoids) cost

Sum, over the target rows, of the Euclidean distance to the nearest
training row: the cost that partitioning around medoids minimizes when
the training rows are the medoids. Smaller is better.

## Usage

``` r
pam_score(X, X0)
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

Kaufman, L., & Rousseeuw, P. J. (1990). Finding groups in data: an
introduction to cluster analysis. Wiley. \<doi:10.1002/9780470316801\>
Guo, T., et al. (2019). Optimal designs for genomic selection in hybrid
crops. Molecular Plant, 12, 390-401. \<doi:10.1016/j.molp.2018.12.022\>
Rio, S., et al. (2021). Genomic prediction and training set optimization
in a structured Mediterranean oat population. Theoretical and Applied
Genetics, 134, 3595-3609. \<doi:10.1007/s00122-021-03916-w\>

## Author

Jen-Hsiang Ou

## Examples

``` r
data(geno)
pam_score(geno[1:50, ], geno[51:100, ])
#> [1] 6588.776
```
