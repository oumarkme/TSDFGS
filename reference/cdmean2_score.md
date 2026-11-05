# CDmean(v2) score

Mean of \\diag(K\_{PT} W K\_{TP}) / diag(K\_{PP})\\, the heuristic form
of CDmean (GBLUP, \\K = XX'/p\\, default \\\lambda = 1\\). Larger is
better.

## Usage

``` r
cdmean2_score(X, X0, lambda = NULL)
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

Chen, S.-P., Sung, W.-H., & Liao, C.-T. (2024). Constructing training
sets for genomic selection to identify superior genotypes in candidate
populations. Theoretical and Applied Genetics, 137, 270.
\<doi:10.1007/s00122-024-04766-y\>

## Author

Jen-Hsiang Ou

## Examples

``` r
data(geno)
cdmean2_score(geno[1:50, ], geno[51:100, ])
#> [1] 0.6036041
```
