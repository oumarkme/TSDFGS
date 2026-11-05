# MSPE score of ridge regression

The approximation of the expected mean squared prediction error of ridge
regression, \\1 + tr\[X_0 A A' X_0'\]/n_0 + tr\[(X_0 - X_0 A X)'(X_0 -
X_0 A X)\]\\ with \\A = X'(XX' + \lambda I)^{-1}\\. Smaller is better.
The default penalty is \\\lambda = 1\\.

## Usage

``` r
mspe_score(X, X0, lambda = NULL)
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
mspe_score(geno[1:50, ], geno[51:100, ])
#> [1] 636411
```
