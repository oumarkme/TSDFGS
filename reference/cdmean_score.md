# CDmean score (GBLUP form)

Mean coefficient of determination of the contrasts between each target
row and the target mean, for a GBLUP model with \\K = XX'/p\\, a general
mean as the only fixed effect, and penalty \\\lambda =
\sigma_e^2/\sigma_A^2\\ (default 1). Larger is better.

## Usage

``` r
cdmean_score(X, X0, lambda = NULL)
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

Rincent, R., et al. (2012). Maximizing the reliability of genomic
selection by optimizing the calibration set of reference individuals.
Genetics, 192, 715-728. \<doi:10.1534/genetics.112.141473\> Chen, S.-P.,
Sung, W.-H., & Liao, C.-T. (2024). Constructing training sets for
genomic selection to identify superior genotypes in candidate
populations. Theoretical and Applied Genetics, 137, 270.
\<doi:10.1007/s00122-024-04766-y\>

## Author

Jen-Hsiang Ou

## Examples

``` r
data(geno)
cdmean_score(geno[1:50, ], geno[51:100, ])
#> [1] 0.6039068
```
