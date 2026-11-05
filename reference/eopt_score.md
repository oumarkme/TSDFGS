# E-optimality score

Largest eigenvalue of \\(X'X + \lambda I)^{-1}\\, that is \\1 /
(\lambda\_{min}(X'X) + \lambda)\\. Smaller is better. It uses the
training set only. When X has more columns than rows, \\X'X\\ is
singular and the score is \\1/\lambda\\ for every training set.

## Usage

``` r
eopt_score(X, lambda = NULL)
```

## Arguments

- X:

  A numeric matrix of the training set (row: individual; column: marker
  or PC).

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
eopt_score(geno[1:50, ])
#> [1] 404
```
