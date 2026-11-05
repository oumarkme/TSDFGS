# A-optimality score

\\tr\[(X'X + \lambda I)^{-1}\]\\: the trace of the inverse ridge
information matrix of the training set. Smaller is better. It uses the
training set only.

## Usage

``` r
aopt_score(X, lambda = NULL)
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
aopt_score(geno[1:50, ])
#> [1] 143016
```
