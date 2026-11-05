# CD-score

This function calculate CD-score \<doi:10.1186/1297-9686-28-4-359\> by
given training set and test set.

## Usage

``` r
cd_score(X, X0, lambda = NULL)
```

## Arguments

- X:

  A numeric matrix. The training set genotypic information matrix can be
  given as genotype matrix (coded as -1, 0, 1) or principal component
  matrix (row: sample; column: marker or PC).

- X0:

  A numeric matrix. The test set genotypic information matrix, with the
  same columns as X.

- lambda:

  A positive number, the ridge penalty, or NULL (default) for the
  default of the criterion (r-score and PEV: 1 / number of columns; CD:
  1).

## Value

A number, the CD score (larger is better).

## Author

Jen-Hsiang Ou

## Examples

``` r
data(geno)
if (FALSE) cd_score(geno[1:50, ], geno[51:100, ]) # \dontrun{}
```
