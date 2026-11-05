# PEV score

This function calculate prediction error variance (PEV) score
\<doi:10.1186/s12711-015-0116-6\> by given training set and test set.

## Usage

``` r
pev_score(X, X0, lambda = NULL)
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

A number, the PEV score (smaller is better).

## Author

Jen-Hsiang Ou

## Examples

``` r
data(geno)
if (FALSE) pev_score(geno[1:50, ], geno[51:100, ]) # \dontrun{}
```
