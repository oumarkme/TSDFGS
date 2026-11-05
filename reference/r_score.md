# r-score

This function calculate r-score \<doi:10.1007/s00122-019-03387-0\> by
given training set and test set.

## Usage

``` r
r_score(X, X0, lambda = NULL)
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

A number, the r-score.

## Details

The value equals the r-score of v2.4.2 (to its single-precision accuracy
of about 1e-7). Since v3 the computation uses p x p matrices when the
number of columns p is at most the training set size n, and n x n
matrices otherwise, and the result is returned in double precision. The
term `1 - 1/n0` (n0: test set size) of v2.4.2 is kept as published in
the software; see NEWS.md.

## References

Ou, J.-H., & Liao, C.-T. (2019). Training set determination for genomic
selection. Theoretical and Applied Genetics, 132, 2781-2792.
\<doi:10.1007/s00122-019-03387-0\>

## Author

Jen-Hsiang Ou

## Examples

``` r
data(geno)
if (FALSE) r_score(geno[1:50, ], geno[51:100, ]) # \dontrun{}
```
