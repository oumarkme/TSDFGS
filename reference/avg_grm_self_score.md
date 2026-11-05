# Avg_GRM_self score

Mean of \\G = XX'/p\\ within the training set. Smaller is better. It
uses the training set only.

## Usage

``` r
avg_grm_self_score(X)
```

## Arguments

- X:

  A numeric matrix of the training set (row: individual; column: marker
  or PC).

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
avg_grm_self_score(geno[1:50, ])
#> [1] 1.412919
```
