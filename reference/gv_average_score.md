# GVaverage score

\\tr(K_T)\\ with \\K_T = XX'/p\\ (average genomic variation). Larger is
better. It uses the training set only.

## Usage

``` r
gv_average_score(X)
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

Chen, S.-P., Sung, W.-H., & Liao, C.-T. (2024). Constructing training
sets for genomic selection to identify superior genotypes in candidate
populations. Theoretical and Applied Genetics, 137, 270.
\<doi:10.1007/s00122-024-04766-y\>

## Author

Jen-Hsiang Ou

## Examples

``` r
data(geno)
gv_average_score(geno[1:50, ])
#> [1] 4012.101
```
