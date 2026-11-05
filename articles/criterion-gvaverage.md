# GVaverage

GVaverage is the mean genetic variance of the training individuals, the
trace of their relationship matrix divided by $`n`$ (TSDFGS returns the
trace).

|  |  |
|----|----|
| Type | Genomic variation |
| Direction | Larger is better |
| Uses a test set | no (a test set given to [`optTrain()`](https://www.oumark.net/TSDFGS/reference/optTrain.md) is ignored) |
| Score a given set | [`gv_average_score()`](https://www.oumark.net/TSDFGS/reference/gv_average_score.md) |
| Search with it | `optTrain(method = "GVaverage")` |

## Definition

$`\mathrm{tr}(G) = \lVert X\rVert_F^2/p`$. Maximizing it favors
individuals far from the population center.

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
gv_average_score(geno[1:50, ])
```

    ## [1] 4012.101

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "GVaverage",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1] 403 348 311 136 339 313 379  65  33 185 399  51 329 301 244  69 179 122 391
    ## [20] 335 350  97 382  18 298 323 216  12  27 296  20  44 312 305  10  40 189  52
    ## [39] 376  16  67 128 162  77  36 389 347 245  32  82

The search is a random-start exchange (hill-climbing) search; see [How
the scores are
computed](https://www.oumark.net/TSDFGS/articles/methods.md). Use
`min.iter` larger than 300 on real data and check `fit$converged`.

## Original work

Chen, S.-P., Sung, W.-H., & Liao, C.-T. (2024). Constructing training
sets for genomic selection to identify superior genotypes in candidate
populations. *Theoretical and Applied Genetics, 137*, 270.
<https://doi.org/10.1007/s00122-024-04766-y>

Fernández-González, J., Akdemir, D., & Isidro y Sánchez, J. (2023). A
comparison of methods for training population optimization in genomic
selection. *Theoretical and Applied Genetics, 136*, 30.
<https://doi.org/10.1007/s00122-023-04265-6>

**Repository or website of the original authors:** R package TSODBG of
Chen et al. (2024):
[github.com/spcspin/TSODBG](https://github.com/spcspin/TSODBG) (GPL \>=
3).

## License review

Formulas and algorithms are ideas and methods, which copyright does not
protect. Where an original package is under a GPL license, we read it to
check our results and credit it; TSDFGS is
[GPL-3](https://github.com/oumarkme/TSDFGS/blob/main/LICENSE.md), which
is compatible with GPL-3 and GPL (\>= 3) code. This is our reading of
the public record, not legal advice; see the [license
review](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).

- **Original software:** TSODBG, GPL (\>= 3); compared, not copied.
- **Verdict:** included.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
