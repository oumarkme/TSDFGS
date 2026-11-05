# GVoverall

GVoverall measures the total genetic variation covered by the training
set: the generalized variance of its relationship matrix.

|  |  |
|----|----|
| Type | Genomic variation |
| Direction | Larger is better |
| Uses a test set | no (a test set given to [`optTrain()`](https://www.oumark.net/TSDFGS/reference/optTrain.md) is ignored) |
| Score a given set | [`gv_overall_score()`](https://www.oumark.net/TSDFGS/reference/gv_overall_score.md) |
| Search with it | `optTrain(method = "GVoverall")` |

## Definition

With $`G = XX'/p`$: $`\log\det G`$ (the log of the determinant in the
benchmark of Fernández-González et al., 2023). It needs $`p \ge n`$; use
all PCs or markers.

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
gv_overall_score(geno[1:50, ])
```

    ## [1] 169.8436

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "GVoverall",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1] 324 244 314 299 402 315 321 253 216 185 398  49 103 328 262 278 296  20 392
    ## [20] 388 297 168 196  27 282 260 325 399 369 227 131 128 197 245  10 312  77 389
    ## [39] 256  16 135 268 149  60 152  62 350 316 298 138

The search is a random-start exchange (hill-climbing) search; see [How
the scores are
computed](https://www.oumark.net/TSDFGS/articles/methods.md). Use
`min.iter` larger than 300 on real data and check `fit$converged`.

## Original work

Fernández-González, J., Akdemir, D., & Isidro y Sánchez, J. (2023). A
comparison of methods for training population optimization in genomic
selection. *Theoretical and Applied Genetics, 136*, 30.
<https://doi.org/10.1007/s00122-023-04265-6>

Chen, S.-P., Sung, W.-H., & Liao, C.-T. (2024). Constructing training
sets for genomic selection to identify superior genotypes in candidate
populations. *Theoretical and Applied Genetics, 137*, 270.
<https://doi.org/10.1007/s00122-024-04766-y>

**Repository or website of the original authors:** Code of
Fernández-González et al. (2023):
[github.com/TheRocinante-lab/Publications](https://github.com/TheRocinante-lab/Publications/tree/main/2022/Fernandez-Gonzalez_et_al_2022_Comparison).
Package TrainSel:
[github.com/TheRocinante-lab/TrainSel](https://github.com/TheRocinante-lab/TrainSel).
R package TSODBG of Chen et al. (2024):
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

- **Original software:** the repository of Fernández-González et
  al. (2023) has no LICENSE file, and TrainSel carries an
  all-rights-reserved license with a patent notice. We did not read or
  copy either; the definition comes from the open-access article.
- **TSODBG:** GPL (\>= 3); compared, not copied.
- **Article license:** Fernández-González et al. (2023) CC BY 4.0.
- **Verdict:** included, written from the article.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
