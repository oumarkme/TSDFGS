# Avg_GRM_MinMax

Combines the two: high relationship to the target and low relationship
within the training set.

|  |  |
|----|----|
| Type | Relationship |
| Direction | Larger is better |
| Uses a test set | yes |
| Score a given set | [`avg_grm_minmax_score()`](https://www.oumark.net/TSDFGS/reference/avg_grm_minmax_score.md) |
| Search with it | `optTrain(method = "AvgGRMminmax")` |

## Definition

$`\mathrm{Avg\_GRM}(T,P) - \mathrm{Avg\_GRM\_self}(T)`$. Without a test
set the target is the candidates left out of the training set.

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
avg_grm_minmax_score(geno[1:50, ], geno[51:100, ])
```

    ## [1] -1.570958

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "AvgGRMminmax",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1] 123 351 129  97 137 187 350  85 143 189 119 104 389 176 160 346 272 217 393
    ## [20] 127 290 288 382  89 289   1 360 331  42 125  20  44 175  70  30  40 163  25
    ## [39] 259  17  39 298 111  76  14 130  45 196 281 230

``` r

# targeted: the test set is known
fit2 <- optTrain(geno, cand = 1:300, n.train = 50, test = 301:404,
                  method = "AvgGRMminmax", min.iter = 300, console = FALSE, seed = 1)
fit2$OPTtrain
```

    ##  [1]  82 129 142 209  85 143 164 147 298 276 287 217  76  89  42 111 268  44 216
    ## [20] 282  40 172  25 116 195  39 160 213 135  45 246 283  60  18 109 254 294 258
    ## [39]  57  17  97 189  84 265 141 252  67 251  48 177

The search is a random-start exchange (hill-climbing) search; see [How
the scores are
computed](https://www.oumark.net/TSDFGS/articles/methods.md). Use
`min.iter` larger than 300 on real data and check `fit$converged`.

## Original work

Fernández-González, J., Akdemir, D., & Isidro y Sánchez, J. (2023). A
comparison of methods for training population optimization in genomic
selection. *Theoretical and Applied Genetics, 136*, 30.
<https://doi.org/10.1007/s00122-023-04265-6>

**Repository or website of the original authors:** Code of
Fernández-González et al. (2023):
[github.com/TheRocinante-lab/Publications](https://github.com/TheRocinante-lab/Publications/tree/main/2022/Fernandez-Gonzalez_et_al_2022_Comparison).
Package TrainSel:
[github.com/TheRocinante-lab/TrainSel](https://github.com/TheRocinante-lab/TrainSel).

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
- **Verdict:** included, written from the article.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
