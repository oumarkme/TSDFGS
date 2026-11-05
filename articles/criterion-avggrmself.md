# Avg_GRM_self

The mean relationship within the training set; minimizing it spreads the
training individuals apart.

|  |  |
|----|----|
| Type | Relationship |
| Direction | Smaller is better |
| Uses a test set | no (a test set given to [`optTrain()`](https://www.oumark.net/TSDFGS/reference/optTrain.md) is ignored) |
| Score a given set | [`avg_grm_self_score()`](https://www.oumark.net/TSDFGS/reference/avg_grm_self_score.md) |
| Search with it | `optTrain(method = "AvgGRMself")` |

## Definition

$`\mathbf{1}_n' G_{TT}\mathbf{1}_n/n^2`$ with $`G = XX'/p`$.

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
avg_grm_self_score(geno[1:50, ])
```

    ## [1] 1.412919

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "AvgGRMself",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1] 123 351 129  97 137 187 350  85 143 189 119 104 389 176 160 346 272 217 393
    ## [20] 127 290 288 382  89 289   1 360 331  42 125  20  44 175  70  30  40 163  25
    ## [39] 259  17  39 298 111  76  14 130  45 196 281 230

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
