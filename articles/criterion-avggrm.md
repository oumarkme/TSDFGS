# Avg_GRM

The mean genomic relationship between the training set and the target.

|  |  |
|----|----|
| Type | Relationship |
| Direction | Larger is better |
| Uses a test set | yes; without a test set, all candidates |
| Score a given set | [`avg_grm_score()`](https://www.oumark.net/TSDFGS/reference/avg_grm_score.md) |
| Search with it | `optTrain(method = "AvgGRM")` |

## Definition

$`\mathrm{Avg\_GRM} = \mathbf{1}_n' G_{TP} \mathbf{1}_{n_P}/(n\,n_P)`$
with $`G = XX'/p`$. Without a test set, the target is the whole
candidate set.

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
avg_grm_score(geno[1:50, ], geno[51:100, ])
```

    ## [1] -0.158039

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "AvgGRM",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1] 403 303  74 339 142  71 307  85 277 304 330  51  12 244 321  60 117 240 392
    ## [20] 306 302 168 216 325 300  96 154  43  65 301  20  44 320 305 121 335 260 388
    ## [39]  82   3 152 298  16  77  36  17 349 246 348  97

``` r

# targeted: the test set is known
fit2 <- optTrain(geno, cand = 1:300, n.train = 50, test = 301:404,
                  method = "AvgGRM", min.iter = 300, console = FALSE, seed = 1)
fit2$OPTtrain
```

    ##  [1]  77 245 143  43 215  31  88 149 116 297 185 268  52 161  27  65 283  44 168
    ## [20] 121  40 244  10 282 296  39 162  20  12 246 142 195  71  18  74 256 216 299
    ## [39]  67  17  96 189 117 126 248 277   5 221 156 101

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
- **Article license:** CC BY 4.0.
- **Verdict:** included, written from the article.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
