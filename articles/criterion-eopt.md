# E-optimality

E-optimality limits the worst-estimated direction: the inverse of the
smallest eigenvalue of the information matrix.

|  |  |
|----|----|
| Type | Optimal design |
| Direction | Smaller is better |
| Uses a test set | no (a test set given to [`optTrain()`](https://www.oumark.net/TSDFGS/reference/optTrain.md) is ignored) |
| Score a given set | [`eopt_score()`](https://www.oumark.net/TSDFGS/reference/eopt_score.md) |
| Search with it | `optTrain(method = "EOpt")` |

## Definition

$`\mathrm{E} = 1/(\lambda_{\min}(X'X) + \lambda)`$. When $`p > n`$ the
smallest eigenvalue is zero and the score is the constant $`1/\lambda`$,
so the criterion does not discriminate between training sets.

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
eopt_score(geno[1:50, ])
```

    ## [1] 404

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "EOpt",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1] 324 167 129 299 270 187 307  85 277 362 330 263 329  79 213  37 105 217 366
    ## [20] 165 290 395 382  89 289 340 326 394  42 111  20  44 343  70 121  40 172  25
    ## [39] 248 198  39 298 280 160  14 130  45  22 206 230

The search is a random-start exchange (hill-climbing) search; see [How
the scores are
computed](https://www.oumark.net/TSDFGS/articles/methods.md). Use
`min.iter` larger than 300 on real data and check `fit$converged`.

## Original work

Akdemir, D., & Isidro-Sánchez, J. (2019). Design of training populations
for selective phenotyping in genomic prediction. *Scientific Reports,
9*, 1446. <https://doi.org/10.1038/s41598-018-38081-6>

**Repository or website of the original authors:** R package STPGA by
Akdemir and colleagues:
[CRAN.R-project.org/package=STPGA](https://CRAN.R-project.org/package=STPGA)
(GPL-3).

## License review

Formulas and algorithms are ideas and methods, which copyright does not
protect. Where an original package is under a GPL license, we read it to
check our results and credit it; TSDFGS is
[GPL-3](https://github.com/oumarkme/TSDFGS/blob/main/LICENSE.md), which
is compatible with GPL-3 and GPL (\>= 3) code. This is our reading of
the public record, not legal advice; see the [license
review](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).

- **Original software:** STPGA, GPL-3 (EOPT); compared, not copied.
- **Verdict:** included.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
