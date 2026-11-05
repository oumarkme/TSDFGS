# D-optimality

D-optimality maximizes the information content of the training set, the
determinant of the ridge information matrix.

|  |  |
|----|----|
| Type | Optimal design |
| Direction | Larger is better |
| Uses a test set | no (a test set given to [`optTrain()`](https://www.oumark.net/TSDFGS/reference/optTrain.md) is ignored) |
| Score a given set | [`dopt_score()`](https://www.oumark.net/TSDFGS/reference/dopt_score.md) |
| Search with it | `optTrain(method = "DOpt")` |

## Definition

$`\mathrm{D} = \log\det(X'X + \lambda I)`$, $`\lambda = 1/p`$, to be
maximized. STPGA returns the negative, to be minimized.

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
dopt_score(geno[1:50, ])
```

    ## [1] -1654.586

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "DOpt",
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

Akdemir, D., & Isidro-Sánchez, J. (2019). Design of training populations
for selective phenotyping in genomic prediction. *Scientific Reports,
9*, 1446. <https://doi.org/10.1038/s41598-018-38081-6>

Akdemir, D., Sanchez, J. I., & Jannink, J.-L. (2015). Optimization of
genomic selection training populations with a genetic algorithm.
*Genetics Selection Evolution, 47*, 38.
<https://doi.org/10.1186/s12711-015-0116-6>

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

- **Original software:** STPGA, GPL-3 (DOPT); compared, not copied.
- **Verdict:** included.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
