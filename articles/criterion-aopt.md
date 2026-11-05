# A-optimality

A-optimality minimizes the average variance of the ridge coefficients.
It uses the training set alone.

|  |  |
|----|----|
| Type | Optimal design |
| Direction | Smaller is better |
| Uses a test set | no (a test set given to [`optTrain()`](https://www.oumark.net/TSDFGS/reference/optTrain.md) is ignored) |
| Score a given set | [`aopt_score()`](https://www.oumark.net/TSDFGS/reference/aopt_score.md) |
| Search with it | `optTrain(method = "AOpt")` |

## Definition

With $`\lambda = 1/p`$:
$`\mathrm{A} = \mathrm{tr}[(X'X + \lambda I)^{-1}]`$. Akdemir and
Isidro-Sánchez (2019) allow a contrast matrix; TSDFGS uses the identity.
For $`p > n`$ the $`p - n`$ zero eigenvalues each add $`1/\lambda`$.

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
aopt_score(geno[1:50, ])
```

    ## [1] 143016

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "AOpt",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1] 324 316 310 299 402 314 379  85 277 332 101  46 196 311 158 132  56 122 268
    ## [20]  73 391 369 382 246 256 110 152 399  65 386 227 136 260  70  10 328 162  44
    ## [39] 253 262 135 309 149  62 370  52 350 245  33 138

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

- **Original software:** STPGA, GPL-3 (AOPT); compared, not copied.
  STPGA reports the trace divided by $`p`$, which does not change the
  optimum.
- **Article license:** Akdemir and Isidro-Sánchez (2019) open access
  (Scientific Reports).
- **Verdict:** included.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
