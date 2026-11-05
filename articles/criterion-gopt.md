# G-optimality (GOPTPEV)

G-optimality minimizes the largest eigenvalue of the PEV matrix of the
target, the variance in the worst direction.

|  |  |
|----|----|
| Type | Optimal design |
| Direction | Smaller is better |
| Uses a test set | yes |
| Score a given set | [`gopt_score()`](https://www.oumark.net/TSDFGS/reference/gopt_score.md) |
| Search with it | `optTrain(method = "GOpt")` |

## Definition

With $`B`$ as for
[PEVmax](https://www.oumark.net/TSDFGS/articles/criterion-pevmax.md):
$`\mathrm{G} = \lambda_{\max}(X_0 B X_0')`$, computed from the
$`p \times p`$ matrix $`L^{-1}X_0'X_0L^{-T}`$ with
$`LL' = X'X + \lambda I`$.

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
gopt_score(geno[1:50, ], geno[51:100, ])
```

    ## [1] 20127415

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "GOpt",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1] 391 347 173 336 127 315 146  62 342 393 330 290 134 153 158 369 297 217 392
    ## [20] 389 390 128 373 244 350 364 222 162  69 304  20 283 343 300 151 376  81  15
    ## [39] 229  44 135 298 149 160  36 388  45 250 281  80

``` r

# targeted: the test set is known
fit2 <- optTrain(geno, cand = 1:300, n.train = 50, test = 301:404,
                  method = "GOpt", min.iter = 300, console = FALSE, seed = 1)
fit2$OPTtrain
```

    ##  [1] 256 126 143 137 140  93  97  84 298 129 196 278  31  89 259 152 127  44 168
    ## [20]  76 108 212  71 117 198 216 162 232  67  39 245  74  10 220  72  24 300  52
    ## [39] 151 150  65 225 260  32 141 116 230 249  43  77

The search is a random-start exchange (hill-climbing) search; see [How
the scores are
computed](https://www.oumark.net/TSDFGS/articles/methods.md). Use
`min.iter` larger than 300 on real data and check `fit$converged`.

## Original work

Akdemir, D., Sanchez, J. I., & Jannink, J.-L. (2015). Optimization of
genomic selection training populations with a genetic algorithm.
*Genetics Selection Evolution, 47*, 38.
<https://doi.org/10.1186/s12711-015-0116-6>

Isidro y Sánchez, J., & Akdemir, D. (2021). Training set optimization
for sparse phenotyping in genomic selection: a conceptual overview.
*Frontiers in Plant Science, 12*, 715910.
<https://doi.org/10.3389/fpls.2021.715910>

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

- **Original software:** STPGA, GPL-3 (GOPTPEV); compared, not copied.
- **Verdict:** included.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
