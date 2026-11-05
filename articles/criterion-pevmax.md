# PEVmax

PEVmax is the largest prediction error variance among the test
individuals, so it protects the worst-predicted one.

|  |  |
|----|----|
| Type | Model-based |
| Direction | Smaller is better |
| Uses a test set | yes |
| Score a given set | [`pevmax_score()`](https://www.oumark.net/TSDFGS/reference/pevmax_score.md) |
| Search with it | `optTrain(method = "PEVmax")` |

## Definition

$`\mathrm{PEVmax} = \max_i x_{0i}' B x_{0i}`$ with an intercept column
and $`B = (X'X + \lambda I)^{-1}`$. Without a test set the candidates
left out of the training set are the target.

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
pevmax_score(geno[1:50, ], geno[51:100, ])
```

    ## [1] 12487101

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "PEVmax",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1] 324 310 171 218 283  73 150  81 226 354 337  88 314 313 165 372 180 243 220
    ## [20] 339 227 373  68 335 292 239 262  78 371 247 268 189 298 312  70 194  79 208
    ## [39]  75  11  56 118  42  82 156 321  30 237  33 230

``` r

# targeted: the test set is known
fit2 <- optTrain(geno, cand = 1:300, n.train = 50, test = 301:404,
                  method = "PEVmax", min.iter = 300, console = FALSE, seed = 1)
fit2$OPTtrain
```

    ##  [1]  78 129  96  28 216 277  35 293 179  43 186  27 152 296 161  68 103 142  70
    ## [20]  32 137 110 113 149 290 201 164  14  80  40 189 283 260 262 101  74 202  38
    ## [39]  61  18  33 230 157 246 141  31 197 256 225 108

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

- **Original software:** STPGA, GPL-3 (PEVMAX); compared, not copied.
- **Article license:** Akdemir et al. (2015) CC BY 4.0; Isidro y Sánchez
  and Akdemir (2021) open access.
- **Verdict:** included.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
