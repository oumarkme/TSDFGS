# r-score

The r-score is the expected Pearson correlation between predicted and
observed values of the test set, written in closed form for ridge
regression. It was introduced with this package.

|  |  |
|----|----|
| Type | Model-based |
| Direction | Larger is better |
| Uses a test set | yes |
| Score a given set | [`r_score()`](https://www.oumark.net/TSDFGS/reference/r_score.md) |
| Search with it | `optTrain(method = "rScore")` |

## Definition

With training matrix $`X`$ ($`n \times p`$, PCs or markers), test matrix
$`X_0`$ with $`n_0`$ rows, and penalty $`\lambda = 1/p`$, the r-score is
$`q_{12}/\sqrt{q_1 q_2}`$: the covariance of predicted and true test
values over the product of their variances. Ou and Liao (2019) derive
the three terms. TSDFGS evaluates them with a $`p \times p`$ matrix when
$`p \le n`$ and an $`n \times n`$ matrix otherwise.

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
r_score(geno[1:50, ], geno[51:100, ])
```

    ## [1] 0.7850248

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "rScore",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1]   3 248 380 281 282 311 198  85 277 333 330 263  32 308 157 370 105  23 392
    ## [20] 389 390 368  73 245  65   1 323 394  42 301  20 134 166 185 278 237 158  15
    ## [39] 288 221 122 298 345 120 153 375 184  97 186 304

``` r

# targeted: the test set is known
fit2 <- optTrain(geno, cand = 1:300, n.train = 50, test = 301:404,
                  method = "rScore", min.iter = 300, console = FALSE, seed = 1)
fit2$OPTtrain
```

    ##  [1] 164  76 179 178 300  71 189 289 162 123  43  44  21 216 267 285 225 190 160
    ## [20] 287  67  69  64 151  47  89 196  85 172  31  65 246 221  39 130 214 157 129
    ## [39] 199 110  97 298  86 244  17 165  34  59 207 241

The search is a random-start exchange (hill-climbing) search; see [How
the scores are
computed](https://www.oumark.net/TSDFGS/articles/methods.md). Use
`min.iter` larger than 300 on real data and check `fit$converged`.

## Original work

Ou, J.-H., & Liao, C.-T. (2019). Training set determination for genomic
selection. *Theoretical and Applied Genetics, 132*, 2781-2792.
<https://doi.org/10.1007/s00122-019-03387-0>

**Repository or website of the original authors:**
[github.com/oumarkme/TSDFGS](https://github.com/oumarkme/TSDFGS), the
package of the method’s authors.

## License review

Formulas and algorithms are ideas and methods, which copyright does not
protect. Where an original package is under a GPL license, we read it to
check our results and credit it; TSDFGS is
[GPL-3](https://github.com/oumarkme/TSDFGS/blob/main/LICENSE.md), which
is compatible with GPL-3 and GPL (\>= 3) code. This is our reading of
the public record, not legal advice; see the [license
review](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).

- **Method and software:** Ou and Liao (2019); TSDFGS is GPL (\>= 3), as
  v2.4.2 was.
- **Article license:** Springer, not open access; cited, not reproduced.
- **Verdict:** included.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
