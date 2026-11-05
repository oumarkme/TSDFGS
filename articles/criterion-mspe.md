# MSPE of ridge regression

The expected mean squared prediction error of ridge regression for the
target individuals (Chen et al., 2024). Without a test set, all
candidates are the target.

|  |  |
|----|----|
| Type | Model-based |
| Direction | Smaller is better |
| Uses a test set | yes; without a test set, all candidates |
| Score a given set | [`mspe_score()`](https://www.oumark.net/TSDFGS/reference/mspe_score.md) |
| Search with it | `optTrain(method = "MSPE")` |

## Definition

With $`A = X'(XX' + I)^{-1}`$ and target rows $`X_c`$ ($`n_c`$ rows):
$`\mathrm{MSPE} = 1 + \mathrm{tr}(X_c A A' X_c')/n_c + \lVert X_c - X_c A X\rVert_F^2`$
(Eq. 4 of Chen et al., 2024).

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
mspe_score(geno[1:50, ], geno[51:100, ])
```

    ## [1] 636411

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "MSPE",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1]   3 324 119 299 379 129 149  85 277 254 330 263  32 162 158 133  97 217 392
    ## [20] 312 391 287 142  89  65 110 329 372  42 308  20 190 199 311  28 376 172  15
    ## [39] 138 221  39 298  82 160  38 390 185 256 304 127

``` r

# targeted: the test set is known
fit2 <- optTrain(geno, cand = 1:300, n.train = 50, test = 301:404,
                  method = "MSPE", min.iter = 300, console = FALSE, seed = 1)
fit2$OPTtrain
```

    ##  [1]  52 129  76 196  85 260  99 148 298  43  74 217 165  89 157  66 127  44 216
    ## [20]  31  40 172 108 162 193 116 160  71  72  97 245 283  65 183 164 254 300  16
    ## [39]  58  17 189  34 202 204 173 212 268 199  78 136

The search is a random-start exchange (hill-climbing) search; see [How
the scores are
computed](https://www.oumark.net/TSDFGS/articles/methods.md). Use
`min.iter` larger than 300 on real data and check `fit$converged`.

## Original work

Chen, S.-P., Sung, W.-H., & Liao, C.-T. (2024). Constructing training
sets for genomic selection to identify superior genotypes in candidate
populations. *Theoretical and Applied Genetics, 137*, 270.
<https://doi.org/10.1007/s00122-024-04766-y>

**Repository or website of the original authors:** R package TSODBG of
Chen et al. (2024):
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

- **Original software:** TSODBG, GPL (\>= 3); we compared with its
  `mspe_score` and did not copy.
- **Article license:** CC BY-NC-ND 4.0; cited, not reproduced.
- **Verdict:** included.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
