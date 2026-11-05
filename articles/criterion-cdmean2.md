# CDmean(v2)

The heuristic form of CDmean from Chen et al. (2024): the contrast is
replaced by the individual itself.

|  |  |
|----|----|
| Type | Model-based |
| Direction | Larger is better |
| Uses a test set | yes; without a test set, all candidates |
| Score a given set | [`cdmean2_score()`](https://www.oumark.net/TSDFGS/reference/cdmean2_score.md) |
| Search with it | `optTrain(method = "CDmeanV2")` |

## Definition

Mean of $`\mathrm{diag}(K_{PT}WK_{TP})/\mathrm{diag}(K_{PP})`$. TSODBG
sums instead of averaging; the optimum is the same for a fixed target.

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
cdmean2_score(geno[1:50, ], geno[51:100, ])
```

    ## [1] 0.6036041

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "CDmeanV2",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1] 181 312 119 184 178 170 251 256 277   5 330 263  32  79 158 135 105 217 199
    ## [20]  73 391 287 373  89 289 370 153 394  42 237  20  44 308  92 121 247 205 389
    ## [39] 262 221  62 298 165 160  38  15 183 196 206 230

``` r

# targeted: the test set is known
fit2 <- optTrain(geno, cand = 1:300, n.train = 50, test = 301:404,
                  method = "CDmeanV2", min.iter = 300, console = FALSE, seed = 1)
fit2$OPTtrain
```

    ##  [1] 162 129 102 183  85 293  65 147 298  43 163 217 165  89 261  34 199  44 168
    ## [20]  76  99 172  71 186 220 216 160  12 130 133 116 118  60 193  70  67 300  16
    ## [39] 151 289  93  31  15  48 245 252 180 221 247  97

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

- **Original software:** TSODBG, GPL (\>= 3); compared, not copied.
- **Article license:** CC BY-NC-ND 4.0.
- **Verdict:** included.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
