# CDmean

The mean coefficient of determination of the GBLUP contrast between each
target individual and the target mean (Rincent et al., 2012).

|  |  |
|----|----|
| Type | Model-based |
| Direction | Larger is better |
| Uses a test set | yes; without a test set, all candidates |
| Score a given set | [`cdmean_score()`](https://www.oumark.net/TSDFGS/reference/cdmean_score.md) |
| Search with it | `optTrain(method = "CDmean")` |

## Definition

With $`K = XX'/p`$, $`M = I - J/n`$ and $`W = (MK_T + \lambda I)^{-1}M`$
($`\lambda = 1`$): for each target $`i`$ with contrast
$`c_i = e_i - \mathbf{1}/n_P`$, $`\mathrm{CD}_i = v_i'Wv_i/(c_i'Kc_i)`$
with $`v_i = K_{TP}c_i`$; the score is the mean. Equal to the CDMEAN of
STPGA in our tests.

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
cdmean_score(geno[1:50, ], geno[51:100, ])
```

    ## [1] 0.6039068

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "CDmean",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1] 181 312 119 184 178 170 251 256 277   5 330 263  32  79 158 135 105 217 199
    ## [20]  73 391 287 373  89 289 370 153 394  42 237  20  44 308  92 121 247 205 389
    ## [39] 262 221  62 298 165 160  38  15 183 196 206 230

``` r

# targeted: the test set is known
fit2 <- optTrain(geno, cand = 1:300, n.train = 50, test = 301:404,
                  method = "CDmean", min.iter = 300, console = FALSE, seed = 1)
fit2$OPTtrain
```

    ##  [1] 126 129  76   1  85  55 264 293 298  43  74 217 165  89 278  67  35  44 281
    ## [20]  34  40 216 280 225 289 226 160 220  65 133 245 283 253 262 300  71 202   9
    ## [39]  13  17  31 124 260  99 240 151   4 162 268 122

The search is a random-start exchange (hill-climbing) search; see [How
the scores are
computed](https://www.oumark.net/TSDFGS/articles/methods.md). Use
`min.iter` larger than 300 on real data and check `fit$converged`.

## Original work

Rincent, R., Laloë, D., Nicolas, S., Altmann, T., Brunel, D., Revilla,
P., et al. (2012). Maximizing the reliability of genomic selection by
optimizing the calibration set of reference individuals: comparison of
methods in two diverse groups of maize inbreds (*Zea mays* L.).
*Genetics, 192*, 715-728. <https://doi.org/10.1534/genetics.112.141473>

Chen, S.-P., Sung, W.-H., & Liao, C.-T. (2024). Constructing training
sets for genomic selection to identify superior genotypes in candidate
populations. *Theoretical and Applied Genetics, 137*, 270.
<https://doi.org/10.1007/s00122-024-04766-y>

**Repository or website of the original authors:** R package STPGA by
Akdemir and colleagues:
[CRAN.R-project.org/package=STPGA](https://CRAN.R-project.org/package=STPGA)
(GPL-3). R package TSODBG of Chen et al. (2024):
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

- **Original software:** STPGA and TSODBG, GPL; compared, not copied.
- **Article license:** Rincent et al. (2012) CHORUS.
- **Verdict:** included.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
