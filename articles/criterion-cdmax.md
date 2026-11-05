# CDmax (ridge, PCs)

CDmax is the largest coefficient of determination among the target
individuals. It is the maximum where
[CD](https://www.oumark.net/TSDFGS/articles/criterion-cd.md) takes the
mean.

|  |  |
|----|----|
| Type | Model-based |
| Direction | Larger is better |
| Uses a test set | yes |
| Score a given set | [`cdmax_score()`](https://www.oumark.net/TSDFGS/reference/cdmax_score.md) |
| Search with it | `optTrain(method = "CDmax")` |

## Definition

With $`B = (X'X + I)^{-1}`$ and test rows $`x_i`$:
$`\mathrm{CDmax} = \max_i \, x_i' B X'X B x_i / (x_i' x_i)`$. The
formula is the CDMAX2 of STPGA. Because it is a maximum over targets, it
depends on the single best-determined individual.

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
cdmax_score(geno[1:50, ], geno[51:100, ])
```

    ## [1] 4.920842e-05

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "CDmax",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1] 228 305 316 339 278 163 286 273 403 251  14 263 318 329 242 279  59 217 212
    ## [20] 233 200 395  57 241 351   1 326 394  42 224 171 133 267 308 396 334 404  25
    ## [39] 284 261  95  58 280 352 328 400  45 247 300 112

``` r

# targeted: the test set is known
fit2 <- optTrain(geno, cand = 1:300, n.train = 50, test = 301:404,
                  method = "CDmax", min.iter = 300, console = FALSE, seed = 1)
fit2$OPTtrain
```

    ##  [1] 223  92 167 187 219 259 263 285 265  41 286 217  75 267 220  56   1 249 281
    ## [20] 229 242 208   9 241 236  95 231  14  29 181 134 206  58 233  30 273 266 160
    ## [39]  96 256 176 115 204 123 159 212 150 252 271 148

The search is a random-start exchange (hill-climbing) search; see [How
the scores are
computed](https://www.oumark.net/TSDFGS/articles/methods.md). Use
`min.iter` larger than 300 on real data and check `fit$converged`.

## Original work

Laloë, D., Phocas, F., & Ménissier, F. (1996). Considerations on
measures of precision and connectedness in mixed linear models of
genetic evaluation. *Genetics Selection Evolution, 28*, 359.
<https://doi.org/10.1186/1297-9686-28-4-359>

Rincent, R., Laloë, D., Nicolas, S., Altmann, T., Brunel, D., Revilla,
P., et al. (2012). Maximizing the reliability of genomic selection by
optimizing the calibration set of reference individuals: comparison of
methods in two diverse groups of maize inbreds (*Zea mays* L.).
*Genetics, 192*, 715-728. <https://doi.org/10.1534/genetics.112.141473>

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

- **Original software:** STPGA, GPL-3 (CDMAX2); compared with the
  definition, not copied.
- **Article license:** see the [CD
  page](https://www.oumark.net/TSDFGS/articles/criterion-cd.md).
- **Verdict:** included.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
