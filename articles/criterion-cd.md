# CD (ridge, PCs)

CD is the squared correlation between true and predicted genetic value
(reliability). The score is the mean CD of the test individuals,
computed on PCs.

|  |  |
|----|----|
| Type | Model-based |
| Direction | Larger is better |
| Uses a test set | yes |
| Score a given set | [`cd_score()`](https://www.oumark.net/TSDFGS/reference/cd_score.md) |
| Search with it | `optTrain(method = "CD")` |

## Definition

With $`B = (X'X + I)^{-1}`$ and $`x_i`$ a test row:

``` math
\mathrm{CD} = \frac{1}{n_0}\sum_i \frac{x_i' B X'X B x_i}{x_i' x_i}.
```

This equals the CDMEAN2 of STPGA in our tests. For the GBLUP form with
contrasts, see
[CDmean](https://www.oumark.net/TSDFGS/articles/criterion-cdmean.md).

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
cd_score(geno[1:50, ], geno[51:100, ])
```

    ## [1] 1.05567e-05

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "CD",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1] 230 193 129 183 325  73 307  85 245 345 101 263 349 298 261  37 301  21 199
    ## [20] 215 304 287 393  96 289 371 214 375  42 306 178 143 342 302 155 357 260  48
    ## [39] 390  13  92 133  81  77 156 130 352 137 281 246

``` r

# targeted: the test set is known
fit2 <- optTrain(geno, cand = 1:300, n.train = 50, test = 301:404,
                  method = "CD", min.iter = 300, console = FALSE, seed = 1)
fit2$OPTtrain
```

    ##  [1] 254  32 163 131  85 252 263 138 273 261 128 217  67 262 151 276 282  10 207
    ## [20]  97  39 259 133 241 272   1 147  14  28 127 204 271  55 193 194 179 206 250
    ## [39] 268 289 123  34 157 234 239 144  61 221 156 117

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

- **Original software:** STPGA, GPL-3; compared, not copied.
- **Article license:** Laloë et al. (1996) Genetics Selection Evolution;
  Rincent et al. (2012) CHORUS.
- **Verdict:** included.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
