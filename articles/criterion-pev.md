# PEVmean (ridge, PCs)

The prediction error variance (PEV) measures how uncertain the model
predictions are for the test set. TSDFGS sums it over the test
individuals.

|  |  |
|----|----|
| Type | Model-based |
| Direction | Smaller is better |
| Uses a test set | yes |
| Score a given set | [`pev_score()`](https://www.oumark.net/TSDFGS/reference/pev_score.md) |
| Search with it | `optTrain(method = "PEV")` |

## Definition

With an intercept column added, $`B = (X'X + \lambda I)^{-1}`$ and
$`\lambda = 1/p`$:

``` math
\mathrm{PEV} = \mathrm{tr}(X_0 B X_0').
```

Dividing by $`n_0`$ gives PEVmean. This is the PC-ridge form of Akdemir
et al. (2015). The relationship-matrix PEV of Rincent et al. (2012) uses
contrasts instead; see
[CDmean](https://www.oumark.net/TSDFGS/articles/criterion-cdmean.md).

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
pev_score(geno[1:50, ], geno[51:100, ])
```

    ## [1] 257145357

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "PEV",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1]   3 324 119 299 379 129 149  85 277 254 330 263  32 162 158 133  97 217 392
    ## [20] 312 391 287 142  89  65 110 329 372  42 308  20 190 199 311  28 376 172  15
    ## [39] 138 221  39 298  82 160  38 390 185 256 304 127

``` r

# targeted: the test set is known
fit2 <- optTrain(geno, cand = 1:300, n.train = 50, test = 301:404,
                  method = "PEV", min.iter = 300, console = FALSE, seed = 1)
fit2$OPTtrain
```

    ##  [1] 183 218  17 129 171 112  41  93  75  71 261  16 143 162 146  76 298 122 289
    ## [20] 163 272  59  13 180 164 190 245 199 157 216 160  69   2  97 130 152 113  47
    ## [39]  44  35  31 151 165 300  65 264  43  89 221 240

The search is a random-start exchange (hill-climbing) search; see [How
the scores are
computed](https://www.oumark.net/TSDFGS/articles/methods.md). Use
`min.iter` larger than 300 on real data and check `fit$converged`.

## Original work

Akdemir, D., Sanchez, J. I., & Jannink, J.-L. (2015). Optimization of
genomic selection training populations with a genetic algorithm.
*Genetics Selection Evolution, 47*, 38.
<https://doi.org/10.1186/s12711-015-0116-6>

Rincent, R., Laloë, D., Nicolas, S., Altmann, T., Brunel, D., Revilla,
P., et al. (2012). Maximizing the reliability of genomic selection by
optimizing the calibration set of reference individuals: comparison of
methods in two diverse groups of maize inbreds (*Zea mays* L.).
*Genetics, 192*, 715-728. <https://doi.org/10.1534/genetics.112.141473>

**Repository or website of the original authors:** R package STPGA by
Akdemir and colleagues:
[CRAN.R-project.org/package=STPGA](https://CRAN.R-project.org/package=STPGA)
(GPL-3). Rincent et al. (2012) give no software link.

## License review

Formulas and algorithms are ideas and methods, which copyright does not
protect. Where an original package is under a GPL license, we read it to
check our results and credit it; TSDFGS is
[GPL-3](https://github.com/oumarkme/TSDFGS/blob/main/LICENSE.md), which
is compatible with GPL-3 and GPL (\>= 3) code. This is our reading of
the public record, not legal advice; see the [license
review](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).

- **Original software:** STPGA, GPL-3; we compared outputs and did not
  copy code.
- **Article license:** Akdemir et al. (2015) CC BY 4.0; Rincent et
  al. (2012) CHORUS public access.
- **Patents:** none in the articles. TrainSel (2021) says its concepts
  may be patented; PEV predates it and TSDFGS does not use TrainSel.
- **Verdict:** included.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
