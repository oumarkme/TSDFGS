# Partitioning around medoids (PAM)

The k-medoids cost: the training individuals act as medoids of the
target.

|  |  |
|----|----|
| Type | Clustering |
| Direction | Smaller is better |
| Uses a test set | yes |
| Score a given set | [`pam_score()`](https://www.oumark.net/TSDFGS/reference/pam_score.md) |
| Search with it | `optTrain(method = "PAM")` |

## Definition

$`\sum_i \min_{j \in T}\lVert x_{0i} - x_j\rVert`$ over target rows
$`i`$. The classical algorithm of Kaufman and Rousseeuw (1990) searches
medoids freely; here the search is the exchange search of
[`optTrain()`](https://www.oumark.net/TSDFGS/reference/optTrain.md).

## Use in TSDFGS

``` r

library(TSDFGS)
data(geno)
pam_score(geno[1:50, ], geno[51:100, ])
```

    ## [1] 6588.776

``` r

fit <- optTrain(geno, cand = 1:404, n.train = 50, method = "PAM",
                 min.iter = 300, console = FALSE, seed = 1)
fit$OPTtrain
```

    ##  [1] 156 309  74 184 108 131 381 245 277 133 119 263 122  33 325 281 100 217 393
    ## [20]  31  41 288  70  96 289 256 368 399  42 305  20 372 262  43 153  60 172 390
    ## [39] 379 221  67 298 162  76 161 375 220 199 301 139

``` r

# targeted: the test set is known
fit2 <- optTrain(geno, cand = 1:300, n.train = 50, test = 301:404,
                  method = "PAM", min.iter = 300, console = FALSE, seed = 1)
fit2$OPTtrain
```

    ##  [1] 221 129 142  34  85  67  93 264 298  43 127 217  71  89 169  72 122  44   3
    ## [20] 199  40 245  58 109  69  39  31 220 130  90 216 283  24 193 161 137 204 103
    ## [39] 263 289 183  29 157  76 178 163  64   4 156  99

The search is a random-start exchange (hill-climbing) search; see [How
the scores are
computed](https://www.oumark.net/TSDFGS/articles/methods.md). Use
`min.iter` larger than 300 on real data and check `fit$converged`.

## Original work

Kaufman, L., & Rousseeuw, P. J. (1990). *Finding groups in data: an
introduction to cluster analysis.* Wiley.
<https://doi.org/10.1002/9780470316801>

Fernández-González, J., Akdemir, D., & Isidro y Sánchez, J. (2023). A
comparison of methods for training population optimization in genomic
selection. *Theoretical and Applied Genetics, 136*, 30.
<https://doi.org/10.1007/s00122-023-04265-6>

**Repository or website of the original authors:** Algorithm in
textbooks and the R package `cluster` (GPL \>= 2).

## License review

Formulas and algorithms are ideas and methods, which copyright does not
protect. Where an original package is under a GPL license, we read it to
check our results and credit it; TSDFGS is
[GPL-3](https://github.com/oumarkme/TSDFGS/blob/main/LICENSE.md), which
is compatible with GPL-3 and GPL (\>= 3) code. This is our reading of
the public record, not legal advice; see the [license
review](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).

- **Original software:** none required; the cost is a textbook
  definition.
- **Verdict:** included.

Reviewed 2026-11-04. [All criteria and the ones left
out](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).
