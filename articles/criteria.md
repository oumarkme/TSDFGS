# Choosing a criterion

A training set criterion is a number computed from genotypes that
predicts how well a model trained on the chosen individuals predicts the
others.
[`optTrain()`](https://www.oumark.net/TSDFGS/reference/optTrain.md)
searches for the training set that optimizes it. TSDFGS v3 has 17
criteria, each with its own page: definition, original reference, the
repository or website of the original authors, and a license review. Run
[`criteria()`](https://www.oumark.net/TSDFGS/reference/criteria.md) for
the same list in R.

## Targeted and untargeted

- **Targeted**: the individuals to predict (the test set) are genotyped
  when the training set is chosen. The criterion compares the training
  set with the test set.
- **Untargeted**: only the candidate pool is known. For the criteria of
  v2 and for PEVmax, G-optimality, and Avg_GRM_MinMax, the candidates
  left out of the training set act as the test set. MSPE, Avg_GRM,
  CDmean, and CDmean(v2) use all candidates as the target.
- **Set-only criteria** (A-, D-, E-optimality, GVoverall, GVaverage,
  Avg_GRM_self) use the training set alone and ignore a test set.

In a benchmark on seven data sets of six species, targeted optimization
beat untargeted optimization, with the largest gap at low heritability
(Fernández-González et al., 2023). Akdemir and Isidro-Sánchez (2019)
reached the same conclusion.

## Criteria in TSDFGS

| Criterion | Type | Direction | Test set | `method` |
|----|----|----|----|----|
| [r-score](https://www.oumark.net/TSDFGS/articles/criterion-rscore.md) | Model-based | larger is better | yes | `rScore` |
| [PEVmean (ridge, PCs)](https://www.oumark.net/TSDFGS/articles/criterion-pev.md) | Model-based | smaller is better | yes | `PEV` |
| [CD (ridge, PCs)](https://www.oumark.net/TSDFGS/articles/criterion-cd.md) | Model-based | larger is better | yes | `CD` |
| [A-optimality](https://www.oumark.net/TSDFGS/articles/criterion-aopt.md) | Optimal design | smaller is better | no | `AOpt` |
| [D-optimality](https://www.oumark.net/TSDFGS/articles/criterion-dopt.md) | Optimal design | larger is better | no | `DOpt` |
| [E-optimality](https://www.oumark.net/TSDFGS/articles/criterion-eopt.md) | Optimal design | smaller is better | no | `EOpt` |
| [CDmax (ridge, PCs)](https://www.oumark.net/TSDFGS/articles/criterion-cdmax.md) | Model-based | larger is better | yes | `CDmax` |
| [PEVmax](https://www.oumark.net/TSDFGS/articles/criterion-pevmax.md) | Model-based | smaller is better | yes | `PEVmax` |
| [G-optimality (GOPTPEV)](https://www.oumark.net/TSDFGS/articles/criterion-gopt.md) | Optimal design | smaller is better | yes | `GOpt` |
| [MSPE of ridge regression](https://www.oumark.net/TSDFGS/articles/criterion-mspe.md) | Model-based | smaller is better | yes | `MSPE` |
| [GVoverall](https://www.oumark.net/TSDFGS/articles/criterion-gvoverall.md) | Genomic variation | larger is better | no | `GVoverall` |
| [GVaverage](https://www.oumark.net/TSDFGS/articles/criterion-gvaverage.md) | Genomic variation | larger is better | no | `GVaverage` |
| [Avg_GRM](https://www.oumark.net/TSDFGS/articles/criterion-avggrm.md) | Relationship | larger is better | yes | `AvgGRM` |
| [Avg_GRM_self](https://www.oumark.net/TSDFGS/articles/criterion-avggrmself.md) | Relationship | smaller is better | no | `AvgGRMself` |
| [Avg_GRM_MinMax](https://www.oumark.net/TSDFGS/articles/criterion-avggrmminmax.md) | Relationship | larger is better | yes | `AvgGRMminmax` |
| [CDmean](https://www.oumark.net/TSDFGS/articles/criterion-cdmean.md) | Model-based | larger is better | yes | `CDmean` |
| [CDmean(v2)](https://www.oumark.net/TSDFGS/articles/criterion-cdmean2.md) | Model-based | larger is better | yes | `CDmeanV2` |
| [Partitioning around medoids (PAM)](https://www.oumark.net/TSDFGS/articles/criterion-pam.md) | Clustering | smaller is better | yes | `PAM` |

Left out for license or scope reasons: see the [license
review](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md).

### What is tested

Only three criteria are tested carefully: the r-score, PEVmean, and CD,
the criteria of v2.4.2. For all 12 routes of
[`optTrain()`](https://www.oumark.net/TSDFGS/reference/optTrain.md)
(with or without sub-populations, with or without a test set), the
training sets, scores, and score histories equal the output recorded
from v2.4.2.

The other 15 criteria are provided as they are described in their
papers, and the code runs without errors on all routes, but we have not
checked them carefully. Several were compared with other implementations
(the [license
review](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md)
lists which), and each has a small test against a direct formula. This
is not a validation. We recommend that you evaluate these criteria
yourself on your own data before you rely on them, for example by
scoring training sets you know with the matching `*_score()` function. A
test also does not show that a criterion picks the best training set for
your trait; the papers on each criterion page discuss that.

## What the comparisons found

- CDmean was the best targeted criterion, and the lowest mean
  relationship the best untargeted one, in the benchmark of
  Fernández-González et al. (2023). A training set of 50-55% of the
  candidates (targeted) or 65-85% (untargeted) reached 95% of the
  maximum accuracy.
- CDmean and PEVmean beat clustering under strong population structure
  in oat (Rio et al., 2021), yet PEVmean was the weakest in a
  multi-parental wheat population (Zhang et al., 2021). No criterion
  wins everywhere.
- The goal matters. Criteria that maximize the mean accuracy are not the
  best for picking the top genotypes (Chen et al., 2024).
- We found no peer-reviewed deep-learning or active-learning method for
  choosing training individuals in genomic prediction. A recent review
  names it as a future direction (Alemu et al., 2024).

## References

Akdemir, D., Sanchez, J. I., & Jannink, J.-L. (2015). Optimization of
genomic selection training populations with a genetic algorithm.
*Genetics Selection Evolution, 47*, 38.
<https://doi.org/10.1186/s12711-015-0116-6>

Akdemir, D., Rio, S., & Isidro y Sánchez, J. (2021). TrainSel: an R
package for selection of training populations. *Frontiers in Genetics,
12*, 655287. <https://doi.org/10.3389/fgene.2021.655287>

Alemu, A., et al. (2024). Genomic selection in plant breeding: key
factors shaping two decades of progress. *Molecular Plant*.

Bustos-Korts, D., et al. (2016). Improvement of predictive ability by
uniform coverage of the target genetic space. *G3, 6*, 3733-3747.
<https://doi.org/10.1534/g3.116.035410>

Chen, S.-P., et al. (2024). Constructing training sets for genomic
selection to identify superior genotypes in candidate populations.
*Theoretical and Applied Genetics*.

De Beukelaer, H., et al. (2018). Core Hunter 3: flexible core subset
selection. *BMC Bioinformatics*.

Fernández-González, J., et al. (2023). A comparison of methods for
training population optimization in genomic selection. *Theoretical and
Applied Genetics*.

Gevartosky, R., et al. (2021). Enviromic-based kernels may optimize
resource allocation with multi-trait multi-environment genomic
prediction for tropical maize. *BMC Plant Biology*.

Isidro, J., et al. (2015). Training set optimization under population
structure in genomic selection. *Theoretical and Applied Genetics, 128*,
145-158. <https://doi.org/10.1007/s00122-014-2418-4>

Laloë, D., Phocas, F., & Ménissier, F. (1996). Considerations on
measures of precision and connectedness in mixed linear models of
genetic evaluation. *Genetics Selection Evolution, 28*, 359.
<https://doi.org/10.1186/1297-9686-28-4-359>

Liu, Z.-J., & Liao, C.-T. (2026). Optimizing training sets for genomic
selection to identify superior genotypes across multiple environments.
*G3, 16*(4). <https://doi.org/10.1093/g3journal/jkag031>

Ou, J.-H., & Liao, C.-T. (2019). Training set determination for genomic
selection. *Theoretical and Applied Genetics, 132*, 2781-2792.
<https://doi.org/10.1007/s00122-019-03387-0>

Rincent, R., et al. (2012). Maximizing the reliability of genomic
selection by optimizing the calibration set of reference individuals.
*Genetics*.

Rio, S., et al. (2021). Genomic prediction and training set optimization
in a structured Mediterranean oat population. *Theoretical and Applied
Genetics, 134*, 3595-3609. <https://doi.org/10.1007/s00122-021-03916-w>

Tibbs-Cortes, L. E., et al. (2022). Genomic prediction of tocochromanols
in exotic-derived maize. *The Plant Genome, 16*, e20286.
<https://doi.org/10.1002/tpg2.20286>

Zhang, W., et al. (2021). Evaluation of genomic prediction for Fusarium
head blight resistance with a multi-parental population. *Biology, 10*,
756. <https://doi.org/10.3390/biology10080756>
