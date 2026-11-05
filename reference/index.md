# Package index

## Find the training set

The main function. It searches the candidates for the training set with
the best criterion value.

- [`optTrain()`](https://www.oumark.net/TSDFGS/reference/optTrain.md) :
  Optimal training set determination

## Find the training set size

Simulate r-scores for several sizes and fit the operating curve.

- [`SSDFGS()`](https://www.oumark.net/TSDFGS/reference/SSDFGS.md) :
  Sample size determination for genomic selection
- [`FGCM()`](https://www.oumark.net/TSDFGS/reference/FGCM.md) : Fit
  logistic growth curve model
- [`nt2r()`](https://www.oumark.net/TSDFGS/reference/nt2r.md) : Simulate
  r-scores of each training set size

## Criteria: list

All criteria, with the name to give to optTrain(method = ).

- [`criteria()`](https://www.oumark.net/TSDFGS/reference/criteria.md) :
  List of training set criteria

## Score a given training set

Criterion values of a training set. One page per criterion under
Criteria in the top bar.

- [`r_score()`](https://www.oumark.net/TSDFGS/reference/r_score.md) :
  r-score
- [`pev_score()`](https://www.oumark.net/TSDFGS/reference/pev_score.md)
  : PEV score
- [`cd_score()`](https://www.oumark.net/TSDFGS/reference/cd_score.md) :
  CD-score
- [`pevmax_score()`](https://www.oumark.net/TSDFGS/reference/pevmax_score.md)
  : PEVmax score
- [`mspe_score()`](https://www.oumark.net/TSDFGS/reference/mspe_score.md)
  : MSPE score of ridge regression
- [`cdmean_score()`](https://www.oumark.net/TSDFGS/reference/cdmean_score.md)
  : CDmean score (GBLUP form)
- [`cdmean2_score()`](https://www.oumark.net/TSDFGS/reference/cdmean2_score.md)
  : CDmean(v2) score
- [`cdmax_score()`](https://www.oumark.net/TSDFGS/reference/cdmax_score.md)
  : CDmax score
- [`pam_score()`](https://www.oumark.net/TSDFGS/reference/pam_score.md)
  : PAM (k-medoids) cost

## Optimal design criteria

- [`aopt_score()`](https://www.oumark.net/TSDFGS/reference/aopt_score.md)
  : A-optimality score
- [`dopt_score()`](https://www.oumark.net/TSDFGS/reference/dopt_score.md)
  : D-optimality score
- [`eopt_score()`](https://www.oumark.net/TSDFGS/reference/eopt_score.md)
  : E-optimality score
- [`gopt_score()`](https://www.oumark.net/TSDFGS/reference/gopt_score.md)
  : G-optimality score (GOPTPEV)

## Genomic variation and relationship criteria

- [`gv_overall_score()`](https://www.oumark.net/TSDFGS/reference/gv_overall_score.md)
  : GVoverall score
- [`gv_average_score()`](https://www.oumark.net/TSDFGS/reference/gv_average_score.md)
  : GVaverage score
- [`avg_grm_score()`](https://www.oumark.net/TSDFGS/reference/avg_grm_score.md)
  : Avg_GRM score
- [`avg_grm_self_score()`](https://www.oumark.net/TSDFGS/reference/avg_grm_self_score.md)
  : Avg_GRM_self score
- [`avg_grm_minmax_score()`](https://www.oumark.net/TSDFGS/reference/avg_grm_minmax_score.md)
  : Avg_GRM_MinMax score

## Example data

Rice genotype principal components (Zhao et al., 2011) and
sub-population labels.

- [`geno`](https://www.oumark.net/TSDFGS/reference/geno.md) : Genotype
  information
- [`subpop`](https://www.oumark.net/TSDFGS/reference/subpop.md) :
  Sub-population information

## Package

- [`TSDFGS`](https://www.oumark.net/TSDFGS/reference/TSDFGS-package.md)
  [`TSDFGS-package`](https://www.oumark.net/TSDFGS/reference/TSDFGS-package.md)
  : TSDFGS: Training Set Determination for Genomic Selection
