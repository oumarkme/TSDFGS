# Articles

### All vignettes

- [License review of the
  criteria](https://www.oumark.net/TSDFGS/articles/criteria-license-review.md):

  Which criteria TSDFGS includes, which it leaves out, and why.

- [Choosing a
  criterion](https://www.oumark.net/TSDFGS/articles/criteria.md):

  The training set criteria in TSDFGS, one page each, and what recent
  comparisons found.

- [A-optimality](https://www.oumark.net/TSDFGS/articles/criterion-aopt.md):

  A-optimality minimizes the average variance of the ridge coefficients.

- [Avg_GRM](https://www.oumark.net/TSDFGS/articles/criterion-avggrm.md):

  The mean genomic relationship between the training set and the target.

- [Avg_GRM_MinMax](https://www.oumark.net/TSDFGS/articles/criterion-avggrmminmax.md):

  Combines the two: high relationship to the target and low relationship
  within the training set.

- [Avg_GRM_self](https://www.oumark.net/TSDFGS/articles/criterion-avggrmself.md):

  The mean relationship within the training set; minimizing it spreads
  the training individuals apart.

- [CD (ridge,
  PCs)](https://www.oumark.net/TSDFGS/articles/criterion-cd.md):

  CD is the squared correlation between true and predicted genetic value
  (reliability).

- [CDmax (ridge,
  PCs)](https://www.oumark.net/TSDFGS/articles/criterion-cdmax.md):

  CDmax is the largest coefficient of determination among the target
  individuals.

- [CDmean](https://www.oumark.net/TSDFGS/articles/criterion-cdmean.md):

  The mean coefficient of determination of the GBLUP contrast between
  each target individual and the target mean (Rincent et al.

- [CDmean(v2)](https://www.oumark.net/TSDFGS/articles/criterion-cdmean2.md):

  The heuristic form of CDmean from Chen et al.

- [D-optimality](https://www.oumark.net/TSDFGS/articles/criterion-dopt.md):

  D-optimality maximizes the information content of the training set,
  the determinant of the ridge information matrix.

- [E-optimality](https://www.oumark.net/TSDFGS/articles/criterion-eopt.md):

  E-optimality limits the worst-estimated direction: the inverse of the
  smallest eigenvalue of the information matrix.

- [G-optimality
  (GOPTPEV)](https://www.oumark.net/TSDFGS/articles/criterion-gopt.md):

  G-optimality minimizes the largest eigenvalue of the PEV matrix of the
  target, the variance in the worst direction.

- [GVaverage](https://www.oumark.net/TSDFGS/articles/criterion-gvaverage.md):

  GVaverage is the mean genetic variance of the training individuals,
  the trace of their relationship matrix divided by $`n`$ (TSDFGS
  returns the trace).

- [GVoverall](https://www.oumark.net/TSDFGS/articles/criterion-gvoverall.md):

  GVoverall measures the total genetic variation covered by the training
  set: the generalized variance of its relationship matrix.

- [MSPE of ridge
  regression](https://www.oumark.net/TSDFGS/articles/criterion-mspe.md):

  The expected mean squared prediction error of ridge regression for the
  target individuals (Chen et al.

- [Partitioning around medoids
  (PAM)](https://www.oumark.net/TSDFGS/articles/criterion-pam.md):

  The k-medoids cost: the training individuals act as medoids of the
  target.

- [PEVmean (ridge,
  PCs)](https://www.oumark.net/TSDFGS/articles/criterion-pev.md):

  The prediction error variance (PEV) measures how uncertain the model
  predictions are for the test set.

- [PEVmax](https://www.oumark.net/TSDFGS/articles/criterion-pevmax.md):

  PEVmax is the largest prediction error variance among the test
  individuals, so it protects the worst-predicted one.

- [r-score](https://www.oumark.net/TSDFGS/articles/criterion-rscore.md):

  The r-score is the expected Pearson correlation between predicted and
  observed values of the test set, written in closed form for ridge
  regression.

- [Getting
  started](https://www.oumark.net/TSDFGS/articles/getting-started.md):

  From a genotype matrix to a training set and a training set size.

- [Installation and
  environments](https://www.oumark.net/TSDFGS/articles/installation.md):

  Requirements, tested versions, and ways to install TSDFGS: R, Docker,
  conda, and Apptainer.

- [How the scores are
  computed](https://www.oumark.net/TSDFGS/articles/methods.md):

  The identities that make the v3 scores fast without changing their
  values.

- [Moving from v2 to
  v3](https://www.oumark.net/TSDFGS/articles/migrating-from-v2.md):

  What changed in v3, what did not, and how much faster it is.
