# How the scores are computed

The three criteria compare a training matrix **X** (n × p) with a test
matrix **X0** (n0 × p). The v2 code formed matrices with n0 rows and n0
columns at every step of the search. The v3 code uses matrices whose
size is the smaller of n and p. The values are the same.

Let λ = 1/p for the r-score and PEV, and λ = 1 for CD, and write **G** =
**X**’**X** (p × p).

## The push-through identity

The ridge solution can be written in a primal or a dual form:

``` math
\mathbf{X}'(\mathbf{X}\mathbf{X}' + \lambda \mathbf{I}_n)^{-1} = (\mathbf{X}'\mathbf{X} + \lambda \mathbf{I}_p)^{-1}\mathbf{X}'
```

The left side needs an n × n inverse, the right side a p × p inverse. v3
picks the smaller one:

- **p ≤ n** (primal): everything is a p × p matrix, and the test set
  enters only through **X0**’**X0**.
- **p \> n** (dual): everything is an n × n matrix.

## r-score

With **M** = (**G** + λ**I**)^-1 **G**, **S0** = c² **X0**’**X0**, and c
= 1 − 1/n0 (the constant of v2.4.2, kept as published):

``` math
q_1 = n_0 - 1 + \operatorname{tr}(\mathbf{S}_0), \quad
q_2 = \operatorname{tr}(\mathbf{M}\mathbf{B}\mathbf{S}_0) + \operatorname{tr}(\mathbf{M}'\mathbf{S}_0\mathbf{M}), \quad
q_{12} = c\,\operatorname{tr}(\mathbf{X}_0'\mathbf{X}_0\mathbf{M})
```

and the r-score is q₁₂ / √(q₁ q₂), with **B** = (**G** + λ**I**)^-1. The
test set enters only through the p × p matrix **X0**’**X0**, which the
search computes once when the test set is fixed.

## PEV and CD

- **PEV** is tr(**B** **X0e**’**X0e**) with an intercept column added
  (**X0e**), a p × p trace in the primal form.
- **CD** is the mean over test rows of **x**\_i’ **B** **G** **B**
  **x**\_i / **x**\_i’**x**\_i, a quadratic form per row, with no n0 ×
  n0 matrix.

## What does not change

- The criterion values equal those of v2.4.2: to about 1e-7, the
  accuracy of the single-precision numbers v2 returned. The tests
  compare against an independent R transcription of the v2 formulas to
  1e-9, and against values recorded from v2.4.2.
- The search moves, the order of the random numbers, the default number
  of iterations, and the stopping rule are those of v2.4.2. With
  `options(TSDFGS.legacy_float = TRUE)` and the same seed,
  [`optTrain()`](https://www.oumark.net/TSDFGS/reference/optTrain.md)
  returns the training sets recorded from v2.4.2 for all twelve routes
  (with or without sub-populations, targeted or not, three criteria).

## Where the time goes

Per move, v2 spent O(n0² p + n³) operations; v3 spends O(n0 p² + p³) for
an untargeted search and O(p³ + n p²) for a targeted one, plus the cost
of copying the training rows. The slow case for both is large p together
with large n; reducing p with the leading PCs is the practical remedy,
and it changes the criterion, so choose it before comparing results.

The C++ code uses Eigen through the RcppEigen package. R is the
interface: it checks inputs, draws the starting set, and returns the
result.
