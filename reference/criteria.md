# List of training set criteria

Returns the criteria that
[`optTrain`](https://www.oumark.net/TSDFGS/reference/optTrain.md) can
optimize, with the name to give to its `method` argument, the function
that scores a given training set, the direction of optimization, and
whether the criterion uses a test set.

## Usage

``` r
criteria()
```

## Value

A data frame with the columns `method`, `name`, `type`, `score_fun`,
`maximize`, and `uses_test`.

## Details

Criteria marked `uses_test = FALSE` depend on the training set only.
When a test set is given to `optTrain`, such a criterion ignores it
(with a warning). One page per criterion on the package website gives
the definition, the original reference, and the license review.

## Author

Jen-Hsiang Ou

## Examples

``` r
criteria()
#>          method                              name              type
#> 1        rScore                           r-score       model-based
#> 2           PEV              PEVmean (ridge, PCs)       model-based
#> 3            CD                   CD (ridge, PCs)       model-based
#> 4          AOpt                      A-optimality    optimal design
#> 5          DOpt                      D-optimality    optimal design
#> 6          EOpt                      E-optimality    optimal design
#> 7        PEVmax                            PEVmax       model-based
#> 8          GOpt            G-optimality (GOPTPEV)    optimal design
#> 9          MSPE          MSPE of ridge regression       model-based
#> 10    GVoverall                         GVoverall genomic variation
#> 11    GVaverage                         GVaverage genomic variation
#> 12       AvgGRM                           Avg_GRM      relationship
#> 13   AvgGRMself                      Avg_GRM_self      relationship
#> 14 AvgGRMminmax                    Avg_GRM_MinMax      relationship
#> 15       CDmean                            CDmean       model-based
#> 16     CDmeanV2                        CDmean(v2)       model-based
#> 17          PAM Partitioning around medoids (PAM)        clustering
#> 18        CDmax                CDmax (ridge, PCs)       model-based
#>               score_fun maximize uses_test
#> 1               r_score     TRUE      TRUE
#> 2             pev_score    FALSE      TRUE
#> 3              cd_score     TRUE      TRUE
#> 4            aopt_score    FALSE     FALSE
#> 5            dopt_score     TRUE     FALSE
#> 6            eopt_score    FALSE     FALSE
#> 7          pevmax_score    FALSE      TRUE
#> 8            gopt_score    FALSE      TRUE
#> 9            mspe_score    FALSE      TRUE
#> 10     gv_overall_score     TRUE     FALSE
#> 11     gv_average_score     TRUE     FALSE
#> 12        avg_grm_score     TRUE      TRUE
#> 13   avg_grm_self_score    FALSE     FALSE
#> 14 avg_grm_minmax_score     TRUE      TRUE
#> 15         cdmean_score     TRUE      TRUE
#> 16        cdmean2_score     TRUE      TRUE
#> 17            pam_score    FALSE      TRUE
#> 18          cdmax_score     TRUE      TRUE
```
