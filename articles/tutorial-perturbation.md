# Tutorial: perturbation heritability

[`schint_pert()`](https://LidaWangPSU.github.io/scHINT/reference/schint_pert.md)
asks, for each gene, how much expression variance is attributable to the
**perturbation** a cell received, and how much of that effect depends on
the **cell state** and/or **cell type**. The genetic kernel is replaced
by “same perturbation”.

``` r

head(d, 3)
#>   perturb perturb_group celltype     state        pc1        pc2     GENE_A
#> 1  PERT01        group1      CT2 1.1188078 0.78925056 -0.1370337 0.07107567
#> 2  PERT01        group1      CT3 1.1119961 1.13785263  0.8233651 0.14960499
#> 3  PERT01        group1      CT1 0.4430754 0.08648318  0.4760821 0.38410583
#>        GENE_B    GENE_C     GENE_D
#> 1  0.70179759 1.6306675 -0.3235857
#> 2  0.78308864 0.0930320 -1.0652990
#> 3 -0.05038749 0.8265435 -1.0450563
table(d$perturb_group)
#> 
#> control  group1  group2  group3 
#>     600    1080    1263    1127
```

`NT` marks non-targeting control cells. `GENE_A` has perturbation ×
state variance, `GENE_B` perturbation × cell-type variance, `GENE_C` a
perturbation effect only and `GENE_D` none; group3 perturbations have no
effect.

## Fit several genes at once

``` r

fit <- schint_pert(d, y = c("GENE_A", "GENE_B", "GENE_C", "GENE_D"),
                   perturb = "perturb", context = c("state", "celltype"),
                   covariates = c("pc1", "pc2"), control = "NT", jackknife = TRUE)
keep <- c("P", "PxC:state", "PxC:celltype")
res <- subset(fit$summary, term %in% keep, select = c(gene, term, variance, se))
res$z <- round(res$variance / res$se, 1)
res
#>            gene         term      variance          se    z
#> GENE_A.1 GENE_A            P  0.0272190898 0.008656187  3.1
#> GENE_A.4 GENE_A    PxC:state  0.0507430857 0.014603063  3.5
#> GENE_A.5 GENE_A PxC:celltype  0.0110774848 0.011340415  1.0
#> GENE_B.1 GENE_B            P  0.0435317075 0.015098355  2.9
#> GENE_B.4 GENE_B    PxC:state -0.0027662529 0.002925695 -0.9
#> GENE_B.5 GENE_B PxC:celltype  0.0272631582 0.012551674  2.2
#> GENE_C.1 GENE_C            P  0.0665561848 0.021558539  3.1
#> GENE_C.4 GENE_C    PxC:state  0.0009753755 0.003137058  0.3
#> GENE_C.5 GENE_C PxC:celltype  0.0069372804 0.009553201  0.7
#> GENE_D.1 GENE_D            P -0.0055383708 0.003054088 -1.8
#> GENE_D.4 GENE_D    PxC:state  0.0021470739 0.003810962  0.6
#> GENE_D.5 GENE_D PxC:celltype  0.0092355169 0.006536935  1.4
```

The expensive part (`X'X`) is computed once; each additional gene costs
only a few quadratic forms. `P` is the perturbation variance shared
across cell states, `PxC:state` the variance of the perturbation
response along the state axis and `PxC:celltype` the cell-type-specific
response.

## Perturbation groups

If perturbations fall into groups (pathways, complexes, targets of one
regulator) each group can get its own components:

``` r

fg <- schint_pert(d, c("GENE_A", "GENE_B"), "perturb", context = c("state", "celltype"),
                  perturb_group = "perturb_group", control = "NT")
subset(fg$summary, grepl("^P_group", term) & gene == "GENE_A",
       select = c(term, variance))
#>                     term     variance
#> GENE_A.6        P_group1  0.013266670
#> GENE_A.7      P_group1:C  0.045680191
#> GENE_A.8  P_group1_total  0.058946861
#> GENE_A.9        P_group2  0.014565664
#> GENE_A.10     P_group2:C  0.014481156
#> GENE_A.11 P_group2_total  0.029046820
#> GENE_A.12       P_group3 -0.001479058
#> GENE_A.13     P_group3:C  0.003793489
#> GENE_A.14 P_group3_total  0.002314431
```

`group3` has no simulated effect, `group1` the largest. A named list
works too:
`perturb_group = list(complexA = c("PERT01", "PERT02"), complexB = c(...))`.

## Choosing the context

Use `context = "state"` for a continuous axis (e.g. principal components
of the control cells), `context = "celltype"` for categorical, or both,
as above. `context = NULL` fits perturbation variance only.

## Practical points

- Remove controls (`control`) or give them their own label and exclude
  them; otherwise controls form one very large “perturbation”.
- Perturbations with fewer than `min_cells` (default 2) cells are
  dropped.
- Standard errors use a leave-perturbation-out jackknife (`n_blocks` for
  fewer blocks).
- Expression should be normalized; each gene is centred and scaled.
