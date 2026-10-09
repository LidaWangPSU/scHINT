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

## One interaction at a time

The perturbation can interact with **either** a continuous cell state
**or** a categorical cell type in one model (`context` takes a single
column); to study both, fit the model twice. Several genes can be fitted
at once, and the expensive part (`X'X`) is computed only once; each
additional gene costs a few quadratic forms.

``` r

genes <- c("GENE_A", "GENE_B", "GENE_C", "GENE_D")
fit_state <- schint_pert(d, y = genes, perturb = "perturb", context = "state",
                         covariates = c("pc1", "pc2"), control = "NT", jackknife = TRUE)
fit_type  <- schint_pert(d, y = genes, perturb = "perturb", context = "celltype",
                         covariates = c("pc1", "pc2"), control = "NT", jackknife = TRUE)
pick <- function(f, label) {
  s <- subset(f$summary, term %in% c("P", "PxC"), select = c(gene, term, variance, se))
  s$interaction <- label
  s$z <- round(s$variance / s$se, 1)
  s
}
rbind(pick(fit_state, "P x cell state"), pick(fit_type, "P x cell type"))
#>             gene term     variance          se    interaction    z
#> GENE_A.1  GENE_A    P  0.031529450 0.008370534 P x cell state  3.8
#> GENE_A.2  GENE_A  PxC  0.051115680 0.014668888 P x cell state  3.5
#> GENE_B.1  GENE_B    P  0.053930007 0.015738307 P x cell state  3.4
#> GENE_B.2  GENE_B  PxC -0.001594474 0.002814178 P x cell state -0.6
#> GENE_C.1  GENE_C    P  0.069294137 0.021355071 P x cell state  3.2
#> GENE_C.2  GENE_C  PxC  0.001161913 0.003183737 P x cell state  0.4
#> GENE_D.1  GENE_D    P -0.001956778 0.002082991 P x cell state -0.9
#> GENE_D.2  GENE_D  PxC  0.002472314 0.003809523 P x cell state  0.6
#> GENE_A.11 GENE_A    P  0.023510186 0.008558527  P x cell type  2.7
#> GENE_A.21 GENE_A  PxC  0.021032288 0.012095408  P x cell type  1.7
#> GENE_B.11 GENE_B    P  0.043718630 0.015134585  P x cell type  2.9
#> GENE_B.21 GENE_B  PxC  0.026748300 0.012433431  P x cell type  2.2
#> GENE_C.11 GENE_C    P  0.066462888 0.021654863  P x cell type  3.1
#> GENE_C.21 GENE_C  PxC  0.007168734 0.009545346  P x cell type  0.8
#> GENE_D.11 GENE_D    P -0.005701530 0.002990322  P x cell type -1.9
#> GENE_D.21 GENE_D  PxC  0.009668077 0.006488487  P x cell type  1.5
```

`P` is the perturbation variance shared across contexts, and `PxC` the
variance of the perturbation response along the chosen axis. `GENE_A`
has a clear **P × cell state** component, `GENE_B` a clear **P × cell
type** component, and `GENE_C`/`GENE_D` neither.

## Perturbation groups

If perturbations fall into groups (pathways, complexes, targets of one
regulator) each group can get its own components:

``` r

fg <- schint_pert(d, c("GENE_A", "GENE_B"), "perturb", context = "state",
                  perturb_group = "perturb_group", control = "NT")
subset(fg$summary, grepl("^P_group", term) & gene == "GENE_A",
       select = c(term, variance))
#>                     term      variance
#> GENE_A.4        P_group1  1.789680e-02
#> GENE_A.5      P_group1:C  3.387848e-02
#> GENE_A.6  P_group1_total  5.177528e-02
#> GENE_A.7        P_group2  1.355287e-02
#> GENE_A.8      P_group2:C  1.720480e-02
#> GENE_A.9  P_group2_total  3.075768e-02
#> GENE_A.10       P_group3 -2.267134e-05
#> GENE_A.11     P_group3:C  1.720171e-04
#> GENE_A.12 P_group3_total  1.493458e-04
```

`group3` has no simulated effect, `group1` the largest. A named list
works too:
`perturb_group = list(complexA = c("PERT01", "PERT02"), complexB = c(...))`.

## Choosing the context

Use `context = "state"` for a continuous axis (e.g. a principal
component of the control cells) and `context = "celltype"` for a
categorical one. `context = NULL` fits perturbation variance only.

## Practical points

- Remove controls (`control`) or give them their own label and exclude
  them; otherwise controls form one very large “perturbation”.
- Perturbations with fewer than `min_cells` (default 2) cells are
  dropped.
- Standard errors use a leave-perturbation-out jackknife (`n_blocks` for
  fewer blocks).
- Expression should be normalized; each gene is centred and scaled.
