# Tutorial: perturbation heritability

[`schint_pert()`](https://LidaWangPSU.github.io/scHINT/reference/schint_pert.md)
asks, for each gene, how much expression variance is attributable to the
**perturbation** a cell received (`P`), and how much of that effect
depends on the cellular context (`PxC`). The genetic kernel is replaced
by “same perturbation”. The context is **one variable at a time**: a
continuous **cell state** or a categorical **cell type**.

``` r

head(d[, c("perturb", "celltype", "state", "pc1", "pc2", "GENE_A")], 3)
#>   perturb celltype     state        pc1        pc2     GENE_A
#> 1  PERT01      CT2 1.1188078 0.78925056 -0.1370337 0.07107567
#> 2  PERT01      CT3 1.1119961 1.13785263  0.8233651 0.14960499
#> 3  PERT01      CT1 0.4430754 0.08648318  0.4760821 0.38410583
table(d$celltype)
#> 
#>  CT1  CT2  CT3 
#> 2030 1259  781
```

`NT` marks non-targeting control cells and is removed with
`control = "NT"`. `GENE_A` has perturbation × cell-state variance,
`GENE_B` perturbation × cell-type variance, `GENE_C` a perturbation
effect only and `GENE_D` none. Several genes can be fitted in one call;
the expensive part (`X'X`) is computed only once.

``` r

genes <- c("GENE_A", "GENE_B", "GENE_C", "GENE_D")
```

## Example 1: perturbation × cell type (P×CT)

Set `context` to the cell-type column. `PxC` is the variance of the
perturbation response that is specific to a cell type (cells with the
same perturbation **and** the same cell type).

``` r

fit_ct <- schint_pert(d, y = genes, perturb = "perturb", context = "celltype",
                      covariates = c("pc1", "pc2"), control = "NT", jackknife = TRUE)
res_ct <- subset(fit_ct$summary, term %in% c("P", "PxC"), select = c(gene, term, variance, se))
res_ct$z <- round(res_ct$variance / res_ct$se, 1)
res_ct
#>            gene term     variance          se    z
#> GENE_A.1 GENE_A    P  0.023510186 0.008558527  2.7
#> GENE_A.2 GENE_A  PxC  0.021032288 0.012095408  1.7
#> GENE_B.1 GENE_B    P  0.043718630 0.015134585  2.9
#> GENE_B.2 GENE_B  PxC  0.026748300 0.012433431  2.2
#> GENE_C.1 GENE_C    P  0.066462888 0.021654863  3.1
#> GENE_C.2 GENE_C  PxC  0.007168734 0.009545346  0.8
#> GENE_D.1 GENE_D    P -0.005701530 0.002990322 -1.9
#> GENE_D.2 GENE_D  PxC  0.009668077 0.006488487  1.5
```

`GENE_B` has the largest P×CT component.

## Example 2: perturbation × cell state (P×CS)

Set `context` to a continuous cell-state column (for example a principal
component of the control cells). `PxC` is the variance of the
perturbation response that changes along that axis.

``` r

fit_cs <- schint_pert(d, y = genes, perturb = "perturb", context = "state",
                      covariates = c("pc1", "pc2"), control = "NT", jackknife = TRUE)
res_cs <- subset(fit_cs$summary, term %in% c("P", "PxC"), select = c(gene, term, variance, se))
res_cs$z <- round(res_cs$variance / res_cs$se, 1)
res_cs
#>            gene term     variance          se    z
#> GENE_A.1 GENE_A    P  0.031529450 0.008370534  3.8
#> GENE_A.2 GENE_A  PxC  0.051115680 0.014668888  3.5
#> GENE_B.1 GENE_B    P  0.053930007 0.015738307  3.4
#> GENE_B.2 GENE_B  PxC -0.001594474 0.002814178 -0.6
#> GENE_C.1 GENE_C    P  0.069294137 0.021355071  3.2
#> GENE_C.2 GENE_C  PxC  0.001161913 0.003183737  0.4
#> GENE_D.1 GENE_D    P -0.001956778 0.002082991 -0.9
#> GENE_D.2 GENE_D  PxC  0.002472314 0.003809523  0.6
```

`GENE_A` has a clear P×CS component, while `GENE_C` (perturbation effect
only) and `GENE_D` (none) do not. In this simulation the cell state
partly tracks cell type, so `GENE_A` also shows a weak cell-type term in
Example 1; with real data fit both and compare.

## Practical points

- Remove controls (`control`) or give them their own label and exclude
  them; otherwise controls form one very large “perturbation”.
- Perturbations with fewer than `min_cells` (default 2) cells are
  dropped.
- Standard errors use a leave-perturbation-out jackknife (`n_blocks` for
  fewer blocks).
- Expression should be normalized; each gene is standardized.
