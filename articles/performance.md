# Performance

## Computational complexity

For \\M\\ cells from \\N\\ donors with \\K\\ context/covariate
variables, directly forming the cell-by-cell kernels costs \\O(M^2K^2)\\
time and \\O(M^2)\\ memory. scHINT aggregates cells into donor-level
sufficient statistics (see
[Model](https://LidaWangPSU.github.io/scHINT/articles/model.md), section
5), so

- time: one linear pass over the cells plus \\O(N^2K^2)\\ for the
  donor-level quadratic forms;
- memory: \\O(N^2)\\ for the GRM (one matrix per GRM) – independent of
  \\M\\;
- jackknife: delete-one-donor or block standard errors come from the
  same donor-level quantities, with no re-fitting.

The number of cells therefore affects only the linear pass; the number
of donors dominates for large cohorts.

## Scaling with cells and donors

``` r

geno <- simulate_geno(600, 200, seed = 1)
spec <- list(G = list(g = 0.2, gxc = 0.1, i = 0.1, cov = 0.05))
time_fit <- function(n_donors, cells_per_donor) {
  s <- simulate_schint(geno[seq_len(n_donors), ], genes = spec,
                       cells = rep(cells_per_donor, 2), seed = 1)
  t0 <- proc.time()[["elapsed"]]
  schint_cell(s$cell, "G", "donor", geno = geno, context = "state",
              covariates = c("pc1", "pc2", "age", "sex", "batch"), jackknife = TRUE)
  proc.time()[["elapsed"]] - t0
}
cells <- data.frame(cells_per_donor = c(10, 20, 40, 80, 160))
cells$n_cells <- 300 * cells$cells_per_donor
cells$seconds <- sapply(cells$cells_per_donor, function(m) time_fit(300, m))
cells
#>   cells_per_donor n_cells seconds
#> 1              10    3000   0.060
#> 2              20    6000   0.060
#> 3              40   12000   0.071
#> 4              80   24000   0.099
#> 5             160   48000   0.155
donors <- data.frame(n_donors = c(100, 200, 400, 600))
donors$seconds <- sapply(donors$n_donors, function(n) time_fit(n, 30))
donors
#>   n_donors seconds
#> 1      100   0.029
#> 2      200   0.049
#> 3      400   0.091
#> 4      600   0.134
```

![](performance_files/figure-html/bench-plot-1.png)

Timings include the delete-one-donor jackknife and are for a single gene
on one core. (The manuscript benchmarks scHINT against scGeneHE on up to
one million cells.)

## Memory

Dense \\N\times N\\ GRMs dominate: 8 bytes × \\N^2\\ per GRM (e.g. 1,000
donors ≈ 8 MB; 10,000 donors ≈ 800 MB). Several GRMs and the
all-ones/identity kernels used internally add a small multiple; for
perturbation screens the “same perturbation” kernels are stored as
vectors, so thousands of perturbations are fine.

## Optimization checklist

1.  Fit genes in parallel
    ([`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html),
    a job array); each fit is independent.
2.  Restrict genotypes to the cis window of each gene (fewer SNPs →
    faster GRM).
3.  Use pre-computed GRMs (`grm =`) if the same donors/SNPs are reused.
4.  Use `n_blocks` (e.g. 50-100) instead of one block per donor for very
    large cohorts.
5.  For perturbation screens, pass many genes to one
    [`schint_pert()`](https://LidaWangPSU.github.io/scHINT/reference/schint_pert.md)
    call to share `X'X`.
6.  Keep `context`/`covariates` to what you need: cost grows with the
    square of their number.
