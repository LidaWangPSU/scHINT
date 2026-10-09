# Users manual

## Installation

``` r

install.packages("remotes")
remotes::install_github("LidaWangPSU/scHINT")
```

scHINT depends only on base R (`stats`, `utils`). `BEDMatrix` is
optional and only needed to read PLINK files with
[`read_plink()`](https://LidaWangPSU.github.io/scHINT/reference/read_plink.md).

## Data preparation

1.  **Choose the cell population.** scHINT is applied within one cell
    type for cell-state analyses
    ([`schint_cell()`](https://LidaWangPSU.github.io/scHINT/reference/schint_cell.md)),
    or across cell types for pseudobulk data
    ([`schint_sample()`](https://LidaWangPSU.github.io/scHINT/reference/schint_sample.md)).
2.  **Normalize expression** (log-normalization, regress out technical
    variation if needed).
3.  **Define the cell context**: the leading principal components
    (e.g. 5) of the normalized expression *within* the cell type are a
    convenient continuous cell-state axis; any numeric or categorical
    per-cell variable can be used.
4.  **Collect covariates**: cell-level (e.g. number of UMIs, percent
    mitochondrial, extra PCs) and donor-level (age, sex, genotype PCs,
    batch).
5.  **Extract cis genotypes** of each gene: SNPs within ±500 kb of the
    gene body, dosage matrix with donor IDs as row names.
6.  Match donor IDs between genotypes and metadata.

## Usage

### One gene

``` r

library(scHINT)

fit <- schint_cell(
  data       = cells,                  # one row per cell
  y          = "GENE",                 # expression column
  id         = "donor",
  geno       = geno_cis,               # donors x cis-SNPs
  context    = paste0("PC", 1:5),      # G×PC1..PC5 interactions
  covariates = c("PC6", "PC7", "pct_mt", "age", "sex"),
  jackknife  = TRUE
)
fit$summary
```

Key choices:

| decision | argument | default |
|----|----|----|
| genetic input / number of GRMs | `geno` (one GRM from all cis SNPs) or `grm` (a GRM or a list of GRMs) | – |
| which context interacts with genetics | `context` | none (no G×C term) |
| other covariates (main effects) | `covariates` | none |
| donor effect | `ind_effect` | `TRUE` |
| standard errors | `jackknife`, `n_blocks` | off; one block per donor |

### Many genes

Loop over genes, subsetting the genotypes to each gene’s cis window (see
the [real-data
workflow](https://LidaWangPSU.github.io/scHINT/articles/tutorial-workflow.md)).
Fits are independent, so they can be parallelized across genes
(e.g. [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)).

### Pseudobulk across cell types

``` r

fit <- schint_sample(pseudobulk, y = "GENE", id = "donor", celltype = "celltype",
                     geno = geno_cis, covariates = c("age", "sex"), jackknife = TRUE)
```

### Perturbation screens

``` r

fit <- schint_pert(cells, y = c("GENE1", "GENE2"), perturb = "target",
                   context = "state", control = "NT")   # or "celltype"
```

See the tutorials for complete, runnable examples.
