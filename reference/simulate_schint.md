# Simulate cell-level and sample-level expression for the example data

Cell level: one cell type, \`cells\` per donor, a continuous cell state,
donor and cell covariates. Sample level: \`n_celltypes\` cell types per
donor. Each gene has its own variance architecture (\`genes\`).

## Usage

``` r
simulate_schint(
  geno,
  genes = list(GENE_A = list(g = 0.2, gxc = 0.15, i = 0.1, cov = 0.05), GENE_B = list(g =
    0.25, gxc = 0, i = 0.1, cov = 0.05), GENE_C = list(g = 0, gxc = 0, i = 0.15, cov =
    0.05)),
  n_causal = 30,
  cells = c(15, 45),
  n_celltypes = 5,
  seed = 1
)
```

## Arguments

- geno:

  Donors x SNPs dosage matrix (row names = donor IDs).

- genes:

  Named list; each element is a list with \`g\` (main genetic variance),
  \`gxc\` (genotype x state / cell-type-specific genetic variance),
  \`i\` (donor variance), \`cov\` (covariate variance); the remainder is
  noise.

- n_causal:

  Number of causal SNPs per gene.

- cells:

  Range of cells per donor.

- n_celltypes:

  Number of cell types for the sample-level data.

- seed:

  Random seed.

## Value

A list with \`cell\`, \`sample\` and \`truth\` (the generating
variances).
