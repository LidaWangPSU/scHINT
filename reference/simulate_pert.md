# Simulate perturb-seq-like data for the perturbation heritability model

Cells carry one of \`n_perturb\` perturbations (in three groups) or a
non-targeting control (\`"NT"\`), have a continuous cell \`state\` and
one of three cell types. Gene expression has perturbation (\`p\`),
perturbation x state (\`pxs\`) and perturbation x cell-type (\`pxc\`)
variance, scaled per perturbation group by \`group_scale\`.

## Usage

``` r
simulate_pert(
  genes = list(GENE_A = list(p = 0.1, pxs = 0.1, pxc = 0, cov = 0.05), GENE_B = list(p =
    0.1, pxs = 0, pxc = 0.1, cov = 0.05), GENE_C = list(p = 0.15, pxs = 0, pxc = 0, cov =
    0.05), GENE_D = list(p = 0, pxs = 0, pxc = 0, cov = 0.05)),
  n_perturb = 60,
  cells = c(30, 80),
  n_control = 600,
  group_scale = c(1, 0.5, 0),
  seed = 1
)
```

## Arguments

- genes:

  Named list of variances \`p\`, \`pxs\`, \`pxc\`, \`cov\` (rest is
  noise).

- n_perturb:

  Number of perturbations.

- cells:

  Range of cells per perturbation.

- n_control:

  Number of control cells.

- group_scale:

  Multiplier on the perturbation variances in each of the three groups.

- seed:

  Random seed.

## Value

A data frame with one row per cell.
