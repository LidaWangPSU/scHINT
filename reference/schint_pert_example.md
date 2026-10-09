# Example perturb-seq data for scHINT

Simulated cells (\`simulate_pert()\`): \`perturb\` (60 perturbations in
three groups plus \`"NT"\` controls), \`perturb_group\`, cell type
\`celltype\`, continuous \`state\`, covariates \`pc1\`, \`pc2\` and
expression of four genes \`GENE_A\`-\`GENE_D\`. GENE_A has perturbation
x state variance, GENE_B perturbation x cell type, GENE_C perturbation
only and GENE_D none; group3 perturbations have no effect.

## Usage

``` r
schint_pert_example
```

## Format

A data frame with one row per cell.

## Source

\`data-raw/make_pert_example.R\`
