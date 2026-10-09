# Example data for scHINT

Simulated genotypes and expression. Expression of three genes with known
variance architecture is simulated (see \`truth\`).

## Usage

``` r
schint_example
```

## Format

A list with

- geno:

  400 x 250 integer dosage matrix, donor IDs \`D001\`... as row names.

- snp:

  SNP IDs and (nominal) positions.

- cell:

  12,000 cells of one cell type: \`donor\`, continuous \`state\`,
  categorical \`subtype\`, cell covariates \`pc1\`, \`pc2\`, donor
  covariates \`age\`, \`sex\`, \`batch\`, and expression \`GENE_A\`,
  \`GENE_B\`, \`GENE_C\`.

- sample:

  2,000 donor x cell-type (5 types) pseudobulk values with the same
  genes.

- truth:

  Variances used in the simulation (\`g\`, \`gxc\`, \`i\`, \`cov\`,
  \`noise\`).

## Source

\`data-raw/make_example.R\`
