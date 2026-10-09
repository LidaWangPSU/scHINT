# Split SNPs into groups (for multiple GRMs)

Split SNPs into groups (for multiple GRMs)

## Usage

``` r
split_snps(geno, n, by = c("maf", "position", "random"), seed = 1)
```

## Arguments

- geno:

  Genotype matrix (donors x SNPs).

- n:

  Number of groups.

- by:

  \`"maf"\` (equal-count minor-allele-frequency bins, rarest first),
  \`"position"\` (contiguous blocks in column order, i.e. genomic order
  if columns are sorted) or \`"random"\`.

- seed:

  Seed used when \`by = "random"\`.

## Value

A list of \`n\` genotype matrices; use \`lapply(split_snps(geno, 3),
make_grm)\` to obtain a list of GRMs for the \`grm\` argument of the
model functions.
