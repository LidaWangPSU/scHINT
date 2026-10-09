# Read a PLINK binary fileset as a dosage matrix

Requires the suggested package \`BEDMatrix\`.

## Usage

``` r
read_plink(
  root,
  snps = NULL,
  impute = c("none", "avg"),
  id = c("IID", "FID:IID")
)
```

## Arguments

- root:

  Path prefix of the \`.bed/.bim/.fam\` files.

- snps:

  Optional character vector of SNP IDs to keep (read lazily).

- impute:

  \`"none"\` or \`"avg"\` (mean imputation).

- id:

  Which ID to use for row names: \`"IID"\` (default) or \`"FID:IID"\`.

## Value

A list with \`geno\` (donors x SNPs), \`bim\` and \`fam\`.
