# Build a genomic relationship matrix from genotypes

Build a genomic relationship matrix from genotypes

## Usage

``` r
make_grm(geno, standardize = TRUE)
```

## Arguments

- geno:

  Numeric matrix of allele dosages (0/1/2), donors in rows (with row
  names) and SNPs in columns. \`NA\` is mean-imputed.

- standardize:

  If \`TRUE\` (default) each SNP is centred and scaled by
  \`sqrt(2p(1-p))\`; monomorphic SNPs are dropped. If \`FALSE\`,
  \`geno\` is assumed to be standardized already.

## Value

A donor-by-donor matrix \`X X' / m\`.
