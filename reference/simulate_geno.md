# Simulate genotypes with local LD

Dosages (0/1/2) from a Gaussian copula with AR(1) correlation between
neighbouring SNPs and uniform minor-allele frequencies.

## Usage

``` r
simulate_geno(n = 400, m = 250, rho = 0.6, maf = c(0.05, 0.5), seed = 1)
```

## Arguments

- n:

  Number of donors.

- m:

  Number of SNPs.

- rho:

  Correlation of adjacent SNPs.

- maf:

  Range of minor-allele frequencies.

- seed:

  Random seed.

## Value

Integer matrix with donor IDs \`D001...\` and SNP IDs \`SNP001...\`.
