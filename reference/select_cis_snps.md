# Select cis SNPs around a gene

Select cis SNPs around a gene

## Usage

``` r
select_cis_snps(bim, chr, start, end, window = 5e+05)
```

## Arguments

- bim:

  Data frame with columns \`chr\`, \`snp\`, \`pos\` (e.g. the \`bim\`
  returned by \[read_plink()\]).

- chr, start, end:

  Gene chromosome and coordinates.

- window:

  Window (bp) added on each side of the gene body (default 500 kb).

## Value

Character vector of SNP IDs.
