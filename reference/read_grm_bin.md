# Read / write a GCTA binary GRM

Read / write a GCTA binary GRM

## Usage

``` r
read_grm_bin(prefix)

write_grm_bin(K, prefix, n_snps = 1L)
```

## Arguments

- prefix:

  Path prefix of the \`.grm.bin\`/\`.grm.id\` files.

- K:

  Symmetric GRM with donor IDs as dimnames.

- n_snps:

  Number of SNPs to record in the \`.grm.N.bin\` file.

## Value

\`read_grm_bin()\` returns a full symmetric matrix with donor IDs
(second column of \`.grm.id\`) as dimnames.
