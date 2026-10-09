# Tutorial: a genome-wide workflow

This page shows how to run scHINT over many genes on real data. Code
chunks are not evaluated here because they need your own files.

## 1. Inputs

- PLINK genotypes (`.bed/.bim/.fam`), donor IDs in the `IID` column.
- A cell metadata table with the donor ID, cell-state PCs and
  covariates.
- A normalized expression matrix (cells × genes) with the same cell
  order as the metadata.
- A gene table with `gene`, `chr`, `start`, `end` (same genome build as
  the genotypes).

``` r

library(scHINT)

meta  <- readRDS("meta.rds")            # data.frame, one row per cell
expr  <- readRDS("expr.rds")            # matrix cells x genes (normalized)
genes <- read.table("genes.tsv", header = TRUE)
```

## 2. State axes

Principal components of the normalized expression **within the cell
type** are used as cell states. Compute them once per cell type and add
to `meta`:

``` r

pcs <- prcomp(expr[, hvg], center = TRUE, scale. = TRUE, rank. = 10)$x
meta[, paste0("PC", 1:10)] <- pcs
```

## 3. One gene

``` r

fit_gene <- function(g, window = 5e5) {
  # cis window of the gene: +/- 500 kb around the gene body
  bim  <- plink$bim
  snps <- select_cis_snps(bim, g$chr, g$start, g$end, window = window)
  if (length(snps) < 10) return(NULL)
  geno <- read_plink(plink_root, snps = snps, impute = "avg")$geno

  d <- cbind(meta, y = expr[, g$gene])
  schint_cell(d, y = "y", id = "donor", geno = geno,
              context    = paste0("PC", 1:5),
              covariates = c(paste0("PC", 6:10), "pct_mt", "age", "sex"),
              jackknife  = TRUE)
}
```

`read_plink(..., snps = )` only reads the requested variants, so the
whole genotype file is never loaded.

## 4. Many genes

``` r

plink <- read_plink(plink_root)         # once, only to get the bim table (or read bim directly)

library(parallel)
res <- mclapply(seq_len(nrow(genes)), function(i) {
  f <- tryCatch(fit_gene(genes[i, ]), error = function(e) NULL)
  if (is.null(f)) return(NULL)
  cbind(gene = genes$gene[i], f$summary)
}, mc.cores = 8)
res <- do.call(rbind, res)
write.table(res, "schint_summary.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
```

On a cluster, split `genes` into chunks and submit each as one job (the
fits are independent).

## 5. Using GCTA GRMs

If you already have per-gene GRMs from GCTA (`--make-grm`), skip the
genotypes:

``` r

K <- read_grm_bin("gene_prefix")        # donors x donors, IDs as dimnames
schint_cell(d, "y", "donor", grm = K, context = paste0("PC", 1:5))
```

Pass a named list of GRMs for a multi-GRM model.

## 6. Summarizing

Typical summaries use the heritabilities (`h2`) from `summary`: the
distribution of `h2` for `G` and `GxC`, the fraction of genes with
\\z=\hat\sigma^2\_{G\times C}/\mathrm{SE}\>\\ a threshold, or the mean
of `GxC` over genes in a gene set. Weight by inverse variance (`1/se^2`)
when averaging across genes of unequal precision.

## 7. Choices that matter

- **Cis window.** Heritability grows with the window; keep it constant
  across genes (±500 kb is a common choice).
- **Number of state PCs** in `context` – more axes give more flexibility
  but more parameters; 5 is a reasonable default.
- **Donor count.** Standard errors scale with the number of donors, not
  cells; hundreds of donors are needed for precise per-gene estimates,
  but averages over many genes are informative with fewer.
