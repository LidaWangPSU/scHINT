#' Read a PLINK binary fileset as a dosage matrix
#'
#' Requires the suggested package `BEDMatrix`.
#'
#' @param root Path prefix of the `.bed/.bim/.fam` files.
#' @param snps Optional character vector of SNP IDs to keep (read lazily).
#' @param impute `"none"` or `"avg"` (mean imputation).
#' @param id Which ID to use for row names: `"IID"` (default) or `"FID:IID"`.
#' @return A list with `geno` (donors x SNPs), `bim` and `fam`.
#' @export
read_plink <- function(root, snps = NULL, impute = c("none", "avg"), id = c("IID", "FID:IID")) {
  impute <- match.arg(impute)
  id <- match.arg(id)
  if (!requireNamespace("BEDMatrix", quietly = TRUE)) {
    stop("Package 'BEDMatrix' is required: install.packages('BEDMatrix')", call. = FALSE)
  }
  root <- path.expand(root)
  bim <- utils::read.table(paste0(root, ".bim"), header = FALSE, stringsAsFactors = FALSE,
                           col.names = c("chr", "snp", "cm", "pos", "a1", "a2"))
  fam <- utils::read.table(paste0(root, ".fam"), header = FALSE, stringsAsFactors = FALSE,
                           col.names = c("fid", "iid", "pat", "mat", "sex", "pheno"))
  bm <- BEDMatrix::BEDMatrix(paste0(root, ".bed"), n = nrow(fam), p = nrow(bim))
  cols <- if (is.null(snps)) seq_len(nrow(bim)) else which(bim$snp %in% snps)
  geno <- as.matrix(bm[, cols, drop = FALSE])
  colnames(geno) <- bim$snp[cols]
  rownames(geno) <- if (id == "IID") fam$iid else paste(fam$fid, fam$iid, sep = ":")
  if (impute == "avg") {
    mu <- colMeans(geno, na.rm = TRUE)
    na <- which(is.na(geno), arr.ind = TRUE)
    geno[na] <- mu[na[, 2]]
  }
  list(geno = geno, bim = bim[cols, ], fam = fam)
}

#' Select cis SNPs around a gene
#'
#' @param bim Data frame with columns `chr`, `snp`, `pos` (e.g. the `bim` returned by [read_plink()]).
#' @param chr,start,end Gene chromosome and coordinates.
#' @param window Window (bp) added on each side of the gene body (default 500 kb).
#' @return Character vector of SNP IDs.
#' @export
select_cis_snps <- function(bim, chr, start, end, window = 5e5) {
  bim$snp[as.character(bim$chr) == as.character(chr) &
            bim$pos >= start - window & bim$pos <= end + window]
}

#' Read / write a GCTA binary GRM
#'
#' @param prefix Path prefix of the `.grm.bin`/`.grm.id` files.
#' @return `read_grm_bin()` returns a full symmetric matrix with donor IDs
#'   (second column of `.grm.id`) as dimnames.
#' @export
read_grm_bin <- function(prefix) {
  id <- utils::read.table(paste0(prefix, ".grm.id"), stringsAsFactors = FALSE)
  n <- nrow(id)
  con <- file(paste0(prefix, ".grm.bin"), "rb")
  on.exit(close(con))
  v <- readBin(con, what = numeric(), n = n * (n + 1) / 2, size = 4)
  K <- matrix(0, n, n)
  K[lower.tri(K, diag = TRUE)] <- v
  K[upper.tri(K)] <- t(K)[upper.tri(K)]
  dimnames(K) <- list(id[[2]], id[[2]])
  K
}

#' @rdname read_grm_bin
#' @param K Symmetric GRM with donor IDs as dimnames.
#' @param n_snps Number of SNPs to record in the `.grm.N.bin` file.
#' @export
write_grm_bin <- function(K, prefix, n_snps = 1L) {
  stopifnot(is.matrix(K), nrow(K) == ncol(K))
  vals <- K[lower.tri(K, diag = TRUE)]
  con <- file(paste0(prefix, ".grm.bin"), "wb")
  writeBin(as.numeric(vals), con, size = 4)
  close(con)
  con <- file(paste0(prefix, ".grm.N.bin"), "wb")
  writeBin(rep(as.numeric(n_snps), length(vals)), con, size = 4)
  close(con)
  ids <- rownames(K)
  utils::write.table(data.frame(ids, ids), paste0(prefix, ".grm.id"),
                     quote = FALSE, row.names = FALSE, col.names = FALSE)
  invisible(prefix)
}
