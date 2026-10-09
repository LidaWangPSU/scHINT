#' Build a genomic relationship matrix from genotypes
#'
#' @param geno Numeric matrix of allele dosages (0/1/2), donors in rows (with
#'   row names) and SNPs in columns. `NA` is mean-imputed.
#' @param standardize If `TRUE` (default) each SNP is centred and scaled by
#'   `sqrt(2p(1-p))`; monomorphic SNPs are dropped. If `FALSE`, `geno` is
#'   assumed to be standardized already.
#' @return A donor-by-donor matrix `X X' / m`.
#' @export
make_grm <- function(geno, standardize = TRUE) {
  geno <- as.matrix(geno)
  if (is.null(rownames(geno))) stop("`geno` needs donor IDs as row names.", call. = FALSE)
  if (standardize) {
    mu <- colMeans(geno, na.rm = TRUE)
    p <- mu / 2
    sdv <- sqrt(2 * p * (1 - p))
    keep <- is.finite(sdv) & sdv > 0
    geno <- geno[, keep, drop = FALSE]
    mu <- mu[keep]
    sdv <- sdv[keep]
    geno <- sweep(sweep(geno, 2, mu, "-"), 2, sdv, "/")
  }
  geno[is.na(geno)] <- 0
  K <- tcrossprod(geno) / ncol(geno)
  dimnames(K) <- list(rownames(geno), rownames(geno))
  K
}

#' Split SNPs into groups (for multiple GRMs)
#'
#' @param geno Genotype matrix (donors x SNPs).
#' @param n Number of groups.
#' @param by `"maf"` (equal-count minor-allele-frequency bins, rarest first),
#'   `"position"` (contiguous blocks in column order, i.e. genomic order if
#'   columns are sorted) or `"random"`.
#' @param seed Seed used when `by = "random"`.
#' @return A list of `n` genotype matrices; use `lapply(split_snps(geno, 3), make_grm)` to
#'   obtain a list of GRMs for the `grm` argument of the model functions.
#' @export
split_snps <- function(geno, n, by = c("maf", "position", "random"), seed = 1) {
  by <- match.arg(by)
  m <- ncol(geno)
  if (n < 1 || n > m) stop("`n` must be between 1 and the number of SNPs.", call. = FALSE)
  grp <- switch(by,
    maf = {
      f <- colMeans(geno, na.rm = TRUE) / 2
      maf <- pmin(f, 1 - f)
      cut(rank(maf, ties.method = "first"), n, labels = FALSE)
    },
    position = cut(seq_len(m), n, labels = FALSE),
    random = {
      set.seed(seed)
      sample(rep_len(seq_len(n), m))
    }
  )
  lapply(split(seq_len(m), grp), function(ix) geno[, ix, drop = FALSE])
}

# Resolve the genetic input of a model into a named list of GRMs: one GRM built from
# all SNPs in `geno`, or the GRM(s) supplied in `grm`.
.resolve_grms <- function(geno, grm) {
  if (is.null(geno) == is.null(grm)) {
    stop("Provide exactly one of `geno` (genotypes) or `grm` (pre-computed GRM(s)).", call. = FALSE)
  }
  if (!is.null(grm)) {
    glist <- if (is.matrix(grm)) list(grm) else grm
    if (!is.list(glist) || !all(vapply(glist, is.matrix, TRUE))) {
      stop("`grm` must be a matrix or a list of matrices.", call. = FALSE)
    }
    for (g in glist) {
      if (nrow(g) != ncol(g) || is.null(rownames(g))) {
        stop("Each GRM must be square with donor IDs as row names.", call. = FALSE)
      }
    }
  } else {
    if (!(is.matrix(geno) || is.data.frame(geno))) {
      stop("`geno` must be a donors x SNPs matrix. For several GRMs pass a list of GRMs in `grm`.", call. = FALSE)
    }
    glist <- list(make_grm(as.matrix(geno)))
  }
  nm <- names(glist)
  if (is.null(nm)) nm <- if (length(glist) == 1L) "G" else paste0("G", seq_along(glist))
  nm[nm == ""] <- paste0("G", seq_along(glist))[nm == ""]
  names(glist) <- nm
  glist
}
