#' Simulate cell-level and sample-level expression for the example data
#'
#' Cell level: one cell type, `cells` per donor, a continuous cell state, donor
#' and cell covariates. Sample level: `n_celltypes` cell types per donor.
#' Each gene has its own variance architecture (`genes`).
#'
#' @param geno Donors x SNPs dosage matrix (row names = donor IDs).
#' @param genes Named list; each element is a list with `g` (main genetic
#'   variance), `gxc` (genotype x state / cell-type-specific genetic variance),
#'   `i` (donor variance), `cov` (covariate variance); the remainder is noise.
#' @param n_causal Number of causal SNPs per gene.
#' @param cells Range of cells per donor.
#' @param n_celltypes Number of cell types for the sample-level data.
#' @param seed Random seed.
#' @return A list with `cell`, `sample` and `truth` (the generating variances).
#' @export
simulate_schint <- function(geno,
                            genes = list(GENE_A = list(g = 0.20, gxc = 0.15, i = 0.10, cov = 0.05),
                                         GENE_B = list(g = 0.25, gxc = 0.00, i = 0.10, cov = 0.05),
                                         GENE_C = list(g = 0.00, gxc = 0.00, i = 0.15, cov = 0.05)),
                            n_causal = 30, cells = c(15, 45), n_celltypes = 5, seed = 1) {
  set.seed(seed)
  n <- nrow(geno)
  donors <- rownames(geno)
  Xs <- scale(apply(geno, 2, function(x) { x[is.na(x)] <- mean(x, na.rm = TRUE); x }))
  Xs[is.na(Xs)] <- 0
  gscore <- function() {
    b <- numeric(ncol(Xs))
    b[sample(ncol(Xs), min(n_causal, ncol(Xs)))] <- stats::rnorm(min(n_causal, ncol(Xs)))
    as.numeric(scale(Xs %*% b))
  }
  z <- function(x) as.numeric(scale(x))

  ## donor covariates
  don <- data.frame(donor = donors, age = round(stats::runif(n, 25, 80)),
                    sex = factor(sample(c("F", "M"), n, TRUE)),
                    batch = factor(sample(paste0("B", 1:4), n, TRUE)),
                    stringsAsFactors = FALSE)
  don_cov <- z(don$age) + z(as.integer(don$sex)) + z(as.integer(don$batch))

  ## ---- cell level -------------------------------------------------------------
  nc <- sample(cells[1]:cells[2], n, TRUE)
  cell <- don[rep(seq_len(n), nc), , drop = FALSE]
  N <- nrow(cell)
  shift <- stats::rnorm(n, 0, 0.3)[match(cell$donor, donors)]
  cell$state <- z(shift + stats::rnorm(N))
  cell$subtype <- factor(ifelse(cell$state + stats::rnorm(N, 0, 0.5) > 0, "B", "A"))
  cell$pc1 <- stats::rnorm(N)
  cell$pc2 <- stats::rnorm(N)
  cell$cell_id <- sprintf("cell%05d", seq_len(N))
  ci <- match(cell$donor, donors)
  cell_cov <- z(cell$pc1) + z(cell$pc2) + z(as.integer(cell$subtype)) + z(don_cov[ci])

  ## ---- sample level -------------------------------------------------------------
  cts <- paste0("CT", seq_len(n_celltypes))
  samp <- merge(data.frame(donor = donors, stringsAsFactors = FALSE),
                data.frame(celltype = factor(cts, levels = cts)), by = NULL)
  samp <- merge(samp, don, by = "donor", sort = FALSE)
  samp <- samp[order(match(samp$donor, donors), samp$celltype), ]
  M <- nrow(samp)
  samp$pc1 <- stats::rnorm(M)
  samp$pc2 <- stats::rnorm(M)
  si <- match(samp$donor, donors)
  ct_i <- as.integer(samp$celltype)
  samp_cov <- z(samp$pc1) + z(samp$pc2) + z(don_cov[si])

  truth <- do.call(rbind, lapply(names(genes), function(gn) {
    p <- genes[[gn]]
    data.frame(gene = gn, g = p$g, gxc = p$gxc, i = p$i, cov = p$cov,
               noise = 1 - p$g - p$gxc - p$i - p$cov)
  }))

  for (gn in names(genes)) {
    p <- genes[[gn]]
    noise <- 1 - p$g - p$gxc - p$i - p$cov
    if (noise <= 0) stop("Variances of gene ", gn, " must sum to < 1.")
    ## cell level: G×Cell-State
    g1 <- gscore(); g2 <- gscore(); di <- stats::rnorm(n)
    cell[[gn]] <- sqrt(p$g) * g1[ci] + sqrt(p$gxc) * g2[ci] * cell$state +
      sqrt(p$i) * z(di)[ci] + sqrt(p$cov) * z(cell_cov) + sqrt(noise) * stats::rnorm(N)
    ## sample level: shared G + cell-type-specific G
    g0 <- gscore()
    gct <- vapply(seq_len(n_celltypes), function(k) gscore(), numeric(n))
    di <- z(stats::rnorm(n))
    mu_ct <- stats::rnorm(n_celltypes)
    samp[[gn]] <- sqrt(p$g) * g0[si] + sqrt(p$gxc) * gct[cbind(si, ct_i)] +
      sqrt(p$i) * di[si] + sqrt(p$cov) * z(samp_cov) + sqrt(noise) * stats::rnorm(M) + mu_ct[ct_i]
  }
  rownames(cell) <- rownames(samp) <- NULL
  list(cell = cell[, c("cell_id", "donor", "state", "subtype", "pc1", "pc2", "age", "sex", "batch", names(genes))],
       sample = samp[, c("donor", "celltype", "pc1", "pc2", "age", "sex", "batch", names(genes))],
       truth = truth)
}

#' Simulate perturb-seq-like data for the perturbation heritability model
#'
#' Cells carry one of `n_perturb` perturbations (in three groups) or a non-targeting
#' control (`"NT"`), have a continuous cell `state` and one of three cell types. Gene
#' expression has perturbation (`p`), perturbation x state (`pxs`) and
#' perturbation x cell-type (`pxc`) variance, scaled per perturbation group by
#' `group_scale`.
#'
#' @param genes Named list of variances `p`, `pxs`, `pxc`, `cov` (rest is noise).
#' @param n_perturb Number of perturbations.
#' @param cells Range of cells per perturbation.
#' @param n_control Number of control cells.
#' @param group_scale Multiplier on the perturbation variances in each of the three groups.
#' @param seed Random seed.
#' @return A data frame with one row per cell.
#' @export
simulate_pert <- function(genes = list(GENE_A = list(p = 0.10, pxs = 0.10, pxc = 0.00, cov = 0.05),
                                       GENE_B = list(p = 0.10, pxs = 0.00, pxc = 0.10, cov = 0.05),
                                       GENE_C = list(p = 0.15, pxs = 0.00, pxc = 0.00, cov = 0.05),
                                       GENE_D = list(p = 0.00, pxs = 0.00, pxc = 0.00, cov = 0.05)),
                          n_perturb = 60, cells = c(30, 80), n_control = 600,
                          group_scale = c(1, 0.5, 0), seed = 1) {
  set.seed(seed)
  z <- function(x) as.numeric(scale(x))
  perts <- sprintf("PERT%02d", seq_len(n_perturb))
  grp <- paste0("group", rep_len(seq_along(group_scale), n_perturb))
  nc <- sample(cells[1]:cells[2], n_perturb, TRUE)
  df <- data.frame(perturb = c(rep(perts, nc), rep("NT", n_control)), stringsAsFactors = FALSE)
  N <- nrow(df)
  pi <- match(df$perturb, perts)                    # NA for controls
  df$perturb_group <- ifelse(is.na(pi), "control", grp[pi])
  df$celltype <- factor(sample(c("CT1", "CT2", "CT3"), N, TRUE, prob = c(.5, .3, .2)))
  df$state <- z(stats::rnorm(N) + 0.5 * as.integer(df$celltype))
  df$pc1 <- stats::rnorm(N)
  df$pc2 <- stats::rnorm(N)
  sc <- c(group_scale, 0)[ifelse(is.na(pi), length(group_scale) + 1L, match(grp[pi], paste0("group", seq_along(group_scale))))]
  ct <- as.integer(df$celltype)
  for (gn in names(genes)) {
    p <- genes[[gn]]
    noise <- 1 - p$p - p$pxs - p$pxc - p$cov
    if (noise <= 0) stop("Variances of ", gn, " must sum to < 1.")
    ap <- stats::rnorm(n_perturb); bp <- stats::rnorm(n_perturb)
    cp <- matrix(stats::rnorm(n_perturb * 3), n_perturb, 3)
    idx <- ifelse(is.na(pi), 1L, pi)
    eff <- sqrt(p$p) * ap[idx] + sqrt(p$pxs) * bp[idx] * df$state + sqrt(p$pxc) * cp[cbind(idx, ct)]
    df[[gn]] <- sqrt(sc) * eff + sqrt(p$cov) * z(df$pc1 + df$pc2) + sqrt(noise) * stats::rnorm(N) +
      c(0, 0.5, -0.3)[ct]
  }
  df
}

#' Simulate genotypes with local LD
#'
#' Dosages (0/1/2) from a Gaussian copula with AR(1) correlation between
#' neighbouring SNPs and uniform minor-allele frequencies.
#'
#' @param n Number of donors.
#' @param m Number of SNPs.
#' @param rho Correlation of adjacent SNPs.
#' @param maf Range of minor-allele frequencies.
#' @param seed Random seed.
#' @return Integer matrix with donor IDs `D001...` and SNP IDs `SNP001...`.
#' @export
simulate_geno <- function(n = 400, m = 250, rho = 0.6, maf = c(0.05, 0.5), seed = 1) {
  set.seed(seed)
  f <- stats::runif(m, maf[1], maf[2])
  geno <- matrix(0L, n, m)
  for (h in 1:2) {                                   # two haplotypes per donor
    z <- matrix(stats::rnorm(n * m), n, m)
    for (j in 2:m) z[, j] <- rho * z[, j - 1] + sqrt(1 - rho^2) * z[, j]
    geno <- geno + (sweep(z, 2, stats::qnorm(1 - f), ">") * 1L)
  }
  dimnames(geno) <- list(sprintf("D%03d", seq_len(n)), sprintf("SNP%03d", seq_len(m)))
  geno
}
