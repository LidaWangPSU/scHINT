# Internal HE-regression engine -------------------------------------------
#
# All scHINT models are Haseman-Elston regressions on pairs of observations
# (cells, or sample x cell-type units). For a pair (a < b) the response is
# y_a * y_b and every variance component r contributes one regressor
#
#     X_r[a, b] = K_r[id_a, id_b] * Q_r[a, b],      Q_r = U_r U_r'
#
# where K_r is a donor-by-donor kernel (all-ones, identity, or a GRM) and U_r is
# an observation-by-q matrix (a column of 1s, a continuous covariate, or a
# one-hot encoding of a factor). The normal equations only need
#
#     sum_{a<b} X_r[a, b] * X_s[a, b]   and   sum_{a<b} X_r[a, b] * y_a y_b,
#
# and both reduce to quadratic forms on donor-aggregated vectors, so the cost
# is O(n_donors^2) per entry rather than O(n_obs^2). Leave-donor-block-out
# versions of every sum are obtained from the same quantities, which gives the
# jackknife at almost no extra cost.

# Row-wise Khatri-Rao product of two N x q matrices.
.row_kr <- function(A, B) {
  qa <- ncol(A)
  qb <- ncol(B)
  if (qa == 1L && qb == 1L) return(A * B)
  A[, rep(seq_len(qa), times = qb), drop = FALSE] *
    B[, rep(seq_len(qb), each = qa), drop = FALSE]
}

# Kernels are dense n x n matrices, or compact objects: .k_one() (all ones) and
# .k_diag(d) (diagonal; identity and "same perturbation" kernels). The compact forms
# keep memory O(n), which matters for screens with thousands of perturbations.
.k_one <- function() structure(list(type = "one"), class = "kern")
.k_diag <- function(d) structure(list(type = "diag", d = as.numeric(d)), class = "kern")
.kdiag <- function(M) if (inherits(M, "kern")) { if (M$type == "diag") M$d else NA } else diag(M)
.kmul <- function(M, W, ix = NULL) {
  if (!inherits(M, "kern")) return((if (is.null(ix)) M else M[ix, ix, drop = FALSE]) %*% W)
  if (M$type == "diag") return((if (is.null(ix)) M$d else M$d[ix]) * W)
  matrix(colSums(W), nrow(W), ncol(W), byrow = TRUE)
}
.kprod <- function(A, B) {
  ka <- inherits(A, "kern"); kb <- inherits(B, "kern")
  if (ka && A$type == "one") return(B)
  if (kb && B$type == "one") return(A)
  if (ka && kb) return(.k_diag(A$d * B$d))
  if (ka) return(.k_diag(A$d * diag(B)))
  if (kb) return(.k_diag(B$d * diag(A)))
  A * B
}

# sum_{a<b} M[id_a, id_b] * sum_c Z[a, c] Z[b, c], plus (optionally) the same
# quantity with every block of donors left out in turn.
#   Z    N x q matrix
#   M    n x n donor kernel
#   didx integer donor index of each row of Z (1..n, every donor present)
#   bid  NULL, or integer jackknife block of each donor (1..nb)
.pair_stat <- function(Z, M, didx, bid = NULL) {
  W <- rowsum(Z, didx, reorder = TRUE)            # n x q, donor-aggregated
  MW <- .kmul(M, W)
  per_don <- rowSums(W * MW)
  md <- if (inherits(M, "kern") && M$type == "one") rep(1, nrow(W)) else .kdiag(M)
  dg_cell <- md[didx] * rowSums(Z * Z)            # a == b terms
  total <- 0.5 * (sum(per_don) - sum(dg_cell))
  if (is.null(bid)) return(list(v = total, loo = NULL))

  nb <- max(bid)
  dg_don <- as.numeric(rowsum(dg_cell, didx, reorder = TRUE))
  s1 <- as.numeric(rowsum(per_don, bid, reorder = TRUE))
  self <- as.numeric(rowsum(md * rowSums(W * W), bid, reorder = TRUE))
  cross <- self
  big <- which(tabulate(bid, nbins = nb) > 1L)
  for (b in big) {                                # off-diagonal terms inside a block
    ix <- which(bid == b)
    Wb <- W[ix, , drop = FALSE]
    cross[b] <- cross[b] + sum(Wb * .kmul(M, Wb, ix)) -
      sum(md[ix] * rowSums(Wb * Wb))
  }
  Tb <- sum(per_don) - 2 * s1 + cross
  Dgb <- sum(dg_cell) - as.numeric(rowsum(dg_don, bid, reorder = TRUE))
  list(v = total, loo = 0.5 * (Tb - Dgb))
}

# X'y (and leave-block-out copies) for one response vector.
.he_xty <- function(comps, kern, didx, y, bid = NULL) {
  p <- length(comps)
  Xty <- numeric(p)
  Xtyb <- if (is.null(bid)) NULL else matrix(0, max(bid), p)
  for (r in seq_len(p)) {
    st <- .pair_stat(comps[[r]]$U * y, kern[[comps[[r]]$kernel]], didx, bid)
    Xty[r] <- st$v
    if (!is.null(bid)) Xtyb[, r] <- st$loo
  }
  list(Xty = Xty, Xtyb = Xtyb)
}

# X'X (and leave-block-out copies); independent of the response.
.he_xtx <- function(comps, kern, didx, bid = NULL) {
  p <- length(comps)
  nb <- if (is.null(bid)) 0L else max(bid)
  XtX <- matrix(0, p, p)
  if (nb > 0L) XtXb <- array(0, c(nb, p, p))
  kid <- vapply(comps, function(cp) cp$kernel, "")
  kpairs <- unique(t(apply(cbind(rep(kid, each = p), rep(kid, times = p)), 1L, sort)))
  for (kk in seq_len(nrow(kpairs))) {
    ka <- kpairs[kk, 1L]
    kb <- kpairs[kk, 2L]
    M <- .kprod(kern[[ka]], kern[[kb]])
    for (r in which(kid == ka | kid == kb)) {
      for (s in r:p) {
        ok <- (kid[r] == ka && kid[s] == kb) || (kid[r] == kb && kid[s] == ka)
        if (!ok) next
        cr <- comps[[r]]
        cs <- comps[[s]]
        Z <- if (!is.na(cr$fid) && !is.na(cs$fid) && cr$fid == cs$fid) {
          cr$U                                    # one-hot x same one-hot
        } else {
          .row_kr(cr$U, cs$U)
        }
        st <- .pair_stat(Z, M, didx, bid)
        XtX[r, s] <- XtX[s, r] <- st$v
        if (nb > 0L) {
          XtXb[, r, s] <- st$loo
          XtXb[, s, r] <- st$loo
        }
      }
    }
  }
  if (nb > 0L) list(XtX = XtX, XtXb = XtXb) else list(XtX = XtX)
}

.he_stats <- function(comps, kern, didx, y, bid = NULL) {
  c(.he_xtx(comps, kern, didx, bid), .he_xty(comps, kern, didx, y, bid))
}

.solve_he <- function(XtX, Xty, warn = TRUE) {
  beta <- tryCatch(solve(XtX, Xty), error = function(e) NULL)
  if (is.null(beta)) {
    if (warn) warning("X'X is singular; adding a small ridge (check for collinear components).",
                      call. = FALSE)
    ridge <- 1e-8 * mean(diag(XtX))
    beta <- tryCatch(solve(XtX + diag(ridge, nrow(XtX)), Xty),
                     error = function(e) rep(NA_real_, length(Xty)))
  }
  as.numeric(beta)
}
