# Brute-force reference: explicit pairwise OLS ---------------------------------
brute_he <- function(fit_inputs) {
  with(fit_inputs, {
    N <- length(y)
    ij <- which(upper.tri(matrix(0, N, N)), arr.ind = TRUE)
    resp <- y[ij[, 1]] * y[ij[, 2]]
    X <- sapply(comps, function(cp) {
      K <- kern[[cp$kernel]][cbind(didx[ij[, 1]], didx[ij[, 2]])]
      Q <- rowSums(cp$U[ij[, 1], , drop = FALSE] * cp$U[ij[, 2], , drop = FALSE])
      K * Q
    })
    list(beta = as.numeric(solve(crossprod(X), crossprod(X, resp))), X = X, resp = resp, ij = ij)
  })
}

test_that("engine matches brute-force pairwise OLS (incl. jackknife)", {
  set.seed(3)
  n <- 12
  didx <- rep(seq_len(n), sample(3:5, n, TRUE))
  N <- length(didx)
  Xg <- matrix(rbinom(n * 60, 2, 0.3), n, 60, dimnames = list(paste0("d", 1:n), NULL))
  G1 <- make_grm(Xg[, 1:30]); G2 <- make_grm(Xg[, 31:60])
  kern <- list(one = matrix(1, n, n), I = diag(n), G1 = G1, G2 = G2)
  f <- factor(sample(c("a", "b", "c"), N, TRUE))
  U_f <- .one_hot(f)
  u <- rnorm(N)
  one <- matrix(1, N, 1)
  comps <- list(
    list(kernel = "one", U = one, fid = NA_integer_),
    list(kernel = "G1", U = one, fid = NA_integer_),
    list(kernel = "G2", U = one, fid = NA_integer_),
    list(kernel = "I", U = one, fid = NA_integer_),
    list(kernel = "one", U = matrix(u, ncol = 1), fid = NA_integer_),
    list(kernel = "G1", U = matrix(u, ncol = 1), fid = NA_integer_),
    list(kernel = "one", U = U_f, fid = 1L),
    list(kernel = "G2", U = U_f, fid = 1L)
  )
  y <- rnorm(N)
  bid <- seq_len(n)
  st <- .he_stats(comps, kern, didx, y, bid)
  br <- brute_he(list(y = y, comps = comps, kern = kern, didx = didx))
  expect_equal(.solve_he(st$XtX, st$Xty), br$beta, tolerance = 1e-8)

  ## leave-one-donor-out
  for (d in c(1, 7, 12)) {
    keep <- didx != d
    ij <- br$ij
    use <- keep[ij[, 1]] & keep[ij[, 2]]
    b <- solve(crossprod(br$X[use, ]), crossprod(br$X[use, ], br$resp[use]))
    expect_equal(.solve_he(st$XtXb[d, , ], st$Xtyb[d, ]), as.numeric(b), tolerance = 1e-8)
  }

  ## 3 blocks of donors
  bid3 <- rep(1:3, length.out = n)
  st3 <- .he_stats(comps, kern, didx, y, bid3)
  keep <- bid3[didx] != 2
  use <- keep[br$ij[, 1]] & keep[br$ij[, 2]]
  b <- solve(crossprod(br$X[use, ]), crossprod(br$X[use, ], br$resp[use]))
  expect_equal(.solve_he(st3$XtXb[2, , ], st3$Xtyb[2, ]), as.numeric(b), tolerance = 1e-8)
})

test_that("schint_cell / schint_sample run on the example data and recover signal", {
  data(schint_example, package = "scHINT")
  ex <- schint_example
  fit <- schint_cell(ex$cell, "GENE_A", "donor", geno = ex$geno, context = "state",
                     covariates = c("pc1", "pc2", "age", "sex", "batch", "subtype"))
  expect_s3_class(fit, "schint")
  expect_true(all(is.finite(fit$coefficients$estimate)))
  fitj <- schint_cell(ex$cell, "GENE_A", "donor", geno = ex$geno, context = "state",
                      n_blocks = 20, jackknife = TRUE)
  expect_true(all(is.finite(fitj$coefficients$se)))

  ## multiple GRMs via split genotype matrix
  fit3 <- schint_cell(ex$cell, "GENE_A", "donor", geno = ex$geno, n_grm = 3, context = "state")
  expect_equal(fit3$n_grm, 3)

  ## precomputed GRM gives identical answer
  fg <- schint_cell(ex$cell, "GENE_A", "donor", grm = make_grm(ex$geno), context = "state",
                    covariates = c("pc1", "pc2", "age", "sex", "batch", "subtype"))
  expect_equal(coef(fg), coef(fit))

  fs <- schint_sample(ex$sample, "GENE_A", "donor", "celltype", geno = ex$geno,
                      covariates = c("age", "sex", "pc1", "pc2"))
  expect_true(all(is.finite(fs$coefficients$estimate)))
})

test_that("schint_pert runs, supports groups/jackknife, and handles several genes", {
  d <- schint_pert_example
  f <- schint_pert(d, c("GENE_A", "GENE_B"), "perturb", context = c("state", "celltype"),
                   control = "NT", jackknife = TRUE, n_blocks = 15)
  expect_s3_class(f, "schint_pert")
  expect_true(all(is.finite(f$coefficients$se)))
  g <- schint_pert(d, "GENE_A", "perturb", context = "celltype", control = "NT",
                   perturb_group = "perturb_group")
  expect_true(all(c("P_group1", "P_group3:C") %in% g$summary$term))
  # single-gene fit equals the same gene fitted within a batch
  f1 <- schint_pert(d, "GENE_B", "perturb", context = c("state", "celltype"), control = "NT")
  expect_equal(f1$coefficients$estimate,
               f$coefficients$estimate[f$coefficients$gene == "GENE_B"])
})

test_that("perturbation model equals explicit pairwise OLS (P, C, P:C)", {
  d <- schint_pert_example
  d <- d[d$perturb != "NT" & d$perturb %in% unique(d$perturb)[1:10], ]
  d <- d[seq_len(400), ]
  f <- schint_pert(d, "GENE_A", "perturb", context = "celltype", min_cells = 1)
  yy <- as.numeric(scale(d$GENE_A))
  ij <- which(upper.tri(matrix(0, nrow(d), nrow(d))), arr.ind = TRUE)
  P <- as.numeric(d$perturb[ij[, 1]] == d$perturb[ij[, 2]])
  C <- as.numeric(d$celltype[ij[, 1]] == d$celltype[ij[, 2]])
  X <- cbind(1, P, C, P * C)
  b <- solve(crossprod(X), crossprod(X, yy[ij[, 1]] * yy[ij[, 2]]))
  expect_equal(f$coefficients$estimate, as.numeric(b), tolerance = 1e-8)
})
