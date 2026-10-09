#' Perturbation heritability with cell-state or cell-type interaction
#'
#' The perturb-seq analogue of [schint_cell()]: the donor/genotype kernel is
#' replaced by "same perturbation" (cells carrying the same perturbation are
#' correlated, cells with different perturbations are not). Haseman-Elston
#' regression on pairs of cells partitions the variance of each gene into a
#' perturbation component `P`, a perturbation-by-context component `P x context`
#' (context = continuous cell state and/or categorical cell type), context main
#' effects and covariates.
#'
#' For cells `a`, `b` the regressors are `1`, `1(pert_a == pert_b)`,
#' `u_a u_b` and `1(pert_a == pert_b) u_a u_b`; for a categorical context
#' `u_a u_b` is the indicator that both cells are in the same level (so
#' `P x celltype` is "same perturbation and same cell type").
#' Several genes can be fitted at once; the design-dependent `X'X` is computed
#' only once.
#'
#' @param data Data frame with one row per cell.
#' @param y Character vector of expression column(s) (one per gene).
#' @param perturb Name of the column holding the perturbation (target) of each cell.
#' @param context Column name(s) of the cell context interacting with the
#'   perturbation: numeric = cell state, factor/character = cell type. `NULL` for
#'   no interaction.
#' @param covariates Other covariates (main effects only).
#' @param perturb_group Optional: a column name giving a perturbation group for each
#'   perturbation, or a named list of perturbation vectors. Each group then gets its
#'   own `P` and `P x context` components.
#' @param control Perturbation labels to drop (e.g. non-targeting controls). `NA` and
#'   `""` are always dropped.
#' @param min_cells Minimum cells per perturbation (perturbations with fewer are dropped).
#' @param jackknife Delete-perturbation-block jackknife standard errors.
#' @param n_blocks Number of jackknife blocks over perturbations (default one per perturbation).
#' @inheritParams schint_cell
#' @return An object of class `"schint_pert"` with long-format `coefficients` and
#'   `summary` tables (one block of rows per gene).
#' @examples
#' data(schint_pert_example)
#' fit <- schint_pert(schint_pert_example, y = c("GENE_A", "GENE_B"),
#'                    perturb = "perturb", context = c("state", "celltype"),
#'                    control = "NT", covariates = c("pc1", "pc2"))
#' fit
#' @export
schint_pert <- function(data, y, perturb, context = NULL, covariates = NULL,
                        perturb_group = NULL, control = NULL, min_cells = 2,
                        cat_mode = c("pooled", "per_level"), context_main = TRUE,
                        jackknife = FALSE, n_blocks = NULL, scale_y = TRUE, scale_x = TRUE,
                        seed = 1) {
  cat_mode <- match.arg(cat_mode)
  data <- as.data.frame(data)
  group_col <- if (is.character(perturb_group) && length(perturb_group) == 1L) perturb_group else NULL
  need <- unique(c(y, perturb, context, covariates, group_col))
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("Columns not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  covariates <- setdiff(covariates, context)

  pert <- as.character(data[[perturb]])
  ok <- stats::complete.cases(data[, setdiff(need, y), drop = FALSE]) & !is.na(pert) & pert != "" &
    !(pert %in% control)
  data <- data[ok, , drop = FALSE]
  pert <- pert[ok]
  cnt <- table(pert)
  big <- names(cnt)[cnt >= min_cells]
  data <- data[pert %in% big, , drop = FALSE]
  pert <- pert[pert %in% big]
  if (length(big) < 3L) stop("Fewer than 3 perturbations remain after filtering.", call. = FALSE)
  units <- sort(unique(pert))
  n <- length(units)
  didx <- match(pert, units)
  N <- nrow(data)

  ## perturbation groups
  grp_of <- rep("P", n)
  if (!is.null(group_col)) {
    gg <- tapply(as.character(data[[group_col]]), pert, function(v) unique(v))
    if (any(lengths(gg) != 1L)) stop("Each perturbation must belong to exactly one group.", call. = FALSE)
    grp_of <- unlist(gg[units])
  } else if (is.list(perturb_group)) {
    map <- stats::setNames(rep(names(perturb_group), lengths(perturb_group)), unlist(perturb_group))
    grp_of <- unname(map[units])
    if (anyNA(grp_of)) {
      message(sum(is.na(grp_of)), " perturbations not in `perturb_group` dropped.")
      return(schint_pert(data[!is.na(map[pert]), ], y, perturb, context, covariates, perturb_group,
                         control, min_cells, cat_mode, context_main, jackknife, n_blocks,
                         scale_y, scale_x, seed))
    }
  }
  groups <- sort(unique(grp_of))
  pname <- if (length(groups) == 1L && is.null(group_col) && !is.list(perturb_group)) "P" else paste0("P_", groups)
  kern <- c(list(one = .k_one()),
            stats::setNames(lapply(groups, function(g) .k_diag(as.numeric(grp_of == g))), pname))

  ones <- matrix(1, N, 1L)
  comps <- list(list(name = "Intercept", group = "intercept", kernel = "one", U = ones,
                     fid = NA_integer_, pgrp = NA_character_, cvar = NA_character_))
  add <- function(name, group, kernel, U, fid = NA_integer_, pgrp = NA_character_, cvar = NA_character_) {
    comps[[length(comps) + 1L]] <<- list(name = name, group = group, kernel = kernel, U = U,
                                         fid = fid, pgrp = pgrp, cvar = cvar)
  }
  for (g in pname) add(g, "P", g, ones, pgrp = g)
  ctx <- .design_terms(data, context, cat_mode, scale_x)
  for (tm in ctx$terms) {
    if (context_main) add(tm$name, "context", "one", tm$U, tm$fid)
    for (g in pname) add(paste0(g, ":", tm$name), "PxC", g, tm$U, tm$fid, g, tm$var)
  }
  cov <- .design_terms(data, covariates, "pooled", scale_x, fid_start = ctx$fid)
  for (tm in cov$terms) add(tm$name, "covariate", "one", tm$U, tm$fid)
  nms <- vapply(comps, function(cp) cp$name, "")
  if (anyDuplicated(nms)) stop("Duplicated component names.", call. = FALSE)
  p <- length(comps)

  bid <- NULL
  if (jackknife) {
    nb <- if (is.null(n_blocks)) n else min(as.integer(n_blocks), n)
    if (nb < 3L) stop("`n_blocks` must be at least 3.", call. = FALSE)
    if (!is.null(seed)) set.seed(seed)
    bid <- if (nb == n) seq_len(n) else sample(rep_len(seq_len(nb), n))
  }

  xx <- .he_xtx(comps, kern, didx, bid)
  group <- vapply(comps, function(cp) cp$group, "")
  pg <- vapply(comps, function(cp) cp$pgrp, "")
  mult <- vapply(comps, function(cp) {
    if (cp$group == "intercept") return(NA_real_)
    kd <- if (cp$kernel == "one") 1 else kern[[cp$kernel]]$d[didx]
    mean(kd * rowSums(cp$U^2))
  }, 0)

  terms <- list(P = which(group == "P"))
  if (any(group == "PxC")) {
    terms$PxC <- which(group == "PxC")
    terms$P_total <- which(group %in% c("P", "PxC"))
    cv <- vapply(comps, function(cp) cp$cvar, "")
    if (length(context) > 1L) for (v in context) terms[[paste0("PxC:", v)]] <- which(group == "PxC" & cv == v)
  }
  if (length(pname) > 1L) for (g in pname) {
    terms[[g]] <- which(group == "P" & pg == g)
    if (any(group == "PxC")) {
      terms[[paste0(g, ":C")]] <- which(group == "PxC" & pg == g)
      terms[[paste0(g, "_total")]] <- which(pg == g)
    }
  }
  if (any(group == "context")) terms$context <- which(group == "context")
  if (any(group == "covariate")) terms$covariate <- which(group == "covariate")
  terms$total_explained <- which(group != "intercept")
  A <- do.call(rbind, lapply(terms, function(ix) { a <- numeric(p); a[ix] <- 1; a }))
  rownames(A) <- names(terms)

  coef_l <- list(); summ_l <- list(); jk_l <- list()
  for (gene in y) {
    yy <- data[[gene]]
    if (anyNA(yy) || !all(is.finite(yy)) || stats::sd(yy) == 0) {
      warning("Skipping gene ", gene, " (missing/non-finite values or zero variance).", call. = FALSE)
      next
    }
    yy <- if (scale_y) as.numeric(scale(yy)) else yy - mean(yy)
    xy <- .he_xty(comps, kern, didx, yy, bid)
    beta <- .solve_he(xx$XtX, xy$Xty, warn = FALSE)
    ve <- beta * mult
    vy <- stats::var(yy)
    cf <- data.frame(gene = gene, component = nms, group = group, estimate = beta,
                     variance = ve, fraction = ve / vy, stringsAsFactors = FALSE)
    sm <- data.frame(gene = gene, term = rownames(A),
                     variance = as.numeric(A %*% replace(ve, is.na(ve), 0)), stringsAsFactors = FALSE)
    sm$fraction <- sm$variance / vy
    if (jackknife) {
      nb <- max(bid)
      Bm <- t(vapply(seq_len(nb), function(b) .solve_he(xx$XtXb[b, , ], xy$Xtyb[b, ], warn = FALSE),
                     numeric(p)))
      Bm <- Bm[stats::complete.cases(Bm), , drop = FALSE]
      ng <- nrow(Bm)
      jse <- function(M) sqrt((ng - 1) / ng * colSums(sweep(M, 2, colMeans(M))^2))
      cf$se <- jse(Bm)
      cf$z <- cf$estimate / cf$se
      cf$p <- 2 * stats::pnorm(-abs(cf$z))
      V <- sweep(Bm, 2, mult, "*"); V[is.na(V)] <- 0
      sm$se <- jse(V %*% t(A))
      sm$se_fraction <- sm$se / vy
      jk_l[[gene]] <- Bm
    }
    coef_l[[gene]] <- cf
    summ_l[[gene]] <- sm
  }
  if (!length(coef_l)) stop("No gene could be fitted.", call. = FALSE)
  structure(list(
    coefficients = do.call(rbind, coef_l), summary = do.call(rbind, summ_l),
    n_obs = N, n_perturb = n, perturb_groups = groups, genes = names(coef_l),
    context = context, covariates = covariates,
    XtX = `dimnames<-`(xx$XtX, list(nms, nms)),
    jackknife = if (jackknife) list(estimates = jk_l, n_blocks = max(bid)) else NULL,
    call = match.call()
  ), class = "schint_pert")
}

#' @export
print.schint_pert <- function(x, digits = 3, max_genes = 5, ...) {
  cat("scHINT perturbation model: ", x$n_obs, " cells, ", x$n_perturb, " perturbations, ",
      length(x$genes), " gene", if (length(x$genes) > 1) "s", "\n\n", sep = "")
  s <- x$summary
  s <- s[s$gene %in% utils::head(x$genes, max_genes), ]
  s <- .compact_summary(s, digits)
  print(s, row.names = FALSE, na.print = "")
  if (length(x$genes) > max_genes) cat("... ", length(x$genes) - max_genes, " more genes in $summary\n", sep = "")
  if (!is.null(x$jackknife)) cat("\nStandard errors: jackknife over ", x$jackknife$n_blocks, " perturbation blocks\n", sep = "")
  invisible(x)
}

#' @export
summary.schint_pert <- function(object, ...) object$coefficients
