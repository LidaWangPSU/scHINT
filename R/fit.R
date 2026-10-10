# Shared fitting routine behind schint_cell() and schint_sample() -----------

.one_hot <- function(f) {
  f <- droplevels(as.factor(f))
  U <- matrix(0, length(f), nlevels(f), dimnames = list(NULL, levels(f)))
  U[cbind(seq_along(f), as.integer(f))] <- 1
  U
}

.is_cat <- function(x) is.factor(x) || is.character(x) || is.logical(x)

# Turn a data column into one or more components' design.
# type "continuous" -> one column; type "categorical" -> one-hot (pooled) or
# one component per level (per_level).
.design_terms <- function(data, vars, cat_mode, scale_x, fid_start = 0L) {
  out <- list()
  fid <- fid_start
  for (v in vars) {
    x <- data[[v]]
    if (.is_cat(x)) {
      U <- .one_hot(x)
      if (ncol(U) < 2L) stop("`", v, "` has a single level.", call. = FALSE)
      if (cat_mode == "pooled") {
        fid <- fid + 1L
        out[[v]] <- list(name = v, U = U, fid = fid, var = v)
      } else {
        for (l in colnames(U)) {
          nm <- paste0(v, "=", l)
          out[[nm]] <- list(name = nm, U = U[, l, drop = FALSE], fid = NA_integer_, var = v)
        }
      }
    } else {
      if (!is.numeric(x)) stop("`", v, "` must be numeric, factor or character.", call. = FALSE)
      x <- if (scale_x) as.numeric(scale(x)) else as.numeric(x)
      out[[v]] <- list(name = v, U = matrix(x, ncol = 1L), fid = NA_integer_, var = v)
    }
  }
  list(terms = out, fid = fid)
}

.schint_fit <- function(data, y, id, context, covariates, geno, grm,
                        cat_mode, gxc, ind_effect, context_main, scale_x,
                        jackknife, n_blocks, seed, call) {
  data <- as.data.frame(data)
  need <- unique(c(y, id, context, covariates))
  miss <- setdiff(need, names(data))
  if (length(miss)) stop("Columns not found in `data`: ", paste(miss, collapse = ", "), call. = FALSE)
  covariates <- setdiff(covariates, context)           # context main terms are added automatically
  if (!is.numeric(data[[y]])) stop("`y` must be a numeric column.", call. = FALSE)

  glist <- .resolve_grms(geno, grm)
  gnames <- names(glist)

  ## align observations and donors
  ok <- stats::complete.cases(data[, need, drop = FALSE]) & is.finite(data[[y]])
  if (any(!ok)) message(sum(!ok), " observations dropped for missing values.")
  data <- data[ok, , drop = FALSE]
  ids <- as.character(data[[id]])
  common <- Reduce(intersect, c(list(unique(ids)), lapply(glist, rownames)))
  if (length(common) < 3L) stop("Fewer than 3 donors overlap between `data` and the genetic input.", call. = FALSE)
  keep <- ids %in% common
  if (any(!keep)) message(sum(!keep), " observations dropped: donor not in genotype/GRM.")
  data <- data[keep, , drop = FALSE]
  ids <- ids[keep]
  donors <- rownames(glist[[1]])[rownames(glist[[1]]) %in% common]
  didx <- match(ids, donors)
  donors <- donors[sort(unique(didx))]                 # keep donors with data
  didx <- match(ids, donors)
  n <- length(donors)
  N <- nrow(data)

  kern <- c(list(one = matrix(1, n, n), I = diag(n)),
            lapply(glist, function(g) g[donors, donors, drop = FALSE]))

  yy <- data[[y]]
  yy <- as.numeric(scale(yy))                # expression is always standardized

  ## ---- components -----------------------------------------------------------
  ones <- matrix(1, N, 1L)
  comps <- list(list(name = "Intercept", group = "intercept", kernel = "one", U = ones,
                     fid = NA_integer_, grm = NA_character_, cvar = NA_character_))
  add <- function(name, group, kernel, U, fid = NA_integer_, grm = NA_character_, cvar = NA_character_) {
    comps[[length(comps) + 1L]] <<- list(name = name, group = group, kernel = kernel,
                                         U = U, fid = fid, grm = grm, cvar = cvar)
  }
  for (g in gnames) add(g, "G", g, ones, grm = g)
  if (ind_effect) add("I", "I", "I", ones)

  ctx <- .design_terms(data, context, cat_mode, scale_x)
  for (tm in ctx$terms) {
    if (context_main) add(tm$name, "context", "one", tm$U, tm$fid)
    if (gxc) for (g in gnames) add(paste0(g, ":", tm$name), "GxC", g, tm$U, tm$fid, g, tm$var)
  }
  cov <- .design_terms(data, covariates, "pooled", scale_x, fid_start = ctx$fid)
  for (tm in cov$terms) add(tm$name, "covariate", "one", tm$U, tm$fid)

  nms <- vapply(comps, function(cp) cp$name, "")
  if (anyDuplicated(nms)) stop("Duplicated component names: ", paste(nms[duplicated(nms)], collapse = ", "), call. = FALSE)
  p <- length(comps)
  if (N * (N - 1) / 2 < 5 * p) stop("Too few observations for this many components.", call. = FALSE)

  ## ---- jackknife blocks -----------------------------------------------------
  bid <- NULL
  if (jackknife) {
    nb <- if (is.null(n_blocks)) n else min(as.integer(n_blocks), n)
    if (nb < 3L) stop("`n_blocks` must be at least 3.", call. = FALSE)
    if (!is.null(seed)) set.seed(seed)
    bid <- if (nb == n) seq_len(n) else sample(rep_len(seq_len(nb), n))
  }

  ## ---- fit --------------------------------------------------------------------
  st <- .he_stats(comps, kern, didx, yy, bid)
  beta <- .solve_he(st$XtX, st$Xty)

  mult <- vapply(comps, function(cp) {
    if (cp$group == "intercept") return(NA_real_)
    mean(diag(kern[[cp$kernel]])[didx] * rowSums(cp$U^2))
  }, 0)
  group <- vapply(comps, function(cp) cp$group, "")
  grm_of <- vapply(comps, function(cp) cp$grm, "")
  ve <- beta * mult

  ## summary terms: linear maps of the per-component variance explained
  terms <- list(G = which(group == "G"))
  if (length(gnames) > 1L) for (g in gnames) terms[[g]] <- which(group == "G" & grm_of == g)
  if (any(group == "GxC")) {
    terms$GxC <- which(group == "GxC")
    if (length(gnames) > 1L) for (g in gnames) terms[[paste0(g, ":C")]] <- which(group == "GxC" & grm_of == g)
    terms$G_total <- which(group %in% c("G", "GxC"))
    cv <- vapply(comps, function(cp) cp$cvar, "")
    if (length(context) > 1L) for (v in context) terms[[paste0("GxC:", v)]] <- which(group == "GxC" & cv == v)
  }
  if (any(group == "I")) terms$I <- which(group == "I")
  if (any(group == "context")) terms$context <- which(group == "context")
  if (any(group == "covariate")) terms$covariate <- which(group == "covariate")
  terms$total_explained <- which(group != "intercept")
  A <- do.call(rbind, lapply(terms, function(ix) { a <- numeric(p); a[ix] <- 1; a }))
  rownames(A) <- names(terms)
  summ <- data.frame(term = names(terms),
                     variance = as.numeric(A %*% replace(ve, is.na(ve), 0)),
                     row.names = NULL, stringsAsFactors = FALSE)

  ## population-level heritability: variance / (total explained variance), i.e. with the
  ## cell-level residual removed from the denominator
  summ$h2 <- summ$variance / sum(replace(ve, is.na(ve), 0))
  summ$h2[summ$term == "total_explained"] <- NA_real_

  coefs <- data.frame(component = nms, group = group, grm = grm_of,
                      estimate = beta, variance = ve,
                      row.names = NULL, stringsAsFactors = FALSE)

  jk <- NULL
  if (jackknife) {
    nb <- max(bid)
    Bm <- matrix(NA_real_, nb, p, dimnames = list(NULL, nms))
    for (b in seq_len(nb)) {
      Bm[b, ] <- .solve_he(st$XtXb[b, , ], st$Xtyb[b, ], warn = FALSE)
    }
    good <- stats::complete.cases(Bm)
    if (sum(good) < nb) warning(nb - sum(good), " jackknife blocks failed and were dropped.", call. = FALSE)
    Bg <- Bm[good, , drop = FALSE]
    ng <- nrow(Bg)
    jse <- function(M) sqrt((ng - 1) / ng * colSums(sweep(M, 2, colMeans(M))^2))
    coefs$se <- jse(Bg)
    coefs$se_variance <- coefs$se * abs(mult)
    coefs$z <- coefs$estimate / coefs$se
    coefs$p <- 2 * stats::pnorm(-abs(coefs$z))
    V <- sweep(Bg, 2, mult, "*")
    V[is.na(V)] <- 0
    Sg <- V %*% t(A)
    summ$se <- jse(Sg)
    summ$se_h2 <- jse(Sg / rowSums(V))
    summ$se_h2[summ$term == "total_explained"] <- NA_real_
    jk <- list(estimates = Bm, n_blocks = nb, block = bid)
  }

  structure(list(
    coefficients = coefs, summary = summ,
    n_obs = N, n_donors = n, n_grm = length(gnames), grm_names = gnames,
    context = context, covariates = covariates, cat_mode = cat_mode,
    XtX = `dimnames<-`(st$XtX, list(nms, nms)), Xty = stats::setNames(st$Xty, nms),
    jackknife = jk, call = call
  ), class = "schint")
}
