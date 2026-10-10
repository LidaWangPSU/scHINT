#' Cell-level G×Cell-State heritability model
#'
#' Haseman-Elston regression on pairs of cells that partitions the variance of
#' one gene's expression into a main genetic component (`G`, one per GRM), genetic
#' effects that change along a cell-state/context variable (`G×context`), a
#' shared-donor component (`I`), context main effects and other covariates.
#'
#' For a pair of cells `a`, `b` from donors `i`, `j` the model regresses
#' `y_a * y_b` on `1`, `GRM[i,j]`, `1(i == j)`, `u_a * u_b` (context and
#' covariate main effects) and `GRM[i,j] * u_a * u_b` (genotype-by-context). The
#' estimate for each component is its variance contribution when `y` and the
#' covariates are standardized.
#'
#' @param data Data frame with one row per cell.
#' @param y Name of the (numeric) expression column.
#' @param id Name of the donor ID column. IDs must match the row names of `geno`/`grm`.
#' @param geno Donors x SNPs dosage matrix of the cis SNPs, with donor IDs as row
#'   names; a single GRM is built from all SNPs. Give either `geno` or `grm`.
#' @param grm A pre-computed donor x donor GRM with donor IDs as dimnames, or a
#'   (named) list of GRMs for a multi-GRM model (one main genetic and one
#'   genotype-by-context component per GRM).
#' @param context Character vector of column names used as the cell context that
#'   interacts with genetics. Numeric columns are continuous contexts (e.g. a
#'   pseudotime or PC); factor/character columns are categorical. `NULL` fits a
#'   model without genotype-by-context terms.
#' @param covariates Character vector of other covariates (main effects only).
#'   Variables also listed in `context` are not duplicated.
#' @param cat_mode For categorical variables in `context`: `"pooled"` fits one
#'   shared component (pairs from the same level) and `"per_level"` fits one per level.
#' @param ind_effect Include the shared-donor component `I` (default `TRUE`).
#' @param context_main Include main-effect components for each context variable (default `TRUE`).
#' @param jackknife If `TRUE`, delete-donor-block jackknife standard errors.
#' @param n_blocks Number of jackknife blocks over donors (default: one per donor).
#' @param scale_x Standardize continuous context and covariates (default `TRUE`).
#'   Expression `y` is always standardized to mean 0 and variance 1, so variance
#'   components are on the scale of standardized expression.
#' @param seed Seed for assigning donors to jackknife blocks.
#' @return An object of class `"schint"` with `coefficients`, `summary` (grouped
#'   variance components), the jackknife replicates and the sufficient statistics
#'   `XtX`, `Xty`.
#' @examples
#' data(schint_example)
#' d <- schint_example$cell
#' fit <- schint_cell(d, y = "GENE_A", id = "donor",
#'                    geno = schint_example$geno,
#'                    context = "state", covariates = c("pc1", "pc2", "age", "sex", "batch"))
#' fit
#' @export
schint_cell <- function(data, y, id, geno = NULL, grm = NULL,
                        context = NULL, covariates = NULL,
                        cat_mode = c("pooled", "per_level"),
                        ind_effect = TRUE, context_main = TRUE,
                        jackknife = FALSE, n_blocks = NULL,
                        scale_x = TRUE, seed = 1) {
  .schint_fit(data, y, id, context, covariates, geno, grm,
              match.arg(cat_mode), gxc = length(context) > 0L, ind_effect, context_main,
              scale_x, jackknife, n_blocks, seed, match.call())
}

#' Sample-level G×Cell-Type heritability model
#'
#' The sample-level (pseudobulk) counterpart of [schint_cell()]: each row is one
#' donor x cell-type expression value. Variance is partitioned into a genetic
#' component shared across cell types (`G`), cell-type-specific genetics
#' (`G×Cell-Type`, pairs from the same cell type), a donor component shared
#' across cell types (`I`), a cell-type component (`celltype`) and covariates.
#'
#' @inheritParams schint_cell
#' @param data Data frame with one row per donor x cell type.
#' @param celltype Name of the cell-type column.
#' @param gxc Fit cell-type-specific genetics (default `TRUE`); `FALSE` gives the
#'   shared-genetics-only model.
#' @param celltype_main Include the cell-type main-effect component (default `TRUE`).
#' @param cat_mode `"pooled"` (one `G×Cell-Type` component, default) or `"per_level"`
#'   (a separate genetic variance for each cell type).
#' @param covariates Other covariates (donor-level or observation-level).
#' @inherit schint_cell return
#' @examples
#' data(schint_example)
#' fit <- schint_sample(schint_example$sample, y = "GENE_A", id = "donor",
#'                      celltype = "celltype", geno = schint_example$geno,
#'                      covariates = c("age", "sex", "pc1", "pc2"))
#' fit
#' @export
schint_sample <- function(data, y, id, celltype, geno = NULL, grm = NULL,
                          covariates = NULL, gxc = TRUE, cat_mode = c("pooled", "per_level"),
                          ind_effect = TRUE, celltype_main = TRUE,
                          jackknife = FALSE, n_blocks = NULL,
                          scale_x = TRUE, seed = 1) {
  data <- as.data.frame(data)
  if (!celltype %in% names(data)) stop("`celltype` column not found in `data`.", call. = FALSE)
  data[[celltype]] <- factor(data[[celltype]])
  if (anyDuplicated(data[, c(id, celltype)])) {
    warning("Several rows per donor x cell type; sample-level data should have one. ",
            "Use schint_cell() for cell-level data.", call. = FALSE)
  }
  .schint_fit(data, y, id, context = celltype, covariates, geno, grm,
              match.arg(cat_mode), gxc = gxc, ind_effect, context_main = celltype_main,
              scale_x, jackknife, n_blocks, seed, match.call())
}

.compact_summary <- function(s, digits) {
  s <- s[, !vapply(s, function(v) all(is.na(v)), TRUE), drop = FALSE]
  num <- vapply(s, is.numeric, TRUE)
  s[num] <- lapply(s[num], function(v) signif(v, digits))
  s
}

#' @export
print.schint <- function(x, digits = 3, ...) {
  cat("scHINT model: ", x$n_obs, " observations, ", x$n_donors, " donors, ",
      x$n_grm, " GRM", if (x$n_grm > 1) "s", "\n\n", sep = "")
  cat("Variance components (variance on the standardized-expression scale):\n")
  print(.compact_summary(x$summary, digits), row.names = FALSE, na.print = "")
  if (!is.null(x$jackknife)) {
    cat("\nStandard errors: jackknife over ", x$jackknife$n_blocks, " donor blocks\n", sep = "")
  }
  invisible(x)
}

#' @export
summary.schint <- function(object, ...) object$coefficients

#' @export
coef.schint <- function(object, ...) {
  stats::setNames(object$coefficients$estimate, object$coefficients$component)
}
