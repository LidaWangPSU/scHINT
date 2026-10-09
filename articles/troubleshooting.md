# Troubleshooting

## Input and data problems

**`Fewer than 3 donors overlap between data and the genetic input`**
Donor IDs in `data[[id]]` do not match `rownames(geno)` / the GRM
dimnames. Compare `head(unique(data$donor))` with
`head(rownames(geno))`; IDs are compared as character strings (watch for
`FID:IID` vs `IID`, suffixes, or numeric IDs with leading zeros).

**`Provide exactly one of geno ... or grm`** Give genotypes **or** a
GRM, not both and not neither.

**`Columns not found in data`** `y`, `id`, `context`, `covariates` must
be column names.

**`... is a single level` / constant column** A categorical context or
covariate with one level, or a numeric one with zero variance (scaling
gives `NaN`), cannot be used; drop it.

**Observations dropped for missing values** Rows with `NA` in any used
column (and non-finite expression) are removed; the number is reported.
Impute covariates or drop them if too many cells are lost.

**`This many components ... too few observations`** Fewer pairs than
components. Reduce `context`/`covariates` or add data.

## Numerical warnings

**`X'X is singular; adding a small ridge`** Two components are
collinear. Common causes: the same variable listed as a covariate and as
a context (handled automatically), covariates that are exact functions
of each other (e.g. all levels of a batch plus donor-level constants), a
context that is constant within each donor, or a GRM (almost)
proportional to the identity. Remove redundant covariates.

**A jackknife block failed** A leave-block-out fit was singular and was
dropped (with a warning); with very few donors use fewer, larger blocks
(`n_blocks`).

## Unexpected results

**Negative variance estimates.** HE regression is unconstrained, so
small negative values are expected when the true component is near zero;
they are not an error. Report them as they are, aggregate over genes, or
truncate only for display.

**Estimates \> 1 or `h2_upper` \> 1.** `h2_upper` divides by a smaller
denominator and can exceed 1 when noisy; use `h2` and standard errors.

**Large standard errors.** Precision is governed by the number of donors
and cis SNPs. Many donors with few cells beat few donors with many
cells.

**Different numbers from another package.** Check the denominator
(`fraction` vs `h2`), the standardization (`scale_y`, `scale_x`), the
GRM definition (standardized genotypes, divided by the number of SNPs)
and whether the other method includes a donor component.

**Context main effect or covariate variance is large.** That is allowed
to be large (it is not genetic); it enters `total_explained` and the
`h2` denominator.

## Perturbation model

**Huge `P` for controls.** Non-targeting cells form one big
pseudo-perturbation. Use `control = "NT"`.

**`Skipping gene`.** The gene has missing/non-finite values or zero
variance.

**Fewer than 3 perturbations remain.** Lower `min_cells` or check the
`perturb` column.
