#' Example data for scHINT
#'
#' Simulated genotypes and expression.
#' Expression of three genes with known variance architecture is
#' simulated (see `truth`).
#'
#' @format A list with
#' \describe{
#'   \item{geno}{400 x 250 integer dosage matrix, donor IDs `D001`... as row names.}
#'   \item{snp}{SNP IDs and (nominal) positions.}
#'   \item{cell}{12,000 cells of one cell type: `donor`, continuous `state`,
#'     categorical `subtype`, cell covariates `pc1`, `pc2`, donor covariates
#'     `age`, `sex`, `batch`, and expression `GENE_A`, `GENE_B`, `GENE_C`.}
#'   \item{sample}{2,000 donor x cell-type (5 types) pseudobulk values with the same genes.}
#'   \item{truth}{Variances used in the simulation (`g`, `gxc`, `i`, `cov`, `noise`).}
#' }
#' @source `data-raw/make_example.R`
"schint_example"

#' Example perturb-seq data for scHINT
#'
#' Simulated cells (`simulate_pert()`): `perturb` (60 perturbations in three groups
#' plus `"NT"` controls), `perturb_group`, cell type `celltype`, continuous `state`,
#' covariates `pc1`, `pc2` and expression of four genes `GENE_A`-`GENE_D`.
#' GENE_A has perturbation x state variance, GENE_B perturbation x cell type, GENE_C
#' perturbation only and GENE_D none; group3 perturbations have no effect.
#'
#' @format A data frame with one row per cell.
#' @source `data-raw/make_pert_example.R`
"schint_pert_example"
