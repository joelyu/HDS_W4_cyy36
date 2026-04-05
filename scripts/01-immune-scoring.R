# =============================================================================
# 01-immune-scoring.R
# Compute Danaher immune cell scores for METABRIC + TCGA-BRCA + I-SPY2
#
# Method: mean log2 expression of cell-type-specific marker genes per patient
# Reference: Danaher P et al. J Immunother Cancer 5, 18 (2017)
#
# Input:  data/processed/{cohort}_expression.csv  (from 00-fetch-data.R)
# Output: data/processed/{cohort}_immune_scores.csv
#         In memory: scores_both (METABRIC+TCGA), scores_all (all three)
#
# Expects: proc_dir defined by parent qmd
# =============================================================================

if (!exists("proc_dir")) proc_dir <- file.path("data", "processed")
source("scripts/_immune_markers.R")
library(dplyr)
library(tidyr)
library(tibble)

# --- Score function (reused from module 3 assignment) -----------------------------------------
compute_cell_score <- function(expr_matrix, genes) {
  available <- genes[genes %in% rownames(expr_matrix)]
  if (length(available) == 0) {
    return(rep(NA_real_, ncol(expr_matrix)))
  }
  if (length(available) == 1) {
    return(expr_matrix[available, ])
  }
  colMeans(expr_matrix[available, , drop = FALSE], na.rm = TRUE)
}

score_cohort <- function(expr_file, label) {
  expr_raw <- read.csv(expr_file, check.names = FALSE)
  gene_col <- colnames(expr_raw)[1]

  # Handle KIAA0125/FAM30A alias
  if ("KIAA0125" %in% expr_raw[[gene_col]] && !"FAM30A" %in% expr_raw[[gene_col]]) {
    expr_raw[[gene_col]][expr_raw[[gene_col]] == "KIAA0125"] <- "FAM30A"
  }

  expr_mat <- expr_raw %>%
    column_to_rownames(var = gene_col) %>%
    as.matrix()

  # Gene coverage report
  all_markers <- unique(unlist(immune_markers))
  found <- all_markers[all_markers %in% rownames(expr_mat)]
  missing <- all_markers[!all_markers %in% rownames(expr_mat)]

  message(sprintf("[%s] %d/%d marker genes found", label, length(found), length(all_markers)))
  if (length(missing) > 0) {
    message(sprintf("  Missing: %s", paste(missing, collapse = ", ")))
  }

  # Compute scores
  scores_list <- lapply(immune_markers, compute_cell_score, expr_matrix = expr_mat)
  scores_list[["CD4 T cells"]] <- scores_list[["T-cells"]] - scores_list[["CD8 T cells"]]

  scores_mat <- do.call(rbind, scores_list)
  colnames(scores_mat) <- colnames(expr_mat)

  scores_df <- scores_mat %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("patient_id") %>%
    mutate(cohort = label)

  # Clean column names
  colnames(scores_df) <- gsub("[ -]", "_", colnames(scores_df))

  out_file <- file.path(proc_dir, paste0(label, "_immune_scores.csv"))
  write.csv(scores_df, out_file, row.names = FALSE)
  message(sprintf(
    "[%s] Scores written: %d patients x %d cell types",
    label, nrow(scores_df), nrow(scores_mat)
  ))

  scores_df
}

# --- Score all three cohorts ---------------------------------------------------
scores_metabric <- score_cohort(file.path(proc_dir, "metabric_expression.csv"), "metabric")
scores_tcga <- score_cohort(file.path(proc_dir, "tcga_expression.csv"), "tcga")
scores_ispy2 <- score_cohort(file.path(proc_dir, "ispy2_expression.csv"), "ispy2")

# --- Combined long-format for easy comparison --------------------------------
scores_both <- bind_rows(scores_metabric, scores_tcga)
scores_all  <- bind_rows(scores_metabric, scores_tcga, scores_ispy2)
message(sprintf(
  "\nCombined: %d patients (%d METABRIC, %d TCGA, %d I-SPY2)",
  nrow(scores_all), nrow(scores_metabric), nrow(scores_tcga), nrow(scores_ispy2)
))
