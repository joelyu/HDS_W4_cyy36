# =============================================================================
# 05-ispy2.R
# Test: does cytotoxic effector score predict immunotherapy response (pCR)
# in the I-SPY2 pembrolizumab arm?
#
# Data: GSE194040 (I-SPY2 TRIAL, 988 breast cancer patients)
#
# Expects: scores_ispy2 from 01-immune-scoring.R
#          proc_dir defined by parent qmd
# =============================================================================

library(dplyr)
library(GEOquery)

if (!exists("proc_dir")) proc_dir <- file.path("data", "processed")

# --- 1. Load pre-computed immune scores (from 01-immune-scoring.R) ------------
if (!exists("scores_ispy2")) {
  scores_ispy2 <- read.csv(file.path(proc_dir, "ispy2_immune_scores.csv"))
}
message(sprintf("=== I-SPY2 immune scores: %d patients ===", nrow(scores_ispy2)))

# --- 2. Get phenotype data from GEO ------------------------------------------
gse <- getGEO("GSE194040", destdir = "data/raw", GSEMatrix = TRUE)
pd1 <- pData(gse[[1]])
pd2 <- pData(gse[[2]])
pd_all <- bind_rows(pd1, pd2) %>%
  mutate(
    patient_id = `patient id:ch1`,
    arm = `arm:ch1`,
    pcr = as.integer(`pcr:ch1`),
    hr = as.integer(`hr:ch1`),
    her2 = as.integer(`her2:ch1`)
  ) %>%
  select(patient_id, arm, pcr, hr, her2)

# --- 3. Merge and filter ------------------------------------------------------
merged <- inner_join(scores_ispy2, pd_all, by = "patient_id")
pembro <- merged %>% filter(grepl("Pembrolizumab", arm, ignore.case = TRUE))

message("\nArm label distribution (sanity check):")
print(sort(table(merged$arm), decreasing = TRUE))
message(sprintf("\n=== Pembrolizumab arm: %d patients ===", nrow(pembro)))
message(sprintf("pCR rate: %d/%d = %.1f%%", sum(pembro$pcr), nrow(pembro),
                100 * mean(pembro$pcr)))

# ============================================================
# 5. Logistic regression: individual cell types vs pCR
# ============================================================
message("\n=== Logistic regression: Immune scores predicting pCR ===")

test_cols <- c("Cytotoxic_cells", "CD8_T_cells", "Exhausted_CD8",
               "Mature_NK", "B_cells", "T_cells", "Macrophages",
               "Neutrophils", "Treg", "DC", "Th1_cells")

results <- lapply(test_cols, function(ct) {
  pembro$score_z <- scale(pembro[[ct]])[, 1]
  fit <- glm(pcr ~ score_z, data = pembro, family = binomial)
  s <- summary(fit)
  or <- exp(coef(fit)[2])
  ci <- exp(confint.default(fit)[2, ])
  data.frame(cell_type = ct, OR = or, OR_lower = ci[1], OR_upper = ci[2],
             p_value = s$coefficients[2, 4], stringsAsFactors = FALSE)
}) %>% bind_rows() %>% arrange(p_value)

message("\nUnadjusted (per SD):")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  sig <- ifelse(r$p_value < 0.05, "***", ifelse(r$p_value < 0.1, "*", ""))
  message(sprintf("  %-20s OR=%.2f [%.2f-%.2f] p=%.4f %s",
                  r$cell_type, r$OR, r$OR_lower, r$OR_upper, r$p_value, sig))
}

# ============================================================
# 6. Adjusted for trial stratification variable (HR status)
# ============================================================
# Note: HER2 not included — pembro and control arms are HER2-negative
# by I-SPY2 trial design (all her2=0), so it cannot vary as a covariate.
message("\n=== HR-adjusted ===")
results_adj <- lapply(test_cols, function(ct) {
  pembro$score_z <- scale(pembro[[ct]])[, 1]
  fit <- glm(pcr ~ score_z + factor(hr), data = pembro, family = binomial)
  s <- summary(fit)
  or <- exp(coef(fit)[2])
  ci <- exp(confint.default(fit)[2, ])
  data.frame(cell_type = ct, OR = or, OR_lower = ci[1], OR_upper = ci[2],
             p_value = s$coefficients[2, 4], stringsAsFactors = FALSE)
}) %>% bind_rows() %>% arrange(p_value)

for (i in seq_len(nrow(results_adj))) {
  r <- results_adj[i, ]
  sig <- ifelse(r$p_value < 0.05, "***", ifelse(r$p_value < 0.1, "*", ""))
  message(sprintf("  %-20s OR=%.2f [%.2f-%.2f] p=%.4f %s",
                  r$cell_type, r$OR, r$OR_lower, r$OR_upper, r$p_value, sig))
}

# ============================================================
# 7. Top-5 composite score
# ============================================================
message("\n=== Top-5 effector score ===")
top5 <- c("Cytotoxic_cells", "Exhausted_CD8", "CD8_T_cells", "Mature_NK", "B_cells")
pembro$top5_mean <- rowMeans(pembro[, top5])  # unscaled composite; scale() applied once in model

fit5 <- glm(pcr ~ scale(top5_mean), data = pembro, family = binomial)
s5 <- summary(fit5)
or5 <- exp(coef(fit5)[2])
ci5 <- exp(confint.default(fit5)[2, ])
message(sprintf("Unadjusted: OR=%.2f [%.2f-%.2f] p=%.4f",
                or5, ci5[1], ci5[2], s5$coefficients[2, 4]))

fit5a <- glm(pcr ~ scale(top5_mean) + factor(hr), data = pembro, family = binomial)
s5a <- summary(fit5a)
or5a <- exp(coef(fit5a)[2])
ci5a <- exp(confint.default(fit5a)[2, ])
message(sprintf("HR-adjusted: OR=%.2f [%.2f-%.2f] p=%.4f",
                or5a, ci5a[1], ci5a[2], s5a$coefficients[2, 4]))

# ============================================================
# 8. Control comparison: Paclitaxel-only arm
# ============================================================
message("\n=== Paclitaxel-only (control) arm ===")
control <- merged %>% filter(arm == "Paclitaxel" | arm == "Paclitaxel only")
message(sprintf("Control: %d patients, pCR rate: %.1f%%",
                nrow(control), 100 * mean(control$pcr)))

control$cyto_z <- scale(control$Cytotoxic_cells)[, 1]
fc <- glm(pcr ~ cyto_z + factor(hr), data = control, family = binomial)
sc <- summary(fc)
orc <- exp(coef(fc)[2])
cic <- exp(confint.default(fc)[2, ])
message(sprintf("Cytotoxic → pCR (control): OR=%.2f [%.2f-%.2f] p=%.4f",
                orc, cic[1], cic[2], sc$coefficients[2, 4]))

control$top5_mean <- rowMeans(control[, top5])  # unscaled composite; scale() applied once in model
fc5 <- glm(pcr ~ scale(top5_mean) + factor(hr), data = control, family = binomial)
sc5 <- summary(fc5)
orc5 <- exp(coef(fc5)[2])
cic5 <- exp(confint.default(fc5)[2, ])
message(sprintf("Top-5 → pCR (control):     OR=%.2f [%.2f-%.2f] p=%.4f",
                orc5, cic5[1], cic5[2], sc5$coefficients[2, 4]))

# ============================================================
# 9. Interaction: immune × pembrolizumab
# ============================================================
message("\n=== Interaction: immune score x pembrolizumab ===")
pc <- merged %>%
  filter(grepl("Pembrolizumab", arm, ignore.case = TRUE) |
           arm == "Paclitaxel" | arm == "Paclitaxel only") %>%
  mutate(is_pembro = ifelse(grepl("Pembrolizumab", arm), 1, 0))

# Note: scale() here uses the pooled (pembro + control) mean/SD — intentional for
# the interaction model so that the predictor is on a common scale across arms,
# making the interaction coefficient interpretable. Per-arm ORs in sections 5 & 8
# use arm-specific scaling so those ORs are not directly comparable to these.
pc$cyto_z <- scale(pc$Cytotoxic_cells)[, 1]
fi <- glm(pcr ~ cyto_z * factor(is_pembro) + factor(hr), data = pc, family = binomial)
si <- summary(fi)
message("Cytotoxic x Pembrolizumab interaction:")
print(round(si$coefficients, 4))

pc$top5_z <- scale(rowMeans(pc[, top5]))[, 1]  # single scale() on unscaled composite
fi5 <- glm(pcr ~ top5_z * factor(is_pembro) + factor(hr), data = pc, family = binomial)
si5 <- summary(fi5)
message("\nTop-5 x Pembrolizumab interaction:")
print(round(si5$coefficients, 4))

# ============================================================
# 10. Save
# ============================================================
pembro_out <- pembro %>%
  select(patient_id, arm, pcr, hr, her2, top5_mean,
         all_of(c("Cytotoxic_cells", "CD8_T_cells", "Exhausted_CD8",
                  "Mature_NK", "B_cells", "T_cells", "Macrophages",
                  "Neutrophils", "Treg", "DC", "Th1_cells")))
write.csv(pembro_out, file.path(proc_dir, "ispy2_pembro_scores.csv"), row.names = FALSE)
write.csv(results,     file.path(proc_dir, "ispy2_pembro_logistic_unadj.csv"), row.names = FALSE)
write.csv(results_adj, file.path(proc_dir, "ispy2_pembro_logistic_adj.csv"),   row.names = FALSE)
message("\n=== Done ===")
