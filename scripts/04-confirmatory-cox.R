# =============================================================================
# 04-confirmatory-cox.R
# Confirmatory Cox models: cytotoxic score as pre-specified predictor
#
# Pre-specification rationale: cytotoxic score selected in 03-pca-cytotoxic.R
#   - Rank #1 PC1 loading in both METABRIC and TCGA
#   - B cells (rank #5) adds no independent signal (LR test p ≈ 0.22)
#
# Models:
#   METABRIC: progressive adjustment f0→f1→f2→f3→f4
#     f0: unadjusted
#     f1: + PAM50
#     f2: + PAM50 + age
#     f3: + PAM50 + age + stage          ← PRIMARY (matches TCGA)
#     f4: + PAM50 + age + stage + grade  ← METABRIC sensitivity (grade available)
#   METABRIC NPI sensitivity: PAM50 + age + NPI (full n=1,756, no stage filter)
#
#   TCGA: progressive adjustment g0→g1→g2→g3
#     g0: unadjusted
#     g1: + PAM50
#     g2: + PAM50 + age
#     g3: + PAM50 + age + stage          ← PRIMARY
#     Note: grade not available in TCGA via cBioPortal
#
# Scaling convention: cyto_z and age_z are re-scaled on each model's own
#   analytical sample. HRs are "per 1 SD in that model's sample." This is
#   standard epidemiological practice for progressive adjustment tables.
#
# Assumption tests (Raquel's feedback):
#   1. Schoenfeld residuals — proportional hazards (HR constant over time?)
#   2. Martingale residuals — linearity of continuous predictors
#
# Expects: surv_metabric, surv_tcga from 02-explore.R
# =============================================================================

library(survival)
library(dplyr)

if (!exists("proc_dir")) proc_dir <- file.path("data", "processed")

if (!exists("surv_metabric")) {
  source("scripts/_immune_markers.R")
  source("scripts/01-immune-scoring.R")
  source("scripts/02-explore.R")
}

# =============================================================================
# 1. Load extended clinical confounders (fetched by 00-fetch-data.R)
# =============================================================================
m_ext <- read.csv(file.path(proc_dir, "metabric_clinical_extended.csv"))
t_ext <- read.csv(file.path(proc_dir, "tcga_clinical_extended.csv"))

# =============================================================================
# 2. Merge extended clinical
# =============================================================================
surv_m <- surv_metabric %>% left_join(m_ext, by = "patient_id")

surv_t <- surv_tcga %>%
  mutate(.pid12 = substr(patient_id, 1, 12)) %>%
  left_join(t_ext %>% mutate(.pid12 = substr(patient_id, 1, 12)) %>%
              select(-.data$patient_id), by = ".pid12") %>%
  select(-.pid12)

message(sprintf("\nMETABRIC: stage %d/%d | grade %d/%d | NPI %d/%d",
                sum(!is.na(surv_m$tumor_stage)), nrow(surv_m),
                sum(!is.na(surv_m$grade)),       nrow(surv_m),
                sum(!is.na(surv_m$npi)),          nrow(surv_m)))
message(sprintf("TCGA:     stage %d/%d",
                sum(!is.na(surv_t$tumor_stage_num)), nrow(surv_t)))

# =============================================================================
# 3. METABRIC: progressive adjustment
# =============================================================================
message("\n=== METABRIC: Cytotoxic score — progressive adjustment ===")

# f0: unadjusted
surv_m$cyto_z <- scale(surv_m$Cytotoxic_cells)[, 1]
f0 <- coxph(Surv(os_months, os_event) ~ cyto_z, data = surv_m)
s0 <- summary(f0)
message(sprintf("  f0 Unadjusted:                    HR=%.3f [%.3f-%.3f] p=%.2e  C=%.3f  n=%d",
                s0$conf.int[1,1], s0$conf.int[1,3], s0$conf.int[1,4],
                s0$coefficients[1,5], concordance(f0)$concordance, s0$n))

# f1: + PAM50
f1 <- coxph(Surv(os_months, os_event) ~ cyto_z + pam50_subtype, data = surv_m)
s1 <- summary(f1)
message(sprintf("  f1 + PAM50:                       HR=%.3f [%.3f-%.3f] p=%.2e  C=%.3f  n=%d",
                s1$conf.int[1,1], s1$conf.int[1,3], s1$conf.int[1,4],
                s1$coefficients[1,5], concordance(f1)$concordance, s1$n))

# f2: + PAM50 + age
surv_m_age       <- surv_m %>% filter(!is.na(age_at_diagnosis))
surv_m_age$cyto_z <- scale(surv_m_age$Cytotoxic_cells)[, 1]
surv_m_age$age_z  <- scale(surv_m_age$age_at_diagnosis)[, 1]
f2 <- coxph(Surv(os_months, os_event) ~ cyto_z + pam50_subtype + age_z,
            data = surv_m_age)
s2 <- summary(f2)
message(sprintf("  f2 + PAM50 + Age:                 HR=%.3f [%.3f-%.3f] p=%.2e  n=%d",
                s2$conf.int[1,1], s2$conf.int[1,3], s2$conf.int[1,4],
                s2$coefficients[1,5], s2$n))

# f3: + PAM50 + age + stage  ← PRIMARY MODEL
surv_m_stg        <- surv_m_age %>% filter(!is.na(tumor_stage) & tumor_stage > 0)
surv_m_stg$cyto_z <- scale(surv_m_stg$Cytotoxic_cells)[, 1]
surv_m_stg$age_z  <- scale(surv_m_stg$age_at_diagnosis)[, 1]
f3 <- coxph(Surv(os_months, os_event) ~ cyto_z + pam50_subtype + age_z +
              factor(tumor_stage), data = surv_m_stg)
s3 <- summary(f3)
message(sprintf("  f3 + PAM50 + Age + Stage [PRIMARY]:HR=%.3f [%.3f-%.3f] p=%.2e  n=%d",
                s3$conf.int[1,1], s3$conf.int[1,3], s3$conf.int[1,4],
                s3$coefficients[1,5], s3$n))

# f4: + PAM50 + age + stage + grade  (METABRIC-only sensitivity)
surv_m_full        <- surv_m_stg %>% filter(!is.na(grade))
surv_m_full$cyto_z <- scale(surv_m_full$Cytotoxic_cells)[, 1]
surv_m_full$age_z  <- scale(surv_m_full$age_at_diagnosis)[, 1]
f4 <- coxph(Surv(os_months, os_event) ~ cyto_z + pam50_subtype + age_z +
              factor(tumor_stage) + factor(grade), data = surv_m_full)
s4 <- summary(f4)
message(sprintf("  f4 + PAM50 + Age + Stage + Grade:  HR=%.3f [%.3f-%.3f] p=%.2e  n=%d",
                s4$conf.int[1,1], s4$conf.int[1,3], s4$conf.int[1,4],
                s4$coefficients[1,5], s4$n))

# NPI sensitivity: PAM50 + age + NPI on full cohort (no stage filter, n≈1,756)
# Scale on the NPI-complete subset — same per-model scaling rule as f2/f3/f4.
surv_m_npi        <- surv_m %>% filter(!is.na(npi))
surv_m_npi$cyto_z <- scale(surv_m_npi$Cytotoxic_cells)[, 1]
surv_m_npi$age_z  <- scale(surv_m_npi$age_at_diagnosis)[, 1]
surv_m_npi$npi_z  <- scale(surv_m_npi$npi)[, 1]
f_npi <- coxph(Surv(os_months, os_event) ~ cyto_z + pam50_subtype + age_z + npi_z,
               data = surv_m_npi)
s_npi <- summary(f_npi)
message(sprintf("  PAM50 + Age + NPI (no stage filter): HR=%.3f [%.3f-%.3f] p=%.2e  n=%d",
                s_npi$conf.int[1,1], s_npi$conf.int[1,3], s_npi$conf.int[1,4],
                s_npi$coefficients[1,5], s_npi$n))

# =============================================================================
# 4. TCGA: progressive adjustment
# =============================================================================
message("\n=== TCGA: Cytotoxic score — progressive adjustment ===")

surv_t$cyto_z <- scale(surv_t$Cytotoxic_cells)[, 1]

# g0: unadjusted
g0  <- coxph(Surv(os_months, os_event) ~ cyto_z, data = surv_t)
sg0 <- summary(g0)
message(sprintf("  g0 Unadjusted:                    HR=%.3f [%.3f-%.3f] p=%.2e  C=%.3f  n=%d",
                sg0$conf.int[1,1], sg0$conf.int[1,3], sg0$conf.int[1,4],
                sg0$coefficients[1,5], concordance(g0)$concordance, sg0$n))

# g1: + PAM50
g1  <- coxph(Surv(os_months, os_event) ~ cyto_z + pam50_subtype, data = surv_t)
sg1 <- summary(g1)
message(sprintf("  g1 + PAM50:                       HR=%.3f [%.3f-%.3f] p=%.2e  C=%.3f  n=%d",
                sg1$conf.int[1,1], sg1$conf.int[1,3], sg1$conf.int[1,4],
                sg1$coefficients[1,5], concordance(g1)$concordance, sg1$n))

# g2: + PAM50 + age
surv_t_age        <- surv_t %>% filter(!is.na(age_at_diagnosis))
surv_t_age$cyto_z <- scale(surv_t_age$Cytotoxic_cells)[, 1]
surv_t_age$age_z  <- scale(surv_t_age$age_at_diagnosis)[, 1]
g2  <- coxph(Surv(os_months, os_event) ~ cyto_z + pam50_subtype + age_z,
             data = surv_t_age)
sg2 <- summary(g2)
message(sprintf("  g2 + PAM50 + Age:                 HR=%.3f [%.3f-%.3f] p=%.2e  n=%d",
                sg2$conf.int[1,1], sg2$conf.int[1,3], sg2$conf.int[1,4],
                sg2$coefficients[1,5], sg2$n))

# g3: + PAM50 + age + stage  ← PRIMARY MODEL (grade unavailable in TCGA)
surv_t_stg        <- surv_t_age %>% filter(!is.na(tumor_stage_num))
surv_t_stg$cyto_z <- scale(surv_t_stg$Cytotoxic_cells)[, 1]
surv_t_stg$age_z  <- scale(surv_t_stg$age_at_diagnosis)[, 1]
g3  <- coxph(Surv(os_months, os_event) ~ cyto_z + pam50_subtype + age_z +
               factor(tumor_stage_num), data = surv_t_stg)
sg3 <- summary(g3)
message(sprintf("  g3 + PAM50 + Age + Stage [PRIMARY]:HR=%.3f [%.3f-%.3f] p=%.2e  n=%d",
                sg3$conf.int[1,1], sg3$conf.int[1,3], sg3$conf.int[1,4],
                sg3$coefficients[1,5], sg3$n))

# =============================================================================
# 5. Assumption tests
# =============================================================================

# --- 5a. Schoenfeld residuals (proportional hazards) -------------------------
# Tests whether the HR is constant over time.
# Applied to f3 (METABRIC primary) and g3 (TCGA primary).
# If p < 0.05 for cytotoxic: the effect is time-varying — report as finding.
message("\n=== Assumption test 1: Proportional hazards (Schoenfeld residuals) ===")

ph_m <- cox.zph(f3)
message("\nMETABRIC f3 (PAM50 + Age + Stage):")
print(ph_m)

ph_t <- cox.zph(g3)
message("\nTCGA g3 (PAM50 + Age + Stage):")
print(ph_t)

# Also test unadjusted cytotoxic for comparison
ph_cyto_m <- cox.zph(f0)
message("\nMETABRIC f0 (cytotoxic unadjusted — for reference):")
print(ph_cyto_m)

# --- 5b. Martingale residuals (linearity) ------------------------------------
# Tests whether continuous predictors (cytotoxic score, age) have a linear
# relationship with log-hazard. Null model is fitted WITHOUT cytotoxic so the
# residuals carry its full unmodelled signal — plotting/correlating residuals
# against raw cytotoxic score reveals the functional shape it should take.
# A flat or monotone-linear shape (no curvature) means the linear term is
# adequate; visual inspection of the loess smoother is the primary check.
message("\n=== Assumption test 2: Linearity (Martingale residuals) ===")

# METABRIC null model (same analytical sample as f3)
null_m    <- coxph(Surv(os_months, os_event) ~ pam50_subtype + age_z + factor(tumor_stage),
                   data = surv_m_stg)
mart_resid_m <- residuals(null_m, type = "martingale")
cor_cyto_m <- cor(mart_resid_m, surv_m_stg$Cytotoxic_cells, method = "spearman")
cor_age_m  <- cor(mart_resid_m, surv_m_stg$age_at_diagnosis, method = "spearman")
message(sprintf("METABRIC Spearman(Martingale, Cytotoxic): r = %.3f", cor_cyto_m))
message(sprintf("METABRIC Spearman(Martingale, Age):       r = %.3f", cor_age_m))

# TCGA null model (same analytical sample as g3)
null_t    <- coxph(Surv(os_months, os_event) ~ pam50_subtype + age_z + factor(tumor_stage_num),
                   data = surv_t_stg)
mart_resid_t <- residuals(null_t, type = "martingale")
cor_cyto_t <- cor(mart_resid_t, surv_t_stg$Cytotoxic_cells, method = "spearman")
cor_age_t  <- cor(mart_resid_t, surv_t_stg$age_at_diagnosis, method = "spearman")
message(sprintf("TCGA Spearman(Martingale, Cytotoxic):     r = %.3f", cor_cyto_t))
message(sprintf("TCGA Spearman(Martingale, Age):           r = %.3f", cor_age_t))

# =============================================================================
# 6. Save outputs
# =============================================================================
# Progressive adjustment table
cox_metabric <- data.frame(
  model    = c("Unadjusted", "+PAM50", "+PAM50+Age",
               "+PAM50+Age+Stage [PRIMARY]", "+PAM50+Age+Stage+Grade",
               "+PAM50+Age+NPI (full n)"),
  cohort   = "METABRIC",
  HR       = c(s0$conf.int[1,1], s1$conf.int[1,1], s2$conf.int[1,1],
               s3$conf.int[1,1], s4$conf.int[1,1], s_npi$conf.int[1,1]),
  HR_lower = c(s0$conf.int[1,3], s1$conf.int[1,3], s2$conf.int[1,3],
               s3$conf.int[1,3], s4$conf.int[1,3], s_npi$conf.int[1,3]),
  HR_upper = c(s0$conf.int[1,4], s1$conf.int[1,4], s2$conf.int[1,4],
               s3$conf.int[1,4], s4$conf.int[1,4], s_npi$conf.int[1,4]),
  p_value  = c(s0$coefficients[1,5], s1$coefficients[1,5], s2$coefficients[1,5],
               s3$coefficients[1,5], s4$coefficients[1,5], s_npi$coefficients[1,5]),
  n        = c(s0$n, s1$n, s2$n, s3$n, s4$n, s_npi$n)
)

cox_tcga <- data.frame(
  model    = c("Unadjusted", "+PAM50", "+PAM50+Age", "+PAM50+Age+Stage [PRIMARY]"),
  cohort   = "TCGA",
  HR       = c(sg0$conf.int[1,1], sg1$conf.int[1,1], sg2$conf.int[1,1], sg3$conf.int[1,1]),
  HR_lower = c(sg0$conf.int[1,3], sg1$conf.int[1,3], sg2$conf.int[1,3], sg3$conf.int[1,3]),
  HR_upper = c(sg0$conf.int[1,4], sg1$conf.int[1,4], sg2$conf.int[1,4], sg3$conf.int[1,4]),
  p_value  = c(sg0$coefficients[1,5], sg1$coefficients[1,5], sg2$coefficients[1,5],
               sg3$coefficients[1,5]),
  n        = c(sg0$n, sg1$n, sg2$n, sg3$n)
)

cox_robustness <- bind_rows(cox_metabric, cox_tcga)
write.csv(cox_robustness, file.path(proc_dir, "cox_robustness_confounders.csv"),
          row.names = FALSE)

# Schoenfeld test results — both cohorts
ph_df_m      <- as.data.frame(ph_m$table)
ph_df_m$term <- rownames(ph_df_m)
ph_df_m$cohort <- "METABRIC"
write.csv(ph_df_m, file.path(proc_dir, "schoenfeld_test_metabric.csv"), row.names = FALSE)

ph_df_t      <- as.data.frame(ph_t$table)
ph_df_t$term <- rownames(ph_df_t)
ph_df_t$cohort <- "TCGA"
write.csv(ph_df_t, file.path(proc_dir, "schoenfeld_test_tcga.csv"), row.names = FALSE)

# Martingale test results — both cohorts
martingale_df <- data.frame(
  cohort     = c("METABRIC", "METABRIC", "TCGA", "TCGA"),
  predictor  = c("Cytotoxic", "Age", "Cytotoxic", "Age"),
  spearman_r = c(cor_cyto_m, cor_age_m, cor_cyto_t, cor_age_t),
  n          = c(nrow(surv_m_stg), nrow(surv_m_stg), nrow(surv_t_stg), nrow(surv_t_stg))
)
write.csv(martingale_df, file.path(proc_dir, "martingale_test.csv"), row.names = FALSE)

# Clean up temporary columns written directly onto surv_m / surv_t
# (surv_m$age_z and surv_m$npi_z are NOT written — only on subsets — so not cleaned here)
surv_m$cyto_z <- NULL
surv_t$cyto_z <- NULL

message("\n=== 04-confirmatory-cox.R complete ===")
message(sprintf("Outputs: %s/cox_robustness_confounders.csv", proc_dir))
message(sprintf("         %s/schoenfeld_test_metabric.csv", proc_dir))
message(sprintf("         %s/schoenfeld_test_tcga.csv", proc_dir))
message(sprintf("         %s/martingale_test.csv", proc_dir))
