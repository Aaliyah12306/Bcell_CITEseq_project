# ============================================================
# 48-WEEK MANUSCRIPT SUMMARY
#
# Output:
# 1. W12/W24/W48 response/remission: mNRI + AO
# 2. Baseline predictors of W48 remission:
#    Firth univariable/multivariable OR, 95% CI, P
# 3. Early and sustained remission
# 4. Early remission predicting W48 remission
#
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(logistf)
library(tibble)


# ============================================================
# 1. Read efficacy and treatment-strategy data
# ============================================================

tx <- read_excel(
  "treatment_strategy.xlsx",
  sheet = "SAE"
)

stopifnot(nrow(dat) == nrow(tx))

cat(
  "\nWARNING:\n",
  "The following merge assumes that the two Excel files have\n",
  "exactly the same patient order.\n\n"
)


# ============================================================
# 2. Rename efficacy variables
# ============================================================




# ============================================================
# 3. Add treatment-strategy information
#
# IMPORTANT:
# row-order merge only
# ============================================================

dat <- bind_cols(
  dat,
  tx %>%
    select(
      drug_time,
      drug_event,
      stop_reason,
      strategy24_plot,
      strategy48_plot
    )
)


# ============================================================
# 4. Helper functions
# ============================================================

as_binary <- function(x) {
  
  z <- trimws(
    as.character(x)
  )
  
  z[
    z %in% c(
      "",
      "NA",
      "N/A",
      "NaN"
    )
  ] <- NA
  
  as.integer(z)
}


extract_week_number <- function(x) {
  as.numeric(
    stringr::str_extract(
      as.character(x),
      "\\d+\\.?\\d*"
    )
  )
}


format_p <- function(p) {
  
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "<0.001",
    TRUE ~ sprintf("%.3f", p)
  )
}


format_ci <- function(low, high) {
  
  sprintf(
    "%.2f–%.2f",
    low,
    high
  )
}


format_or <- function(or, low, high) {
  
  sprintf(
    "%.2f (%.2f–%.2f)",
    or,
    low,
    high
  )
}


prop_ci <- function(x, n) {
  
  if (n == 0) {
    
    return(
      tibble(
        Events = 0,
        N = 0,
        Percent = NA_real_,
        CI_low = NA_real_,
        CI_high = NA_real_
      )
    )
  }
  
  bt <- binom.test(
    x = x,
    n = n
  )
  
  tibble(
    Events = x,
    N = n,
    Percent = 100 * x / n,
    CI_low = 100 * bt$conf.int[1],
    CI_high = 100 * bt$conf.int[2]
  )
}


# ============================================================
# 7. Define treatment failures for mNRI
# ============================================================

dat <- dat %>%
  mutate(
    # 使用 extract_week_number 函数处理 drug_week 列
    drug_week_numeric = extract_week_number(drug_time),  # 数据框列 drug_week
    
    # 现在可以进行数值比较
    failure_by12 = drug_event == 1 & drug_week_numeric <= 12,
    failure_by24 = drug_event == 1 & drug_week_numeric <= 24,
    failure_by48 = drug_event == 1 & drug_week_numeric <= 48,
    advanced_rescue_by48 = strategy24_plot == 2
  )


# ============================================================
# 8. Generate mNRI endpoints
# ============================================================

dat <- dat %>%
  mutate(
    
    # ========================================================
    # Week 12
    # ========================================================
    
    w12_response_mnri = case_when(
      failure_by12 ~ 0L,
      is.na(w12_response) ~ 0L,
      TRUE ~ w12_response
    ),
    
    w12_remission_mnri = case_when(
      failure_by12 ~ 0L,
      is.na(w12_remission) ~ 0L,
      TRUE ~ w12_remission
    ),
    
    
    # ========================================================
    # Week 24
    # ========================================================
    
    w24_response_mnri = case_when(
      failure_by24 ~ 0L,
      is.na(w24_response) ~ 0L,
      TRUE ~ w24_response
    ),
    
    w24_remission_mnri = case_when(
      failure_by24 ~ 0L,
      is.na(w24_remission) ~ 0L,
      TRUE ~ w24_remission
    ),
    
    
    # ========================================================
    # Week 48
    #
    # Failure:
    # 1. permanent discontinuation through W48
    # 2. advanced therapy rescue during W24–48
    # 3. missing W48 outcome
    #
    # UPA dose escalation alone is NOT failure.
    # ========================================================
    
    w48_response_mnri = case_when(
      failure_by48 ~ 0L,
      advanced_rescue_by48 ~ 0L,
      is.na(w48_response) ~ 0L,
      TRUE ~ w48_response
    ),
    
    w48_remission_mnri = case_when(
      failure_by48 ~ 0L,
      advanced_rescue_by48 ~ 0L,
      is.na(w48_remission) ~ 0L,
      TRUE ~ w48_remission
    )
  )


# ============================================================
# 9. Function for mNRI / AO summary
# ============================================================

endpoint_summary <- function(
    raw_var,
    mnri_var,
    visit,
    endpoint
) {
  
  # ----------------------------
  # mNRI
  # ----------------------------
  
  x_mnri <-
    sum(
      dat[[mnri_var]] == 1,
      na.rm = TRUE
    )
  
  mnri <-
    prop_ci(
      x_mnri,
      nrow(dat)
    ) %>%
    mutate(
      Visit = visit,
      Endpoint = endpoint,
      Analysis = "mNRI"
    )
  
  
  # ----------------------------
  # As observed
  # ----------------------------
  
  x_raw <- dat[[raw_var]]
  
  n_ao <-
    sum(
      !is.na(x_raw)
    )
  
  events_ao <-
    sum(
      x_raw == 1,
      na.rm = TRUE
    )
  
  ao <-
    prop_ci(
      events_ao,
      n_ao
    ) %>%
    mutate(
      Visit = visit,
      Endpoint = endpoint,
      Analysis = "AO"
    )
  
  
  bind_rows(
    mnri,
    ao
  )
}


# ============================================================
# 10. SUMMARY 1
# W12 / W24 / W48 clinical response and remission
# ============================================================

summary_efficacy <- bind_rows(
  
  endpoint_summary(
    "w12_response",
    "w12_response_mnri",
    "Week 12",
    "Clinical response"
  ),
  
  endpoint_summary(
    "w12_remission",
    "w12_remission_mnri",
    "Week 12",
    "Clinical remission"
  ),
  
  endpoint_summary(
    "w24_response",
    "w24_response_mnri",
    "Week 24",
    "Clinical response"
  ),
  
  endpoint_summary(
    "w24_remission",
    "w24_remission_mnri",
    "Week 24",
    "Clinical remission"
  ),
  
  endpoint_summary(
    "w48_response",
    "w48_response_mnri",
    "Week 48",
    "Clinical response"
  ),
  
  endpoint_summary(
    "w48_remission",
    "w48_remission_mnri",
    "Week 48",
    "Clinical remission"
  )
  
) %>%
  
  mutate(
    
    `n/N (%)` =
      sprintf(
        "%d/%d (%.1f%%)",
        Events,
        N,
        Percent
      ),
    
    `95% CI` =
      sprintf(
        "%.1f–%.1f%%",
        CI_low,
        CI_high
      )
  ) %>%
  
  select(
    Visit,
    Endpoint,
    Analysis,
    `n/N (%)`,
    `95% CI`,
    Events,
    N,
    Percent,
    CI_low,
    CI_high
  )


# ============================================================
# 11. Early and sustained remission
# ============================================================

dat <- dat %>%
  mutate(
    
    early_remission =
      as.integer(
        w12_remission_mnri == 1
      ),
    
    # Strict sustained clinical remission:
    # remission at W12 + W24 + W48
    # with W48 treatment-failure rules already incorporated
    sustained_remission =
      as.integer(
        w12_remission_mnri == 1 &
          w24_remission_mnri == 1 &
          w48_remission_mnri == 1
      ),
    
    early_remission_factor =
      factor(
        early_remission,
        levels = c(0, 1),
        labels = c(
          "No early remission",
          "Early remission"
        )
      )
  )


# ============================================================
# 12. SUMMARY 2
# Early and sustained remission proportions
# ============================================================

early_prop <- prop_ci(
  sum(dat$early_remission == 1),
  nrow(dat)
) %>%
  mutate(
    Outcome =
      "Early clinical remission"
  )


sustained_prop <- prop_ci(
  sum(dat$sustained_remission == 1),
  nrow(dat)
) %>%
  mutate(
    Outcome =
      "Sustained clinical remission"
  )


sustained_among_early <- prop_ci(
  sum(dat$sustained_remission == 1),
  sum(dat$early_remission == 1)
) %>%
  mutate(
    Outcome =
      "Sustained remission among early remitters"
  )


summary_durability <- bind_rows(
  early_prop,
  sustained_prop,
  sustained_among_early
) %>%
  
  mutate(
    
    `n/N (%)` =
      sprintf(
        "%d/%d (%.1f%%)",
        Events,
        N,
        Percent
      ),
    
    `95% CI` =
      sprintf(
        "%.1f–%.1f%%",
        CI_low,
        CI_high
      )
  ) %>%
  
  select(
    Outcome,
    `n/N (%)`,
    `95% CI`,
    Events,
    N,
    Percent,
    CI_low,
    CI_high
  )


# ============================================================
# 13. Firth regression helper
# ============================================================

extract_firth <- function(
    fit,
    model_name
) {
  
  tibble(
    Term =
      names(
        fit$coefficients
      ),
    
    OR =
      exp(
        as.numeric(
          fit$coefficients
        )
      ),
    
    CI_low =
      exp(
        as.numeric(
          fit$ci.lower
        )
      ),
    
    CI_high =
      exp(
        as.numeric(
          fit$ci.upper
        )
      ),
    
    P_value =
      as.numeric(
        fit$prob
      ),
    
    Model =
      model_name
  ) %>%
    filter(
      Term != "(Intercept)"
    )
}


fit_firth <- function(
    outcome,
    predictors,
    data = dat
) {
  
  f <- as.formula(
    paste(
      outcome,
      "~",
      paste(
        predictors,
        collapse = " + "
      )
    )
  )
  
  logistf(
    f,
    data = data
  )
}


# ============================================================
# 14. Baseline predictor variables
# ============================================================

baseline_predictors <- c(
  "age",
  "disease_duration",
  "sex",
  "location_group",
  "l4",
  "behavior_group",
  "complex_fistula",
  "surgery_any",
  "prior_biologics_cat"
)


# ============================================================
# 15. Univariable Firth models
# ============================================================

baseline_uni <- map_dfr(
  
  baseline_predictors,
  
  function(v) {
    
    fit <- fit_firth(
      outcome =
        "w48_remission_mnri",
      
      predictors =
        v
    )
    
    extract_firth(
      fit,
      "Univariable"
    )
  }
)


# ============================================================
# 16. Multivariable Firth model
# ============================================================

baseline_multi_fit <- fit_firth(
  
  outcome =
    "w48_remission_mnri",
  
  predictors =
    baseline_predictors
)


baseline_multi <- extract_firth(
  baseline_multi_fit,
  "Multivariable"
)


# ============================================================
# 17. Academic labels
# ============================================================

term_labels <- c(
  
  "age" =
    "Age, per year",
  
  "disease_duration" =
    "Disease duration, per year",
  
  "sexFemale" =
    "Female vs male",
  
  "location_groupL3" =
    "Ileocolonic disease (L3) vs L1/L2",
  
  "l4Yes" =
    "L4 involvement: yes vs no",
  
  "behavior_groupB2" =
    "Stricturing phenotype (B2) vs B1",
  
  "behavior_groupPenetrating" =
    "Penetrating phenotype vs B1",
  
  "complex_fistulaYes" =
    "Complex fistula: yes vs no",
  
  "surgery_anyYes" =
    "Prior bowel surgery: yes vs no",
  
  "prior_biologics_cat3" =
    "3 prior biologic classes vs ≤2",
  
  "prior_biologics_cat≥4" =
    "≥4 prior biologic classes vs ≤2"
)


# ============================================================
# 18. SUMMARY 3
# Baseline predictors of W48 clinical remission
# ============================================================

summary_baseline_firth <- full_join(
  baseline_uni %>%
    select(Term, Uni_OR = OR, Uni_low = CI_low,
           Uni_high = CI_high, Uni_P = P_value),
  baseline_multi %>%
    select(Term, Adj_OR = OR, Adj_low = CI_low,
           Adj_high = CI_high, Adj_P = P_value),
  by = "Term"
) %>%
  mutate(
    Predictor = case_when(
      Term == "age" ~ "Age, per year",
      Term == "disease_duration" ~ "Disease duration, per year",
      Term == "sexFemale" ~ "Female vs male",
      Term == "location_groupL3" ~ "Ileocolonic disease (L3) vs L1/L2",
      Term == "l4Yes" ~ "L4 involvement: yes vs no",
      Term == "behavior_groupB2" ~ "Stricturing phenotype (B2) vs B1",
      Term == "behavior_groupPenetrating" ~ "Penetrating phenotype vs B1",
      Term == "complex_fistulaYes" ~ "Complex fistula: yes vs no",
      Term == "surgery_anyYes" ~ "Prior bowel surgery: yes vs no",
      Term == "prior_biologics_cat3" ~ "3 prior biologic classes vs ≤2",
      Term == "prior_biologics_cat≥4" ~ "≥4 prior biologic classes vs ≤2",
      TRUE ~ Term
    ),
    `Univariable OR (95% CI)` = ifelse(
      is.na(Uni_OR),
      "",
      format_or(Uni_OR, Uni_low, Uni_high)
    ),
    `Univariable P` = format_p(Uni_P),
    `Adjusted OR (95% CI)` = ifelse(
      is.na(Adj_OR),
      "",
      format_or(Adj_OR, Adj_low, Adj_high)
    ),
    `Adjusted P` = format_p(Adj_P)
  ) %>%
  select(
    Predictor,
    `Univariable OR (95% CI)`,
    `Univariable P`,
    `Adjusted OR (95% CI)`,
    `Adjusted P`,
    everything()
  )
# ============================================================
# 19. Early remission -> W48 remission
# Crude Firth
# ============================================================

early_crude_fit <- logistf(
  
  w48_remission_mnri ~
    early_remission_factor,
  
  data = dat
)


early_crude <- extract_firth(
  early_crude_fit,
  "Crude"
) %>%
  
  filter(
    Term ==
      "early_remission_factorEarly remission"
  )


# ============================================================
# 20. Early remission -> W48 remission
# Adjusted Firth
# ============================================================

early_adjusted_fit <- logistf(
  
  w48_remission_mnri ~
    
    early_remission_factor +
    
    age +
    
    disease_duration +
    
    sex +
    
    location_group +
    
    l4 +
    
    behavior_group +
    
    complex_fistula +
    
    surgery_any +
    
    prior_biologics_cat,
  
  data = dat
)


early_adjusted <- extract_firth(
  early_adjusted_fit,
  "Adjusted"
) %>%
  
  filter(
    Term ==
      "early_remission_factorEarly remission"
  )


# ============================================================
# 21. Actual W48 remission rate according to early remission
# ============================================================

early_group_rates <- dat %>%
  
  group_by(
    early_remission_factor
  ) %>%
  
  summarise(
    
    N = n(),
    
    Events =
      sum(
        w48_remission_mnri == 1
      ),
    
    .groups = "drop"
  ) %>%
  
  rowwise() %>%
  
  mutate(
    
    tmp =
      list(
        prop_ci(
          Events,
          N
        )
      ),
    
    Percent =
      tmp$Percent,
    
    CI_low =
      tmp$CI_low,
    
    CI_high =
      tmp$CI_high
  ) %>%
  
  ungroup() %>%
  
  select(
    -tmp
  ) %>%
  
  mutate(
    
    `Week-48 remission, n/N (%)` =
      sprintf(
        "%d/%d (%.1f%%)",
        Events,
        N,
        Percent
      ),
    
    `95% CI` =
      sprintf(
        "%.1f–%.1f%%",
        CI_low,
        CI_high
      )
  )


# ============================================================
# 22. SUMMARY 4
# Early remission predictive value
# ============================================================

summary_early_prediction <- bind_rows(
  early_crude,
  early_adjusted
) %>%
  
  mutate(
    
    Predictor =
      "Week-12 early clinical remission",
    
    `OR (95% CI)` =
      format_or(
        OR,
        CI_low,
        CI_high
      ),
    
    `P value` =
      format_p(
        P_value
      )
  ) %>%
  
  select(
    Model,
    Predictor,
    `OR (95% CI)`,
    `P value`,
    OR,
    CI_low,
    CI_high,
    P_value
  )


# ============================================================
# 23. Create one final summary object
# ============================================================

summary_48w <- list(
  
  `1. Clinical efficacy: mNRI and AO` =
    summary_efficacy,
  
  `2. Baseline Firth predictors of Week-48 clinical remission` =
    summary_baseline_firth,
  
  `3. Early and sustained clinical remission` =
    summary_durability,
  
  `4. Week-48 remission according to early remission status` =
    early_group_rates,
  
  `5. Early remission predicting Week-48 clinical remission` =
    summary_early_prediction
)


# ============================================================
# 24. PRINT EVERYTHING DIRECTLY IN R
# ============================================================

cat(
  "\n\n",
  "============================================================\n",
  "48-WEEK EFFICACY MANUSCRIPT SUMMARY\n",
  "============================================================\n"
)


# ------------------------------------------------------------
# SECTION 1
# ------------------------------------------------------------

cat(
  "\n\n",
  "1. CLINICAL RESPONSE AND REMISSION: mNRI AND AO\n",
  "------------------------------------------------------------\n"
)

print(
  summary_efficacy,
  n = Inf
)


# ------------------------------------------------------------
# SECTION 2
# ------------------------------------------------------------

cat(
  "\n\n",
  "2. BASELINE PREDICTORS OF WEEK-48 CLINICAL REMISSION\n",
  "   FIRTH PENALIZED LOGISTIC REGRESSION\n",
  "------------------------------------------------------------\n"
)

print(
  summary_baseline_firth %>%
    select(
      Predictor,
      `Univariable OR (95% CI)`,
      `Univariable P`,
      `Adjusted OR (95% CI)`,
      `Adjusted P`
    ),
  n = Inf
)


# ------------------------------------------------------------
# SECTION 3
# ------------------------------------------------------------

cat(
  "\n\n",
  "3. EARLY AND SUSTAINED CLINICAL REMISSION\n",
  "------------------------------------------------------------\n"
)

print(
  summary_durability,
  n = Inf
)


# ------------------------------------------------------------
# SECTION 4
# ------------------------------------------------------------

cat(
  "\n\n",
  "4. WEEK-48 REMISSION BY EARLY REMISSION STATUS\n",
  "------------------------------------------------------------\n"
)

print(
  early_group_rates %>%
    select(
      early_remission_factor,
      `Week-48 remission, n/N (%)`,
      `95% CI`
    ),
  n = Inf
)


# ------------------------------------------------------------
# SECTION 5
# ------------------------------------------------------------

cat(
  "\n\n",
  "5. EARLY REMISSION PREDICTING WEEK-48 CLINICAL REMISSION\n",
  "------------------------------------------------------------\n"
)

print(
  summary_early_prediction %>%
    select(
      Model,
      Predictor,
      `OR (95% CI)`,
      `P value`
    ),
  n = Inf
)


cat(
  "\n============================================================\n"
)


# ============================================================
# 25. Optional: directly open tables in RStudio Viewer
# ============================================================

 View(summary_efficacy)
 View(summary_baseline_firth)
 View(summary_durability)
 View(early_group_rates)
 View(summary_early_prediction)
 # 保存为CSV文件
 write.csv(summary_efficacy, "summary_efficacy.csv", row.names = FALSE)
 write.csv(summary_baseline_firth, "summary_baseline_firth.csv", row.names = FALSE)
 write.csv(summary_durability, "summary_durability.csv", row.names = FALSE)
 write.csv(early_group_rates, "early_group_rates.csv", row.names = FALSE)
 write.csv(summary_early_prediction, "summary_early_prediction.csv", row.names = FALSE)