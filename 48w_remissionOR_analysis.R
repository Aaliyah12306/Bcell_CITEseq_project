# ============================================================
# 48-WEEK EFFICACY ANALYSIS
# Prospective DTT-CD cohort
#
# Main analyses:
# 1. Week-48 efficacy: NRI + AO sensitivity analysis
# 2. Longitudinal clinical response/remission: Weeks 12/24/48
# 3. Early and sustained clinical remission
# 4. Baseline predictors of Week-48 clinical remission
#    using Firth penalized logistic regression
# 5. Week-12 early remission predicting Week-48 remission
# 6. Forest plots
#
# All figure/table labels are in academic English.
# ============================================================


# ============================================================
# 0. Packages
# ============================================================

packages <- c(
  "readxl",
  "dplyr",
  "tidyr",
  "stringr",
  "purrr",
  "ggplot2",
  "logistf",
  "car",
  "readr",
  "scales"
)
lapply(packages, library, character.only = TRUE)

# ============================================================
# 1. File paths
# ============================================================

input_file <- "efficacy_analysis.xlsx"
sheet_name <- "DTT-new"

output_dir <- "efficacy_48w_results"
table_dir  <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")

dir.create(output_dir, showWarnings = FALSE)
dir.create(table_dir, showWarnings = FALSE)
dir.create(figure_dir, showWarnings = FALSE)


# ============================================================
# 2. Helper functions
# ============================================================

# Convert mixed 0/1/"NA" columns into true binary numeric variables
as_binary <- function(x) {
  
  z <- trimws(as.character(x))
  
  z[
    z %in% c(
      "",
      "NA",
      "N/A",
      "na",
      "Na",
      "NaN"
    )
  ] <- NA_character_
  
  out <- suppressWarnings(as.integer(z))
  
  out[
    !is.na(out) &
      !out %in% c(0L, 1L)
  ] <- NA_integer_
  
  return(out)
}


# Exact binomial proportion + 95% CI
binom_summary <- function(events, total) {
  
  if (total == 0) {
    return(
      tibble(
        events = 0,
        total = 0,
        percent = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_
      )
    )
  }
  
  bt <- binom.test(
    x = events,
    n = total
  )
  
  tibble(
    events = events,
    total = total,
    percent = 100 * events / total,
    ci_low = 100 * bt$conf.int[1],
    ci_high = 100 * bt$conf.int[2]
  )
}


# P-value formatting
format_p <- function(p) {
  
  ifelse(
    is.na(p),
    "",
    ifelse(
      p < 0.001,
      "<0.001",
      sprintf("%.3f", p)
    )
  )
}


# OR formatting
format_or_ci <- function(or, low, high) {
  
  ifelse(
    is.na(or),
    "",
    sprintf(
      "%.2f (%.2f–%.2f)",
      or,
      low,
      high
    )
  )
}


# ============================================================
# 3. Read raw data
# ============================================================

raw <- read_excel(
  input_file,
  sheet = sheet_name,
  na = c("", "NA", "N/A")
)
dat<-raw

# ============================================================

# ============================================================


# ============================================================
# 6. Generate unique patient ID
#
# Cohort number alone is NOT globally unique across centers.
# ============================================================

# ============================================================
# 7. Clean baseline variables
# ============================================================

dat <- dat %>%
  mutate(
    
    # --------------------------------------------------------
    # Numeric variables
    # --------------------------------------------------------
    
    age =
      as.numeric(age),
    
    disease_duration =
      as.numeric(disease_duration),
    
    surgery_number =
      as.numeric(surgery_number),
    
    prior_biologics_n =
      as.numeric(prior_biologics_n),
    
    
    # --------------------------------------------------------
    # Sex
    # 1 = male; 2 = female
    # --------------------------------------------------------
    
    sex = factor(
      sex_code,
      levels = c(1, 2),
      labels = c(
        "Male",
        "Female"
      )
    ),
    
    
    # --------------------------------------------------------
    # Disease location
    #
    # Raw:
    # L1, L2, L3, L1+L4, L2+L4, L3+L4
    #
    # Separate primary location from L4 involvement.
    # --------------------------------------------------------
    
    location_raw =
      str_trim(
        as.character(location_raw)
      ),
    
    location_detail =
      str_extract(
        location_raw,
        "L[123]"
      ),
    
    location_detail = factor(
      location_detail,
      levels = c(
        "L1",
        "L2",
        "L3"
      )
    ),
    
    # Collapsed variable used in multivariable model
    # because isolated L2 group is very small.
    location_group = case_when(
      
      location_detail == "L1" ~ "L1",
      location_detail == "L2" ~ "L2",
      location_detail == "L3" ~ "L3",
      TRUE ~ NA_character_
    ),
    
    location_group = factor(
      location_group,
      levels = c(
        "L1",
        "L2",
        "L3"
      )
    ),
    
    
    # --------------------------------------------------------
    # L4
    # --------------------------------------------------------
    
    l4 = factor(
      l4_code,
      levels = c(0, 1),
      labels = c(
        "No",
        "Yes"
      )
    ),
    
    
    # --------------------------------------------------------
    # Disease behavior
    #
    # IMPORTANT:
    # trim trailing spaces:
    # B2 / B2  -> same category
    # B3 / B3  -> same category
    # --------------------------------------------------------
    
    behavior_clean =
      str_trim(
        as.character(behavior_raw)
      ),
    
    # Use safe model coding
    behavior_detail = case_when(
      
      behavior_clean == "B1" ~ "B1",
      
      behavior_clean == "B2" ~ "B2",
      
      behavior_clean == "B3" ~ "B3",
      
      behavior_clean == "B2+B3" ~ "B2B3",
      
      TRUE ~ NA_character_
    ),
    
    behavior_detail = factor(
      behavior_detail,
      levels = c(
        "B1",
        "B2",
        "B3",
        "B2B3"
      )
    ),
    
    # Parsimonious variable used in multivariable analysis
    behavior_group = case_when(
      
      behavior_detail == "B1" ~ "B1",
      
      behavior_detail == "B2" ~ "B2",
      
      behavior_detail %in%
        c("B3", "B2B3") ~ "B3",
      
      TRUE ~ NA_character_
    ),
    
    behavior_group = factor(
      behavior_group,
      levels = c(
        "B1",
        "B2",
        "B3"
      )
    ),
    
    
    # --------------------------------------------------------
    # Perianal disease
    # --------------------------------------------------------
    
    perianal = factor(
      perianal_code,
      levels = c(0, 1),
      labels = c(
        "No",
        "Yes"
      )
    ),
    
    
    # --------------------------------------------------------
    # Complex fistula
    # --------------------------------------------------------
    
    complex_fistula = factor(
      complex_fistula_code,
      levels = c(0, 1),
      labels = c(
        "No",
        "Yes"
      )
    ),
    
    
    # --------------------------------------------------------
    # Previous surgery
    # --------------------------------------------------------
    
    surgery_any = factor(
      surgery_any_code,
      levels = c(0, 1),
      labels = c(
        "No",
        "Yes"
      )
    ),
    
    surgery_cat = case_when(
      
      surgery_number == 0 ~ "S0",
      
      surgery_number == 1 ~ "S1",
      
      surgery_number >= 2 ~ "S2plus",
      
      TRUE ~ NA_character_
    ),
    
    surgery_cat = factor(
      surgery_cat,
      levels = c(
        "S0",
        "S1",
        "S2plus"
      )
    ),
    
    
    # --------------------------------------------------------
    # Number of prior biologic classes
    # Consistent with previous manuscript grouping
    # --------------------------------------------------------
    
    prior_biologics_cat = case_when(
      
      prior_biologics_n <= 2 ~ "LE2",
      
      prior_biologics_n == 3 ~ "N3",
      
      prior_biologics_n >= 4 ~ "GE4",
      
      TRUE ~ NA_character_
    ),
    
    prior_biologics_cat = factor(
      prior_biologics_cat,
      levels = c(
        "LE2",
        "N3",
        "GE4"
      )
    )
  )


# ============================================================
# 8. Clean efficacy variables
# ============================================================

dat <- dat %>%
  mutate(
    
    w12_response =
      as_binary(w12_response),
    
    w12_remission =
      as_binary(w12_remission),
    
    w24_response =
      as_binary(w24_response),
    
    w24_remission =
      as_binary(w24_remission),
    
    w48_response =
      as_binary(w48_response),
    
    w48_remission =
      as_binary(w48_remission)
  )


# ============================================================
# 9. Missingness check
# ============================================================

missing_table <- tibble(
  
  Variable = c(
    "Week 12 clinical response",
    "Week 12 clinical remission",
    "Week 24 clinical response",
    "Week 24 clinical remission",
    "Week 48 clinical response",
    "Week 48 clinical remission"
  ),
  
  Missing = c(
    sum(is.na(dat$w12_response)),
    sum(is.na(dat$w12_remission)),
    sum(is.na(dat$w24_response)),
    sum(is.na(dat$w24_remission)),
    sum(is.na(dat$w48_response)),
    sum(is.na(dat$w48_remission))
  ),
  
  Total = nrow(dat)
  
) %>%
  mutate(
    `Missing, %` =
      round(100 * Missing / Total, 1)
  )

print(missing_table)

readr::write_excel_csv(
  missing_table,
  file.path(
    table_dir,
    "QC_Missing_Efficacy_Data.csv"
  )
)


# ============================================================
# 10. Create NRI variables
#
# NRI:
# missing outcome = non-response / non-remission
# ============================================================

dat <- dat %>%
  mutate(
    
    w12_response_nri =
      replace_na(
        w12_response,
        0L
      ),
    
    w12_remission_nri =
      replace_na(
        w12_remission,
        0L
      ),
    
    w24_response_nri =
      replace_na(
        w24_response,
        0L
      ),
    
    w24_remission_nri =
      replace_na(
        w24_remission,
        0L
      ),
    
    w48_response_nri =
      replace_na(
        w48_response,
        0L
      ),
    
    w48_remission_nri =
      replace_na(
        w48_remission,
        0L
      )
  )


# ============================================================
# 11. Early and sustained remission
# ============================================================

dat <- dat %>%
  mutate(
    
    # Week-12 clinical remission
    early_remission =
      as.integer(
        w12_remission_nri == 1
      ),
    
    early_remission_factor = factor(
      early_remission,
      levels = c(0, 1),
      labels = c(
        "No",
        "Yes"
      )
    ),
    
    # Strict visit-based sustained remission:
    # remission at W12 + W24 + W48
    sustained_remission_nri =
      as.integer(
        w12_remission_nri == 1 &
          w24_remission_nri == 1 &
          w48_remission_nri == 1
      ),
    
    # Delayed remission (optional exploratory outcome)
    delayed_remission_nri =
      as.integer(
        w12_remission_nri == 0 &
          w48_remission_nri == 1
      ),
    
    # Loss of remission by W48
    loss_of_remission_nri =
      as.integer(
        w12_remission_nri == 1 &
          w48_remission_nri == 0
      )
  )


# ============================================================
# 12. Save cleaned analysis dataset
# ============================================================

readr::write_excel_csv(
  dat,
  file.path(
    output_dir,
    "analysis_dataset_clean.csv"
  )
)


# ============================================================
# 13. Supplementary Table S3
# Week-48 efficacy: NRI and AO
# ============================================================

make_endpoint_summary <- function(
    raw_var,
    nri_var,
    endpoint_name
) {
  
  # NRI
  x_nri <- sum(
    dat[[nri_var]] == 1,
    na.rm = TRUE
  )
  
  n_nri <- nrow(dat)
  
  nri <- binom_summary(
    x_nri,
    n_nri
  ) %>%
    mutate(
      Endpoint = endpoint_name,
      Analysis = "NRI"
    )
  
  # AO
  observed <- dat[[raw_var]]
  
  n_ao <- sum(
    !is.na(observed)
  )
  
  x_ao <- sum(
    observed == 1,
    na.rm = TRUE
  )
  
  ao <- binom_summary(
    x_ao,
    n_ao
  ) %>%
    mutate(
      Endpoint = endpoint_name,
      Analysis = "As observed"
    )
  
  bind_rows(
    nri,
    ao
  )
}


week48_table <- bind_rows(
  
  make_endpoint_summary(
    raw_var = "w48_response",
    nri_var = "w48_response_nri",
    endpoint_name = "Clinical response"
  ),
  
  make_endpoint_summary(
    raw_var = "w48_remission",
    nri_var = "w48_remission_nri",
    endpoint_name = "Clinical remission"
  )
  
) %>%
  mutate(
    
    `n/N (%)` = sprintf(
      "%d/%d (%.1f%%)",
      events,
      total,
      percent
    ),
    
    `95% CI` = sprintf(
      "%.1f–%.1f%%",
      ci_low,
      ci_high
    )
  ) %>%
  select(
    Endpoint,
    Analysis,
    events,
    total,
    percent,
    ci_low,
    ci_high,
    `n/N (%)`,
    `95% CI`
  )

print(week48_table)

readr::write_excel_csv(
  week48_table,
  file.path(
    table_dir,
    "Supplementary_Table_S3_Week48_Efficacy_NRI_AO.csv"
  )
)


# ============================================================
# 14. Supplementary Table S4
# Longitudinal clinical effectiveness
# ============================================================

efficacy_specs <- tibble::tribble(
  
  ~Visit, ~Endpoint, ~raw_var, ~nri_var,
  
  "Week 12",
  "Clinical response",
  "w12_response",
  "w12_response_nri",
  
  "Week 12",
  "Clinical remission",
  "w12_remission",
  "w12_remission_nri",
  
  "Week 24",
  "Clinical response",
  "w24_response",
  "w24_response_nri",
  
  "Week 24",
  "Clinical remission",
  "w24_remission",
  "w24_remission_nri",
  
  "Week 48",
  "Clinical response",
  "w48_response",
  "w48_response_nri",
  
  "Week 48",
  "Clinical remission",
  "w48_remission",
  "w48_remission_nri"
)


longitudinal_table <- purrr::pmap_dfr(
  
  efficacy_specs,
  
  function(
    Visit,
    Endpoint,
    raw_var,
    nri_var
  ) {
    
    # NRI
    nri <- binom_summary(
      events = sum(
        dat[[nri_var]] == 1
      ),
      total = nrow(dat)
    ) %>%
      mutate(
        Visit = Visit,
        Endpoint = Endpoint,
        Analysis = "NRI"
      )
    
    # AO
    raw_outcome <- dat[[raw_var]]
    
    ao <- binom_summary(
      events = sum(
        raw_outcome == 1,
        na.rm = TRUE
      ),
      total = sum(
        !is.na(raw_outcome)
      )
    ) %>%
      mutate(
        Visit = Visit,
        Endpoint = Endpoint,
        Analysis = "As observed"
      )
    
    bind_rows(
      nri,
      ao
    )
  }
  
) %>%
  mutate(
    
    Visit = factor(
      Visit,
      levels = c(
        "Week 12",
        "Week 24",
        "Week 48"
      )
    ),
    
    Endpoint = factor(
      Endpoint,
      levels = c(
        "Clinical response",
        "Clinical remission"
      )
    ),
    
    `n/N (%)` = sprintf(
      "%d/%d (%.1f%%)",
      events,
      total,
      percent
    ),
    
    `95% CI` = sprintf(
      "%.1f–%.1f%%",
      ci_low,
      ci_high
    )
  ) %>%
  arrange(
    Visit,
    Endpoint,
    Analysis
  )


readr::write_excel_csv(
  longitudinal_table,
  file.path(
    table_dir,
    "Supplementary_Table_S4_Longitudinal_Clinical_Effectiveness.csv"
  )
)


# ============================================================
# 15. Figure 2A
# Clinical response and remission through Week 48
# Primary figure uses NRI
# ============================================================

plot_long <- longitudinal_table %>%
  filter(
    Analysis == "NRI"
  )


efficacy_colors <- c(
  "Clinical response" = "#AFCBE8",
  "Clinical remission" = "#4D83B8"
)


p_efficacy <- ggplot(
  plot_long,
  aes(
    x = Visit,
    y = percent,
    fill = Endpoint
  )
) +
  
  geom_col(
    position = position_dodge(
      width = 0.76
    ),
    width = 0.66
  ) +
  
  geom_errorbar(
    aes(
      ymin = ci_low,
      ymax = ci_high
    ),
    position = position_dodge(
      width = 0.76
    ),
    width = 0.15,
    linewidth = 0.45
  ) +
  
  geom_text(
    aes(
      label = sprintf(
        "%d/%d\n(%.1f%%)",
        events,
        total,
        percent
      )
    ),
    position = position_dodge(
      width = 0.76
    ),
    vjust = -0.45,
    size = 3.4
  ) +
  
  scale_fill_manual(
    values = efficacy_colors
  ) +
  
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(
      0,
      100,
      by = 20
    ),
    labels = function(x) {
      paste0(x, "%")
    },
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  
  labs(
    title = "Clinical Effectiveness Through Week 48",
    x = NULL,
    y = "Patients, %",
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 12,
    base_family = "Arial"
  ) +
  
  theme(
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    
    legend.position = "top",
    
    legend.text = element_text(
      size = 10
    ),
    
    axis.text.x = element_text(
      face = "bold"
    )
  )


print(p_efficacy)


ggsave(
  file.path(
    figure_dir,
    "Figure_2A_Clinical_Effectiveness_Through_Week48.pdf"
  ),
  p_efficacy,
  width = 7,
  height = 5.5,
  device = cairo_pdf
)

ggsave(
  file.path(
    figure_dir,
    "Figure_2A_Clinical_Effectiveness_Through_Week48.tiff"
  ),
  p_efficacy,
  width = 7,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)

# ============================================================
# Figure 2A-AO
# Clinical response and remission through Week 48
# As-observed analysis
# No error bars
# ============================================================

plot_long_ao <- longitudinal_table %>%
  filter(
    Analysis == "As observed"
  )


efficacy_colors <- c(
  "Clinical response"  = "#AFCBE8",
  "Clinical remission" = "#4D83B8"
)


p_efficacy_ao <- ggplot(
  plot_long_ao,
  aes(
    x = Visit,
    y = percent,
    fill = Endpoint
  )
) +
  
  geom_col(
    position = position_dodge(
      width = 0.76
    ),
    width = 0.66
  ) +
  
  geom_text(
    aes(
      label = sprintf(
        "%d/%d\n(%.1f%%)",
        events,
        total,
        percent
      )
    ),
    position = position_dodge(
      width = 0.76
    ),
    vjust = -0.45,
    size = 3.4,
    family = "Arial"
  ) +
  
  scale_fill_manual(
    values = efficacy_colors
  ) +
  
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(
      0,
      100,
      by = 20
    ),
    labels = function(x) {
      paste0(x, "%")
    },
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  
  labs(
    title = "Clinical Effectiveness Through Week 48\n(As-Observed Analysis)",
    x = NULL,
    y = "Patients, %",
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 12,
    base_family = "Arial"
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13,
      lineheight = 1.1
    ),
    
    legend.position = "top",
    
    legend.text = element_text(
      size = 10
    ),
    
    axis.text.x = element_text(
      face = "bold"
    )
  )


print(p_efficacy_ao)


# ============================================================
# Save AO figure
# ============================================================

ggsave(
  filename = file.path(
    figure_dir,
    "Figure_Clinical_Effectiveness_Through_Week48_AO.pdf"
  ),
  plot = p_efficacy_ao,
  width = 7,
  height = 5.5,
  device = cairo_pdf
)


ggsave(
  filename = file.path(
    figure_dir,
    "Figure_Clinical_Effectiveness_Through_Week48_AO.tiff"
  ),
  plot = p_efficacy_ao,
  width = 7,
  height = 5.5,
  dpi = 600,
  compression = "lzw"
)
# ============================================================
# 16. Longitudinal remission statistical comparison
#
# NRI provides complete binary outcome data at all 3 visits.
#
# Overall:
# Cochran's Q
#
# Pairwise:
# exact McNemar tests with Holm adjustment
# ============================================================

remission_matrix <- as.matrix(
  dat %>%
    select(
      w12_remission_nri,
      w24_remission_nri,
      w48_remission_nri
    )
)


# ------------------------------------------------------------
# Cochran Q test calculated directly
# ------------------------------------------------------------

cochran_q_test <- function(mat) {
  
  k <- ncol(mat)
  
  column_totals <-
    colSums(mat)
  
  row_totals <-
    rowSums(mat)
  
  total_success <-
    sum(column_totals)
  
  denominator <-
    k * total_success -
    sum(row_totals^2)
  
  numerator <-
    (k - 1) *
    (
      k * sum(column_totals^2) -
        total_success^2
    )
  
  Q <- numerator / denominator
  
  df <- k - 1
  
  p_value <- pchisq(
    Q,
    df = df,
    lower.tail = FALSE
  )
  
  tibble(
    Test = "Cochran's Q test",
    Statistic = Q,
    df = df,
    P_value = p_value,
    P_formatted = format_p(p_value)
  )
}


cochran_result <-
  cochran_q_test(
    remission_matrix
  )


# ------------------------------------------------------------
# Exact paired McNemar
# ------------------------------------------------------------

exact_mcnemar <- function(
    x,
    y,
    comparison
) {
  
  discordant_01 <-
    sum(
      x == 0 &
        y == 1
    )
  
  discordant_10 <-
    sum(
      x == 1 &
        y == 0
    )
  
  discordant_total <-
    discordant_01 +
    discordant_10
  
  if (discordant_total == 0) {
    
    p_value <- 1
    
  } else {
    
    p_value <- binom.test(
      discordant_01,
      discordant_total,
      p = 0.5
    )$p.value
  }
  
  tibble(
    
    Comparison = comparison,
    
    Remission_from =
      100 * mean(x),
    
    Remission_to =
      100 * mean(y),
    
    Difference_percentage_points =
      100 * (
        mean(y) -
          mean(x)
      ),
    
    `Non-remission to remission` =
      discordant_01,
    
    `Remission to non-remission` =
      discordant_10,
    
    P_value =
      p_value
  )
}


pairwise_remission <- bind_rows(
  
  exact_mcnemar(
    dat$w12_remission_nri,
    dat$w24_remission_nri,
    "Week 12 vs Week 24"
  ),
  
  exact_mcnemar(
    dat$w24_remission_nri,
    dat$w48_remission_nri,
    "Week 24 vs Week 48"
  ),
  
  exact_mcnemar(
    dat$w12_remission_nri,
    dat$w48_remission_nri,
    "Week 12 vs Week 48"
  )
  
) %>%
  mutate(
    
    Holm_adjusted_P =
      p.adjust(
        P_value,
        method = "holm"
      ),
    
    P_formatted =
      format_p(
        P_value
      ),
    
    Holm_P_formatted =
      format_p(
        Holm_adjusted_P
      )
  )


readr::write_excel_csv(
  cochran_result,
  file.path(
    table_dir,
    "Supplementary_Table_S4A_Overall_Clinical_Remission_Comparison.csv"
  )
)

readr::write_excel_csv(
  pairwise_remission,
  file.path(
    table_dir,
    "Supplementary_Table_S4B_Pairwise_Clinical_Remission_Comparisons.csv"
  )
)


# ============================================================
# 17. Supplementary Table S5
# Early and sustained clinical remission
# ============================================================

n_total <- nrow(dat)

early_n <-
  sum(
    dat$early_remission == 1
  )

sustained_n <-
  sum(
    dat$sustained_remission_nri == 1
  )

early_summary <-
  binom_summary(
    early_n,
    n_total
  ) %>%
  mutate(
    Outcome = "Early clinical remission"
  )

sustained_summary <-
  binom_summary(
    sustained_n,
    n_total
  ) %>%
  mutate(
    Outcome = "Sustained clinical remission"
  )

# Conditional:
# sustained among early remitters
sustained_among_early <-
  binom_summary(
    sustained_n,
    early_n
  ) %>%
  mutate(
    Outcome =
      "Sustained remission among early remitters"
  )


durability_table <- bind_rows(
  
  early_summary,
  sustained_summary,
  sustained_among_early
  
) %>%
  mutate(
    
    `n/N (%)` = sprintf(
      "%d/%d (%.1f%%)",
      events,
      total,
      percent
    ),
    
    `95% CI` = sprintf(
      "%.1f–%.1f%%",
      ci_low,
      ci_high
    )
  )


readr::write_excel_csv(
  durability_table,
  file.path(
    table_dir,
    "Supplementary_Table_S5_Early_and_Sustained_Remission.csv"
  )
)


# ============================================================
# 18. Figure 2B
# Early and sustained remission
# Both percentages use the total cohort denominator
# ============================================================

durability_plot_dat <-
  bind_rows(
    early_summary,
    sustained_summary
  ) %>%
  mutate(
    Outcome = factor(
      Outcome,
      levels = c(
        "Early clinical remission",
        "Sustained clinical remission"
      )
    )
  )


p_durability <- ggplot(
  durability_plot_dat,
  aes(
    x = Outcome,
    y = percent,
    fill = Outcome
  )
) +
  
  geom_col(
    width = 0.58
  ) +
  
  geom_errorbar(
    aes(
      ymin = ci_low,
      ymax = ci_high
    ),
    width = 0.13,
    linewidth = 0.45
  ) +
  
  geom_text(
    aes(
      label = sprintf(
        "%d/%d\n(%.1f%%)",
        events,
        total,
        percent
      )
    ),
    vjust = -0.45,
    size = 3.5
  ) +
  
  scale_fill_manual(
    values = c(
      "#7EADD8",
      "#4F76A3"
    )
  ) +
  
  scale_y_continuous(
    limits = c(0, 70),
    labels = function(x) {
      paste0(x, "%")
    }
  ) +
  
  labs(
    title = "Durability of Clinical Remission",
    x = NULL,
    y = "Patients, %",
    fill = NULL
  ) +
  
  theme_classic(
    base_size = 12,
    base_family = "Arial"
  ) +
  
  theme(
    
    legend.position = "none",
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    
    axis.text.x = element_text(
      face = "bold"
    )
  )


print(p_durability)


ggsave(
  file.path(
    figure_dir,
    "Figure_2B_Durability_of_Clinical_Remission.pdf"
  ),
  p_durability,
  width = 6,
  height = 5,
  device = cairo_pdf
)

ggsave(
  file.path(
    figure_dir,
    "Figure_2B_Durability_of_Clinical_Remission.tiff"
  ),
  p_durability,
  width = 6,
  height = 5,
  dpi = 600,
  compression = "lzw"
)


# ============================================================
# 19. Firth logistic regression functions
# ============================================================

extract_firth <- function(
    fit,
    model_name
) {
  
  tibble(
    
    term =
      names(
        fit$coefficients
      ),
    
    log_OR =
      as.numeric(
        fit$coefficients
      ),
    
    log_CI_low =
      as.numeric(
        fit$ci.lower
      ),
    
    log_CI_high =
      as.numeric(
        fit$ci.upper
      ),
    
    P_value =
      as.numeric(
        fit$prob
      )
    
  ) %>%
    
    filter(
      term != "(Intercept)"
    ) %>%
    
    mutate(
      
      OR =
        exp(log_OR),
      
      CI_low =
        exp(log_CI_low),
      
      CI_high =
        exp(log_CI_high),
      
      Model =
        model_name
    )
}


fit_firth_model <- function(
    data,
    outcome,
    predictors
) {
  
  formula_string <- paste(
    outcome,
    "~",
    paste(
      predictors,
      collapse = " + "
    )
  )
  
  logistf(
    as.formula(
      formula_string
    ),
    data = data
  )
}


# ============================================================
# 20. Baseline predictor model specification
#
# IMPORTANT:
# Parsimonious model:
#
# - age
# - disease duration
# - sex
# - L3 vs L1/L2
# - L4
# - disease behavior
# - complex fistula
# - prior surgery
# - prior biologic classes
#
# We intentionally do NOT include both:
# perianal disease + complex fistula
# surgery history + surgery count
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
# 21. Human-readable regression labels
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
    "Upper gastrointestinal involvement (L4): yes vs no",
  
  "behavior_groupB2" =
    "Stricturing phenotype (B2) vs inflammatory phenotype (B1)",
  
  "behavior_groupPenetrating" =
    "Penetrating phenotype vs inflammatory phenotype (B1)",
  
  "complex_fistulaYes" =
    "Complex fistula: yes vs no",
  
  "surgery_anyYes" =
    "Prior bowel surgery: yes vs no",
  
  "prior_biologics_catN3" =
    "3 prior biologic classes vs ≤2",
  
  "prior_biologics_catGE4" =
    "≥4 prior biologic classes vs ≤2"
)


# ============================================================
# 22. Univariable Firth regression
# Primary outcome = Week-48 remission by NRI
# ============================================================

univariable_results <- purrr::map_dfr(
  
  baseline_predictors,
  
  function(variable) {
    
    fit <- fit_firth_model(
      data = dat,
      outcome = "w48_remission_nri",
      predictors = variable
    )
    
    extract_firth(
      fit,
      model_name = "Univariable"
    )
  }
)


# ============================================================
# 23. Multivariable Firth regression
# ============================================================

fit_baseline_multi <- fit_firth_model(
  
  data = dat,
  
  outcome = "w48_remission_nri",
  
  predictors = baseline_predictors
)


multivariable_results <- extract_firth(
  
  fit_baseline_multi,
  
  model_name =
    "Multivariable"
)


# ============================================================
# 24. Supplementary Table S6A
# Baseline categorical remission rates
# ============================================================

categorical_variables <- list(
  
  list(
    var = "sex",
    label = "Sex"
  ),
  
  list(
    var = "location_group",
    label = "Disease location"
  ),
  
  list(
    var = "l4",
    label = "L4 involvement"
  ),
  
  list(
    var = "behavior_group",
    label = "Disease behavior"
  ),
  
  list(
    var = "complex_fistula",
    label = "Complex fistula"
  ),
  
  list(
    var = "surgery_any",
    label = "Prior bowel surgery"
  ),
  
  list(
    var = "prior_biologics_cat",
    label = "Prior biologic classes"
  )
)


make_group_rate <- function(
    data,
    variable,
    variable_label
) {
  
  data %>%
    
    group_by(
      Level =
        .data[[variable]]
    ) %>%
    
    summarise(
      
      N = n(),
      
      Remission =
        sum(
          w48_remission_nri == 1
        ),
      
      .groups = "drop"
    ) %>%
    
    mutate(
      
      Variable =
        variable_label,
      
      Percent =
        100 * Remission / N,
      
      CI =
        purrr::map2(
          Remission,
          N,
          ~binom.test(.x, .y)$conf.int
        ),
      
      CI_low =
        100 *
        purrr::map_dbl(
          CI,
          1
        ),
      
      CI_high =
        100 *
        purrr::map_dbl(
          CI,
          2
        ),
      
      `Remission, n/N (%)` =
        sprintf(
          "%d/%d (%.1f%%)",
          Remission,
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
      Variable,
      Level,
      N,
      Remission,
      Percent,
      CI_low,
      CI_high,
      `Remission, n/N (%)`,
      `95% CI`
    )
}


baseline_rates <- purrr::map_dfr(
  
  categorical_variables,
  
  function(x) {
    
    make_group_rate(
      data = dat,
      variable = x$var,
      variable_label = x$label
    )
  }
)

readr::write_excel_csv(
  baseline_rates,
  file.path(
    table_dir,
    "Supplementary_Table_S6A_Baseline_Remission_Rates.csv"
  )
)


# ============================================================
# 25. Supplementary Table S6B
# Univariable + multivariable Firth regression
# ============================================================

regression_table <- full_join(
  
  univariable_results %>%
    select(
      term,
      uni_OR = OR,
      uni_CI_low = CI_low,
      uni_CI_high = CI_high,
      uni_P = P_value
    ),
  
  multivariable_results %>%
    select(
      term,
      adj_OR = OR,
      adj_CI_low = CI_low,
      adj_CI_high = CI_high,
      adj_P = P_value
    ),
  
  by = "term"
) %>%
  
  mutate(
    
    Predictor =
      dplyr::recode(
        term,
        !!!term_labels,
        .default = term
      ),
    
    `Univariable OR (95% CI)` =
      format_or_ci(
        uni_OR,
        uni_CI_low,
        uni_CI_high
      ),
    
    `Univariable P` =
      format_p(
        uni_P
      ),
    
    `Adjusted OR (95% CI)` =
      format_or_ci(
        adj_OR,
        adj_CI_low,
        adj_CI_high
      ),
    
    `Adjusted P` =
      format_p(
        adj_P
      )
  ) %>%
  
  select(
    
    Predictor,
    
    uni_OR,
    uni_CI_low,
    uni_CI_high,
    uni_P,
    
    adj_OR,
    adj_CI_low,
    adj_CI_high,
    adj_P,
    
    `Univariable OR (95% CI)`,
    `Univariable P`,
    `Adjusted OR (95% CI)`,
    `Adjusted P`
  )


readr::write_excel_csv(
  regression_table,
  file.path(
    table_dir,
    "Supplementary_Table_S6B_Firth_Regression_Week48_Clinical_Remission.csv"
  )
)


# ============================================================
# 26. Multicollinearity diagnostic
#
# Use ordinary logistic glm ONLY for VIF calculation.
# Effect estimates still come from Firth regression.
# ============================================================

vif_formula <- as.formula(
  paste(
    "w48_remission_nri ~",
    paste(
      baseline_predictors,
      collapse = " + "
    )
  )
)


glm_vif <- glm(
  vif_formula,
  family = binomial(),
  data = dat
)


vif_raw <- car::vif(
  glm_vif
)


if (is.matrix(vif_raw)) {
  
  vif_table <- data.frame(
    Term = rownames(vif_raw),
    vif_raw,
    row.names = NULL
  )
  
  if (
    all(
      c(
        "GVIF",
        "Df"
      ) %in% colnames(vif_table)
    )
  ) {
    
    vif_table$Adjusted_GVIF <-
      vif_table$GVIF^(
        1 /
          (
            2 *
              vif_table$Df
          )
      )
  }
  
} else {
  
  vif_table <- tibble(
    Term = names(vif_raw),
    VIF = as.numeric(vif_raw)
  )
}


readr::write_excel_csv(
  vif_table,
  file.path(
    table_dir,
    "Supplementary_Table_S6C_Multicollinearity_Diagnostics.csv"
  )
)


# ============================================================
# 27. Figure 5A
# Forest plot:
# adjusted baseline predictors of Week-48 remission
# ============================================================

forest_baseline <- multivariable_results %>%
  
  mutate(
    
    Predictor =
      dplyr::recode(
        term,
        !!!term_labels,
        .default = term
      ),
    
    Predictor = factor(
      Predictor,
      levels = rev(
        dplyr::recode(
          term,
          !!!term_labels,
          .default = term
        )
      )
    )
  )


p_forest_baseline <- ggplot(
  forest_baseline,
  aes(
    x = OR,
    y = Predictor
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = 2,
    linewidth = 0.5,
    color = "grey50"
  ) +
  
  geom_segment(
    aes(
      x = CI_low,
      xend = CI_high,
      yend = Predictor
    ),
    linewidth = 0.65
  ) +
  
  geom_point(
    size = 2.7
  ) +
  
  scale_x_log10() +
  
  labs(
    title = "Baseline Factors Associated With Week-48 Clinical Remission",
    x = "Adjusted Odds Ratio (95% CI)",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 11,
    base_family = "Arial"
  ) +
  
  theme(
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 12
    ),
    
    axis.text.y = element_text(
      size = 9
    )
  )


print(p_forest_baseline)


ggsave(
  file.path(
    figure_dir,
    "Figure_5A_Forest_Baseline_Predictors_Week48_Remission.pdf"
  ),
  p_forest_baseline,
  width = 8.5,
  height = 5.8,
  device = cairo_pdf,
)

ggsave(
  file.path(
    figure_dir,
    "Figure_5A_Forest_Baseline_Predictors_Week48_Remission.tiff"
  ),
  p_forest_baseline,
  width = 8.5,
  height = 5.8,
  dpi = 600,
  compression = "lzw"
)


# ============================================================
# 28. AO sensitivity analysis
# Baseline predictors of Week-48 remission
# ============================================================

dat_w48_ao <- dat %>%
  filter(
    !is.na(
      w48_remission
    )
  )


# Univariable
univariable_ao <- purrr::map_dfr(
  
  baseline_predictors,
  
  function(variable) {
    
    fit <- fit_firth_model(
      data = dat_w48_ao,
      outcome = "w48_remission",
      predictors = variable
    )
    
    extract_firth(
      fit,
      model_name = "AO univariable"
    )
  }
)


# Multivariable
fit_multi_ao <- fit_firth_model(
  
  data = dat_w48_ao,
  
  outcome = "w48_remission",
  
  predictors = baseline_predictors
)


multi_ao <- extract_firth(
  fit_multi_ao,
  "AO multivariable"
)


ao_regression_table <- full_join(
  univariable_ao %>%
    select(term, uni_OR = OR, uni_CI_low = CI_low, 
           uni_CI_high = CI_high, uni_P = P_value),
  multi_ao %>%
    select(term, adj_OR = OR, adj_CI_low = CI_low,
           adj_CI_high = CI_high, adj_P = P_value),
  by = "term"
) %>%
  mutate(
    Predictor = case_when(
      term == "age" ~ "Age, per year",
      term == "disease_duration" ~ "Disease duration, per year",
      term == "sexFemale" ~ "Female vs male",
      term == "location_groupL3" ~ "Ileocolonic disease (L3) vs L1/L2",
      term == "l4Yes" ~ "Upper gastrointestinal involvement (L4): yes vs no",
      term == "behavior_groupB2" ~ "Stricturing phenotype (B2) vs inflammatory phenotype (B1)",
      term == "behavior_groupPenetrating" ~ "Penetrating phenotype vs inflammatory phenotype (B1)",
      term == "complex_fistulaYes" ~ "Complex fistula: yes vs no",
      term == "surgery_anyYes" ~ "Prior bowel surgery: yes vs no",
      term == "prior_biologics_catN3" ~ "3 prior biologic classes vs ≤2",
      term == "prior_biologics_catGE4" ~ "≥4 prior biologic classes vs ≤2",
      TRUE ~ term
    ),
    `Univariable OR (95% CI)` = format_or_ci(uni_OR, uni_CI_low, uni_CI_high),
    `Univariable P` = format_p(uni_P),
    `Adjusted OR (95% CI)` = format_or_ci(adj_OR, adj_CI_low, adj_CI_high),
    `Adjusted P` = format_p(adj_P)
  )


readr::write_excel_csv(
  ao_regression_table,
  file.path(
    table_dir,
    "Supplementary_Table_S6D_Firth_Regression_AO_Sensitivity.csv"
  )
)


# ============================================================
# 29. Early remission -> Week-48 remission
#
# Model 1: crude
# Model 2: adjusted
# ============================================================

fit_early_crude <- logistf(
  
  w48_remission_nri ~
    early_remission_factor,
  
  data = dat
)


fit_early_adjusted <- logistf(
  
  w48_remission_nri ~
    
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


early_crude_results <-
  extract_firth(
    fit_early_crude,
    "Crude"
  ) %>%
  filter(
    term ==
      "early_remission_factorYes"
  )


early_adjusted_results <-
  extract_firth(
    fit_early_adjusted,
    "Adjusted"
  ) %>%
  filter(
    term ==
      "early_remission_factorYes"
  )


early_prediction <- bind_rows(
  
  early_crude_results,
  early_adjusted_results
  
) %>%
  
  mutate(
    
    Predictor =
      "Week-12 early clinical remission",
    
    `OR (95% CI)` =
      format_or_ci(
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
    OR,
    CI_low,
    CI_high,
    P_value,
    `OR (95% CI)`,
    `P value`
  )


readr::write_excel_csv(
  early_prediction,
  file.path(
    table_dir,
    "Supplementary_Table_S7A_Firth_Early_Remission_Prediction.csv"
  )
)


# ============================================================
# 30. Week-48 remission proportion by early remission status
# ============================================================

early_group_rates <- dat %>%
  
  group_by(
    early_remission_factor
  ) %>%
  
  summarise(
    
    N = n(),
    
    Remission =
      sum(
        w48_remission_nri == 1
      ),
    
    .groups = "drop"
  ) %>%
  
  mutate(
    
    Percent =
      100 * Remission / N,
    
    CI =
      purrr::map2(
        Remission,
        N,
        ~binom.test(.x, .y)$conf.int
      ),
    
    CI_low =
      100 *
      purrr::map_dbl(
        CI,
        1
      ),
    
    CI_high =
      100 *
      purrr::map_dbl(
        CI,
        2
      ),
    
    `Week-48 remission, n/N (%)` =
      sprintf(
        "%d/%d (%.1f%%)",
        Remission,
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


readr::write_excel_csv(
  early_group_rates,
  file.path(
    table_dir,
    "Supplementary_Table_S7B_Week48_Remission_by_Early_Remission_Status.csv"
  )
)


# ============================================================
# 31. AO sensitivity:
# Early remission predicting W48 remission
# ============================================================

fit_early_crude_ao <- logistf(
  
  w48_remission ~
    early_remission_factor,
  
  data = dat_w48_ao
)


fit_early_adjusted_ao <- logistf(
  
  w48_remission ~
    
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
  
  data = dat_w48_ao
)


early_ao <- bind_rows(
  
  extract_firth(
    fit_early_crude_ao,
    "AO crude"
  ) %>%
    filter(
      term ==
        "early_remission_factorYes"
    ),
  
  extract_firth(
    fit_early_adjusted_ao,
    "AO adjusted"
  ) %>%
    filter(
      term ==
        "early_remission_factorYes"
    )
  
) %>%
  
  mutate(
    
    Predictor =
      "Week-12 early clinical remission",
    
    `OR (95% CI)` =
      format_or_ci(
        OR,
        CI_low,
        CI_high
      ),
    
    `P value` =
      format_p(
        P_value
      )
  )


readr::write_excel_csv(
  early_ao,
  file.path(
    table_dir,
    "Supplementary_Table_S7C_Early_Remission_AO_Sensitivity.csv"
  )
)


# ============================================================
# 32. Figure 5B
# Forest plot:
# Early remission -> Week-48 remission
# Crude and adjusted estimates
# ============================================================

forest_early <- bind_rows(
  
  early_crude_results,
  early_adjusted_results
  
) %>%
  
  mutate(
    
    Model = factor(
      Model,
      levels = c(
        "Adjusted",
        "Crude"
      )
    )
  )


p_forest_early <- ggplot(
  forest_early,
  aes(
    x = OR,
    y = Model
  )
) +
  
  geom_vline(
    xintercept = 1,
    linetype = 2,
    linewidth = 0.5,
    color = "grey50"
  ) +
  
  geom_segment(
    aes(
      x = CI_low,
      xend = CI_high,
      yend = Model
    ),
    linewidth = 0.7
  ) +
  
  geom_point(
    size = 3
  ) +
  
  scale_x_log10() +
  
  labs(
    title = "Association Between Early and Week-48 Clinical Remission",
    x = "Odds Ratio (95% CI)",
    y = NULL
  ) +
  
  theme_classic(
    base_size = 12,
    base_family = "Arial"
  ) +
  
  theme(
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 12
    )
  )


print(p_forest_early)


ggsave(
  file.path(
    figure_dir,
    "Figure_5B_Forest_Early_Remission_Prediction.pdf"
  ),
  p_forest_early,
  width = 6.5,
  height = 3.5,
  device = cairo_pdf
)

ggsave(
  file.path(
    figure_dir,
    "Figure_5B_Forest_Early_Remission_Prediction.tiff"
  ),
  p_forest_early,
  width = 6.5,
  height = 3.5,
  dpi = 600,
  compression = "lzw"
)


# ============================================================
# 33. Simple remission transition table
# Useful for manuscript / supplementary material
# ============================================================

transition_remission <- dat %>%
  
  count(
    
    W12 =
      w12_remission_nri,
    
    W24 =
      w24_remission_nri,
    
    W48 =
      w48_remission_nri,
    
    name = "n"
  ) %>%
  
  mutate(
    
    W12 =
      if_else(
        W12 == 1,
        "Remission",
        "Non-remission"
      ),
    
    W24 =
      if_else(
        W24 == 1,
        "Remission",
        "Non-remission"
      ),
    
    W48 =
      if_else(
        W48 == 1,
        "Remission",
        "Non-remission"
      ),
    
    Percent =
      100 * n / nrow(dat)
  ) %>%
  
  arrange(
    desc(n)
  )


readr::write_excel_csv(
  transition_remission,
  file.path(
    table_dir,
    "Supplementary_Table_S5B_Clinical_Remission_Transitions.csv"
  )
)


# ============================================================
# 34. Print core manuscript results
# ============================================================

cat(
  "\n========================================\n"
)

cat(
  "CORE 48-WEEK EFFICACY RESULTS\n"
)

cat(
  "========================================\n"
)


print(
  week48_table
)


cat(
  "\nEarly remission:\n"
)

print(
  early_summary
)


cat(
  "\nSustained remission:\n"
)

print(
  sustained_summary
)


cat(
  "\nSustained remission among early remitters:\n"
)

print(
  sustained_among_early
)


cat(
  "\nLongitudinal clinical remission comparison:\n"
)

print(
  cochran_result
)


cat(
  "\nPairwise comparisons:\n"
)

print(
  pairwise_remission
)


cat(
  "\nBaseline Firth regression:\n"
)

print(
  regression_table
)


cat(
  "\nEarly remission prediction:\n"
)

print(
  early_prediction
)


# ============================================================
# 35. Save model summaries
# ============================================================

capture.output(
  
  summary(
    fit_baseline_multi
  ),
  
  file = file.path(
    output_dir,
    "Firth_Model_Baseline_Predictors_Summary.txt"
  )
)


capture.output(
  
  summary(
    fit_early_adjusted
  ),
  
  file = file.path(
    output_dir,
    "Firth_Model_Early_Remission_Adjusted_Summary.txt"
  )
)


# ============================================================
# 36. Save R session information
# Essential for reproducibility
# ============================================================

capture.output(
  
  sessionInfo(),
  
  file = file.path(
    output_dir,
    "sessionInfo.txt"
  )
)


cat(
  "\n========================================\n"
)

cat(
  "Analysis completed successfully.\n"
)

cat(
  paste0(
    "Results saved to: ",
    normalizePath(
      output_dir
    ),
    "\n"
  )
)

cat(
  "========================================\n"
)