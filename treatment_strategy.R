# ============================================================
# 0. Packages
# ============================================================

library(readxl)
library(dplyr)
library(stringr)
library(tidyr)
library(ggplot2)
library(survival)
library(survminer)
library(ggalluvial)
library(scales)

# ============================================================
# 1. Read data
# ============================================================

file_path <- "treatment_strategy.xlsx"

dat <- read_excel(
  file_path,
  sheet = "SAE"
)

# Check
dim(dat)
head(dat)
str(dat)
# ============================================================
# 2. Data preparation
# ============================================================

dat <- dat %>%
  mutate(
    # "48w" -> 48
    drug_time_num = as.numeric(
      str_extract(drug_time, "\\d+\\.?\\d*")
    ),
    
    drug_event = as.integer(drug_event),
    stop_reason = as.integer(stop_reason),
    
    strategy24_plot = as.integer(strategy24_plot),
    strategy48_plot = as.integer(strategy48_plot)
  )

# Check
summary(dat$drug_time_num)
table(dat$drug_event)
table(dat$stop_reason, useNA = "ifany")

# ============================================================
# 3. Data consistency checks
# ============================================================

# Event = 1 should have a stop reason
dat %>%
  filter(drug_event == 1 & is.na(stop_reason))

# Event = 0 should not have a stop reason
dat %>%
  filter(drug_event == 0 & !is.na(stop_reason))

# All discontinued patients should be "Discontinued" at Week 48
dat %>%
  filter(drug_event == 1 & strategy48_plot != 3)

# Patients censored before Week 48 should generally be Lost/unknown
dat %>%
  filter(
    drug_event == 0,
    drug_time_num < 48,
    strategy48_plot != 4
  )

# ============================================================
# 4. Kaplan-Meier analysis of drug persistence
# ============================================================

km_drug <- survfit(
  Surv(
    time = drug_time_num,
    event = drug_event
  ) ~ 1,
  data = dat
)

summary(km_drug)
# ============================================================
# 5. Kaplan-Meier plot
# ============================================================

km_plot <- ggsurvplot(
  km_drug,
  data = dat,
  
  # 95% CI
  conf.int = TRUE,
  
  # Show censor marks
  censor = TRUE,
  censor.shape = "|",
  censor.size = 3,
  
  # Risk table
  risk.table = TRUE,
  risk.table.height = 0.23,
  risk.table.y.text = FALSE,
  
  # Axis
  xlim = c(0, 48),
  break.time.by = 12,
  
  ylim = c(0, 1),
  
  xlab = "Weeks Since Treatment Initiation",
  ylab = "Probability of Treatment Persistence",
  
  # Figure
  title = "Upadacitinib Persistence Through Week 48",
  
  # No legend needed for single-arm KM
  legend = "none",
  
  # Academic formatting
  surv.scale = "percent",
  
  ggtheme = theme_classic(base_size = 12),
  
  tables.theme = theme_cleantable()
)

km_plot
# ============================================================
# Week 48 KM estimate
# ============================================================

km48 <- summary(
  km_drug,
  times = 48
)

km48_est <- km48$surv
km48_low <- km48$lower
km48_high <- km48$upper

cat(
  sprintf(
    "Week 48 drug persistence: %.1f%% (95%% CI %.1f–%.1f%%)\n",
    km48_est * 100,
    km48_low * 100,
    km48_high * 100
  )
)
# Save combined KM plot + risk table
pdf(
  "Figure_Drug_Persistence_KM.pdf",
  width = 8,
  height = 8
)

print(km_plot)

dev.off()
tiff(
  "Figure_Drug_Persistence_KM.tiff",
  width = 7,
  height = 6,
  units = "in",
  res = 600,
  compression = "lzw"
)

print(km_plot)

dev.off()

##sanley
# ============================================================
# ============================================================
# 6. Treatment strategy labels
# ============================================================

strategy_labels <- c(
  "0" = "UPA 15 mg",
  "1" = "UPA 30 mg",
  "2" = "UPA + advanced therapy",
  "3" = "Discontinued",
  "4" = "Lost/unknown"
)

dat <- dat %>%
  mutate(
    # Week 12: all patients are on UPA 15 mg maintenance
    strategy12_plot = 0,
    
    strategy12 = factor(
      strategy12_plot,
      levels = 0:4,
      labels = strategy_labels
    ),
    
    # Week 24: strategy applied during Weeks 24–48
    strategy24 = factor(
      strategy24_plot,
      levels = 0:4,
      labels = strategy_labels
    ),
    
    # Week 48: strategy applied beyond Week 48
    strategy48 = factor(
      strategy48_plot,
      levels = 0:4,
      labels = strategy_labels
    ),
    
    # A duplicated column only for visual extension
    strategy48_ext = factor(
      strategy48_plot,
      levels = 0:4,
      labels = strategy_labels
    )
  )
# ============================================================
# 7. Sankey / alluvial plot
# ============================================================

sankey_dat <- dat %>%
  select(
    strategy12,
    strategy24,
    strategy48
  )
strategy_colors <- c(
  "UPA 15 mg" = "#76A9E0",                # creamy blue
  "UPA 30 mg" = "#F2AD5C",                # creamy orange
  "UPA + advanced therapy" = "#A987D2",   # creamy purple
  "Discontinued" = "#E77F82",             # creamy coral
  "Lost/unknown" = "#C8BBAA"              # creamy warm gray
)
p_sankey <- ggplot(
  sankey_dat,
  aes(
    axis1 = strategy12,
    axis2 = strategy24,
    axis3 = strategy48,
    y = 1
  )
) +
  
  # ==========================================================
# Flows
# ==========================================================
geom_alluvium(
  aes(fill = strategy24),
  width = 0.08,
  alpha = 0.82,
  knot.pos = 0.45,
  color = NA
) +
  
  # ==========================================================
# Nodes
# ==========================================================
geom_stratum(
  aes(fill = after_stat(stratum)),
  width = 0.08,
  color = "white",
  linewidth = 0.35
) +
  
  # ==========================================================
# Colors
# ==========================================================
scale_fill_manual(
  values = strategy_colors,
  drop = FALSE
) +
  
  # ==========================================================
# X-axis = time points, not stages
# ==========================================================
scale_x_discrete(
  limits = c("Week 24", "Week 48", ""),
  expand = c(0.03, 0.03)
) +
  
  # ==========================================================
# Labels
# ==========================================================
labs(
  title = "Maintenance Treatment Strategy Transitions",
  x = NULL,
  y = NULL,
  fill = "Treatment Strategy"
) +
  
  theme_classic(base_size = 12) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 13
    ),
    
    axis.text.x = element_text(
      face = "bold",
      size = 11,
      margin = margin(t = 8)
    ),
    
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    
    axis.ticks.x = element_blank(),
    
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 10
    ),
    
    legend.text = element_text(
      size = 9
    ),
    
    plot.margin = margin(20, 25, 35, 20)
  ) +
  
  # ==========================================================
# Stage labels between time points
# ==========================================================
annotate(
  "text",
  x = 1.5, y = -8,
  label = "Weeks 24–48",
  size = 3.6,
  fontface = "bold"
) +
  annotate(
    "text",
    x = 2.5, y = -8,
    label = "Beyond Week 48",
    size = 3.6,
    fontface = "bold"
  ) +
  
  coord_cartesian(clip = "off")

p_sankey

ggsave(
  filename = "Figure_Treatment_Strategy_Sankey.pdf",
  plot = p_sankey,
  width = 8,
  height = 6
)

ggsave(
  filename = "Figure_Treatment_Strategy_Sankey.tiff",
  plot = p_sankey,
  width = 8,
  height = 6,
  units = "in",
  dpi = 600,
  compression = "lzw"
)
# ============================================================
# 8. Supplementary table:
#    Reasons for permanent discontinuation
# ============================================================

n_total <- nrow(dat)

n_discontinued <- sum(
  dat$drug_event == 1,
  na.rm = TRUE
)

stop_reason_table <- dat %>%
  filter(drug_event == 1) %>%
  
  mutate(
    stop_reason_label = case_when(
      stop_reason == 1 ~ "Inefficacy",
      stop_reason == 2 ~ "Adverse events",
      stop_reason == 3 ~ "Pregnancy",
      TRUE ~ "Other"
    )
  ) %>%
  
  count(
    stop_reason_label,
    name = "n"
  ) %>%
  
  mutate(
    pct_discontinuations =
      100 * n / n_discontinued,
    
    pct_total_cohort =
      100 * n / n_total,
    
    `n (% of discontinuations)` =
      sprintf(
        "%d (%.1f%%)",
        n,
        pct_discontinuations
      ),
    
    `n (% of total cohort)` =
      sprintf(
        "%d (%.1f%%)",
        n,
        pct_total_cohort
      )
  ) %>%
  
  select(
    `Reason for Discontinuation` = stop_reason_label,
    `n (% of discontinuations)`,
    `n (% of total cohort)`
  )

stop_reason_table
stop_reason_table_final <- bind_rows(
  stop_reason_table,
  
  tibble(
    `Reason for Discontinuation` = "Total",
    `n (% of discontinuations)` =
      sprintf(
        "%d (100.0%%)",
        n_discontinued
      ),
    
    `n (% of total cohort)` =
      sprintf(
        "%d (%.1f%%)",
        n_discontinued,
        100 * n_discontinued / n_total
      )
  )
)

stop_reason_table_final
write.csv(
  stop_reason_table_final,
  "Supplementary_Table_S1_Discontinuation_Reasons.csv",
  row.names = FALSE
)

# ============================================================
# Supplementary Table S2
# Maintenance Treatment Strategy Transitions Through Week 48
# ============================================================

n_total <- nrow(dat)

supp_table_s2 <- dat %>%
  count(
    strategy12,
    strategy24,
    strategy48,
    name = "n"
  ) %>%
  mutate(
    percentage = 100 * n / n_total
  ) %>%
  arrange(
    strategy12,
    strategy24,
    strategy48
  ) %>%
  transmute(
    `Weeks 12–24 Strategy` = as.character(strategy12),
    `Weeks 24–48 Strategy` = as.character(strategy24),
    `Beyond Week 48 Strategy` = as.character(strategy48),
    `Patients, n` = n,
    `Percentage of Total Cohort, %` = round(percentage, 1)
  )

# View table
supp_table_s2
write.csv(
  supp_table_s2,
  file = "Supplementary_Table_S2_Treatment_Strategy_Transitions.csv",
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

# ============================================================
# Treatment strategy distribution
# ============================================================

strategy_labels <- c(
  "0" = "UPA 15 mg",
  "1" = "UPA 30 mg",
  "2" = "UPA + advanced therapy",
  "3" = "Discontinued",
  "4" = "Lost/unknown"
)

# ------------------------------------------------------------
# 1. Distribution in the total cohort
# ------------------------------------------------------------

strategy_distribution_total <- bind_rows(
  
  dat %>%
    count(strategy24_plot, name = "n") %>%
    mutate(
      Period = "Weeks 24–48",
      percentage = 100 * n / nrow(dat)
    ) %>%
    rename(strategy = strategy24_plot),
  
  dat %>%
    count(strategy48_plot, name = "n") %>%
    mutate(
      Period = "Beyond Week 48",
      percentage = 100 * n / nrow(dat)
    ) %>%
    rename(strategy = strategy48_plot)
  
) %>%
  mutate(
    `Treatment Strategy` =
      recode(as.character(strategy), !!!strategy_labels),
    
    `n (%)` =
      sprintf("%d (%.1f%%)", n, percentage)
  ) %>%
  select(
    Period,
    `Treatment Strategy`,
    n,
    percentage,
    `n (%)`
  )

strategy_distribution_total
strategy_distribution_total_wide <- strategy_distribution_total %>%
  select(
    Period,
    `Treatment Strategy`,
    `n (%)`
  ) %>%
  tidyr::pivot_wider(
    names_from = Period,
    values_from = `n (%)`
  )

strategy_distribution_total_wide