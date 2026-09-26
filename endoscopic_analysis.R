# ============================================================
# 0. 加载包
# ============================================================
library(ggplot2)

# 时间点配色
time_cols <- c(
  "Week 12" = "#F28E2B",
  "Week 24" = "#4C78A8",
  "Week 48" = "#59A14F"
)

# 统一主题
academic_theme <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      size = 15,
      face = "bold",
      hjust = 0.5,
      margin = margin(b = 12)
    ),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 13),
    axis.text.x = element_text(
      size = 12,
      colour = "black",
      margin = margin(t = 6)
    ),
    axis.text.y = element_text(
      size = 11,
      colour = "black"
    ),
    axis.line = element_line(
      colour = "black",
      linewidth = 0.6
    ),
    axis.ticks = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    panel.grid.major.y = element_line(
      colour = "#D9D9D9",
      linewidth = 0.4
    ),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 11),
    plot.margin = margin(12, 18, 10, 12)
  )


# ============================================================
# Figure 1. Clinical outcomes under mNRI
# ============================================================
clinical_data <- data.frame(
  Outcome = rep(
    c("Clinical response", "Clinical remission"),
    each = 3
  ),
  Time = rep(
    c("Week 12", "Week 24", "Week 48"),
    times = 2
  ),
  Percentage = c(
    90.7, 82.1, 60.9,  # Clinical response
    47.0, 62.9, 44.4   # Clinical remission
  )
)

clinical_data$Outcome <- factor(
  clinical_data$Outcome,
  levels = c("Clinical response", "Clinical remission")
)

clinical_data$Time <- factor(
  clinical_data$Time,
  levels = c("Week 12", "Week 24", "Week 48")
)

# ============================================================
# Figure 1
# ============================================================
pd1 <- position_dodge(width = 0.80)

p1 <- ggplot(
  clinical_data,
  aes(
    x = Outcome,
    y = Percentage,
    fill = Time,
    group = Time
  )
) +
  geom_col(
    position = pd1,
    width = 0.72,
    colour = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", Percentage)),
    position = pd1,
    vjust = -0.35,
    size = 4.2
  ) +
  scale_fill_manual(values = time_cols) +
  scale_x_discrete(
    expand = expansion(add = 0.55)
  ) +
  scale_y_continuous(
    limits = c(0, 105),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Clinical outcomes under modified non-responder imputation",
    y = "Patients (%)"
  ) +
  academic_theme

p1


# 保存 Figure 1
ggsave(
  filename = "Figure_1_Clinical_outcomes_mNRI.png",
  plot = p1,
  width = 7.2,
  height = 5.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "Figure_1_Clinical_outcomes_mNRI.pdf",
  plot = p1,
  width = 7.2,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)


# ============================================================
# Figure 2. Endoscopic outcomes
# ============================================================
endoscopic_data <- data.frame(
  Outcome = rep(
    c(
      "Endoscopic remission",
      "Endoscopic response",
      "Non-response"
    ),
    each = 2
  ),
  Time = rep(
    c("Week 24", "Week 48"),
    times = 3
  ),
  Percentage = c(
    28.1, 35.5,  # Endoscopic remission
    42.7, 32.9,  # Endoscopic response
    29.2, 31.6   # Non-response
  )
)

endoscopic_data$Outcome <- factor(
  endoscopic_data$Outcome,
  levels = c(
    "Endoscopic remission",
    "Endoscopic response",
    "Non-response"
  )
)

endoscopic_data$Time <- factor(
  endoscopic_data$Time,
  levels = c("Week 24", "Week 48")
)

# ============================================================
# Figure 2
# ============================================================
pd2 <- position_dodge(width = 0.80)

p2 <- ggplot(
  endoscopic_data,
  aes(
    x = Outcome,
    y = Percentage,
    fill = Time,
    group = Time
  )
) +
  geom_col(
    position = pd2,
    width = 0.72,
    colour = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", Percentage)),
    position = pd2,
    vjust = -0.35,
    size = 4.2
  ) +
  scale_fill_manual(values = time_cols) +
  scale_x_discrete(
    expand = expansion(add = 0.55)
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "Endoscopic outcomes",
    y = "Patients (%)"
  ) +
  academic_theme

p2


# 保存 Figure 2
ggsave(
  filename = "Figure_2_Endoscopic_outcomes.png",
  plot = p2,
  width = 7.5,
  height = 5.2,
  units = "in",
  dpi = 300,
  bg = "white"
)

ggsave(
  filename = "Figure_2_Endoscopic_outcomes.pdf",
  plot = p2,
  width = 7.5,
  height = 5.2,
  units = "in",
  device = cairo_pdf
)