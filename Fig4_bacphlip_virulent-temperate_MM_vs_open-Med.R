# =========================
# 0) Load Libraries
# =========================
library(tidyverse)
library(openxlsx)

# =========================
# 1) Paths
# =========================
input_file <- "Med-MM_vOTU_abund_table_BACPHLIP.csv"
output_dir <- "LYSOGENY"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
output_xlsx <- file.path(output_dir, "bacphlip_all_results.xlsx")

# =========================
# 2) Read & preprocess (DO NOT exclude shared vOTUs)
# =========================
df <- readr::read_delim(input_file, delim = ";", locale = readr::locale(decimal_mark = ","))

# Sample columns and groups
sample_cols <- colnames(df)[-(1:2)]
mar_menor_samples <- sample_cols[str_detect(sample_cols, "^[A-Z][a-z]{2}_[0-9]{4}")]
medit_sea_samples <- sample_cols[str_detect(sample_cols, "^TARA_|^Med")]

# Long format + groups
df_long <- df %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "Abundance_raw") %>%
  mutate(
    Abundance = as.numeric(str_replace(Abundance_raw, ",", ".")),
    Detected  = Abundance > 0,
    Group = case_when(
      Sample %in% mar_menor_samples ~ "Mar_Menor",
      Sample %in% medit_sea_samples  ~ "Medit_Sea",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(Group %in% c("Mar_Menor", "Medit_Sea"))

# Replicative strategy level order
if (!"replicative_strategy" %in% names(df_long)) {
  stop("Column 'replicative_strategy' not found in the input file.")
}
df_long <- df_long %>%
  mutate(replicative_strategy = factor(replicative_strategy,
                                       levels = c("Virulent", "Undetermined", "Temperate")))

# =========================
# 3) Sample × Strategy grid with ZEROS
# =========================
# Counts per sample and strategy (detected only)
counts_present <- df_long %>%
  filter(Detected) %>%
  group_by(Sample, Group, replicative_strategy) %>%
  summarise(n_vOTU = n_distinct(vOTU), .groups = "drop")

# Totals per sample (denominator)
totals_per_sample <- df_long %>%
  filter(Detected) %>%
  group_by(Sample, Group) %>%
  summarise(total_vOTUs_sample = n_distinct(vOTU), .groups = "drop")

# Complete Sample × Strategy grid
samples_df <- df_long %>% distinct(Sample, Group)
strategies  <- levels(df_long$replicative_strategy)

summary_df <- tidyr::crossing(
  samples_df,
  tibble(replicative_strategy = factor(strategies, levels = strategies))
) %>%
  left_join(counts_present, by = c("Sample", "Group", "replicative_strategy")) %>%
  left_join(totals_per_sample, by = c("Sample", "Group")) %>%
  mutate(
    n_vOTU = coalesce(n_vOTU, 0L),
    total_vOTUs_sample = coalesce(total_vOTUs_sample, 0L),
    perc = if_else(total_vOTUs_sample > 0, (n_vOTU / total_vOTUs_sample) * 100, 0)
  )

# Labels n= (unique vOTUs per group and strategy across the dataset)
labels_n <- df_long %>%
  filter(Detected) %>%
  distinct(Group, replicative_strategy, vOTU) %>%
  count(Group, replicative_strategy, name = "n_vOTUs") %>%
  mutate(label = paste0("n=", n_vOTUs))

# =========================
# 4) Boxplot (ALL vOTUs; ZEROS included)
# =========================
p <- ggplot(summary_df, aes(x = replicative_strategy, y = perc, fill = Group)) +
  geom_boxplot(outlier.shape = NA, coef = 1.5,
               position = position_dodge(width = 0.75)) +
  # solid circles (as in plot_box())
  geom_point(
    aes(color = Group, size = Group),
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
    alpha = 0.65,
    shape = 16
  ) +
  geom_text(
    data = labels_n,
    aes(x = replicative_strategy, y = 95, label = label, group = Group),
    position = position_dodge(width = 0.75),
    size = 3.5, vjust = 0, angle = 45
  ) +
  scale_size_manual(values = c("Mar_Menor" = 0.8, "Medit_Sea" = 0.8)) +
  scale_fill_manual(
    values = c("Mar_Menor" = "lightgreen", "Medit_Sea" = "skyblue"),
    labels = c("Mar_Menor" = "Mar Menor", "Medit_Sea" = "Medit Sea")
  ) +
  scale_color_manual(
    values = c("Mar_Menor" = "darkgreen", "Medit_Sea" = "deepskyblue4"),
    labels = c("Mar_Menor" = "Mar Menor", "Medit_Sea" = "Medit Sea")
  ) +
  scale_y_continuous(name = "Percent of hq-vOTUs (%)", limits = c(0, 100)) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    panel.background = element_blank()
  ) +
  labs(
    title = "Replicative strategy (ALL vOTUs; zeros included)",
    x = "Replicative strategy", fill = "Environment", color = "Environment"
  )

ggsave(file.path(output_dir, "BACPHLIP_boxplot_ALL_zeros_included.png"),
       p, width = 8, height = 6)
print(p)

ggsave(file.path(output_dir, "bacphlip_boxplot.svg"), p, width = 10, height = 10, units = "cm", dpi = 300)

# =========================
# 5) Wilcoxon per strategy (ZEROS included)
# =========================
wilcox_results <- summary_df %>%
  group_by(replicative_strategy) %>%
  summarise(
    test = list(tryCatch(
      stats::wilcox.test(perc ~ Group, data = cur_data(), exact = FALSE),
      error = function(e) NULL
    )),
    n_Mar_Menor = sum(cur_data()$Group == "Mar_Menor"),
    n_Medit_Sea = sum(cur_data()$Group == "Medit_Sea"),
    .groups = "drop"
  ) %>%
  filter(!sapply(test, is.null)) %>%
  mutate(
    p_value  = purrr::map_dbl(test, "p.value"),
    statistic = purrr::map_dbl(test, ~ unname(.$statistic)),
    p_adj = p.adjust(p_value, method = "BH"),
    signif = dplyr::case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01  ~ "**",
      p_adj <= 0.05  ~ "*",
      TRUE ~ "ns"
    )
  ) %>%
  select(-test) %>%
  arrange(p_adj)

# =========================
# 6) Save EVERYTHING to a single Excel
# =========================
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Perc_by_sample_ALL")
openxlsx::writeData(wb, "Perc_by_sample_ALL", summary_df)

summary_grouped <- summary_df %>%
  group_by(Group, replicative_strategy) %>%
  summarise(
    mean_perc   = mean(perc, na.rm = TRUE),
    sd_perc     = sd(perc, na.rm = TRUE),
    n_samples   = n_distinct(Sample),
    total_vOTUs = sum(n_vOTU),
    .groups = "drop"
  )

openxlsx::addWorksheet(wb, "Perc_grouped_ALL")
openxlsx::writeData(wb, "Perc_grouped_ALL", summary_grouped)

openxlsx::addWorksheet(wb, "Wilcoxon_ALL")
openxlsx::writeData(wb, "Wilcoxon_ALL", wilcox_results)

openxlsx::addWorksheet(wb, "Labels_n")
openxlsx::writeData(wb, "Labels_n", labels_n)

# Parameters/metadata sheet for traceability
params <- tibble::tibble(
  Item  = c(
    "Shared vOTUs excluded?",
    "Zeros included (Sample×Strategy)?",
    "Statistical test",
    "Multiple testing correction",
    "Input file",
    "Generated figure"
  ),
  Value = c(
    "No (ALL vOTUs included)",
    "Yes (complete grid with zeros)",
    "Wilcoxon rank-sum two-sided, exact=FALSE",
    "BH (FDR)",
    input_file,
    file.path(output_dir, "BACPHLIP_boxplot_ALL_zeros_included.png")
  )
)
openxlsx::addWorksheet(wb, "Params")
openxlsx::writeData(wb, "Params", params)

openxlsx::saveWorkbook(wb, output_xlsx, overwrite = TRUE)

cat("\nDone. Excel with all tabular outputs at:\n", output_xlsx, "\n")
