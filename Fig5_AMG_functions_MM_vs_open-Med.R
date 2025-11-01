# ==== Load libraries ====
library(tidyverse)
library(rstatix)
library(writexl)
library(glue)
library(forcats)

# ==== Set paths ====
data_dir   <- "DATA/AMGs"
output_dir <- "OUTPUTS/AMGs"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ==== Parameters ====
EFFECT_THRESHOLD_PP <- 0.1   # Minimum |Δ median| in percentage points

# ==== Load data ====
df_main <- readr::read_delim(
  file.path(data_dir, "Med-MM_vOTU_abund_table_con_funciones.csv"),
  delim = ";", locale = readr::locale(decimal_mark = ",")
)

# ==== Identify metadata and sample columns ====
metadata_cols <- c("vOTU", "kegg_hit", "FUNCTION_ID", "AMG_count")
sample_cols   <- setdiff(names(df_main), metadata_cols)

# ==== Reshape to long format and filter detected vOTUs ====
df_long <- df_main %>%
  pivot_longer(cols = all_of(sample_cols), names_to = "Sample", values_to = "Abundance") %>%
  filter(Abundance > 0) %>%
  mutate(FUNCTION_ID = as.character(FUNCTION_ID))

# ==== Classify samples by ecosystem ====
df_long <- df_long %>%
  mutate(Group = ifelse(grepl("^[A-Z][a-z]{2}_\\d{4}$", Sample), "Mar_Menor", "Medit_Sea"))

# ==== Keep only vOTUs exclusive to one ecosystem ====
vOTU_group_presence <- df_long %>%
  distinct(vOTU, Group) %>%
  count(vOTU, name = "n_groups") %>%
  filter(n_groups == 1) %>%
  pull(vOTU)

df_long <- df_long %>% filter(vOTU %in% vOTU_group_presence)

# ==== % of vOTUs with any AMG per sample (global) ====
df_amg_fraction <- df_long %>%
  mutate(has_amg = !is.na(FUNCTION_ID) & FUNCTION_ID != "") %>%
  group_by(Sample, Group) %>%
  summarise(
    total_vOTUs      = n_distinct(vOTU),
    vOTUs_with_AMG   = n_distinct(vOTU[has_amg]),
    percent_with_AMG = (vOTUs_with_AMG / total_vOTUs) * 100,
    .groups = "drop"
  )

# ==== Expand FUNCTION_IDs (split by " OR ") & drop non-metabolic ====
df_amg_expanded <- df_long %>%
  separate_rows(FUNCTION_ID, sep = " OR ") %>%
  filter(!is.na(FUNCTION_ID), FUNCTION_ID != "", FUNCTION_ID != "Non_Metabolic_Function")

# ==== Totals per sample (denominator) ====
total_detected_vOTUs <- df_long %>%
  group_by(Sample) %>%
  summarise(total_vOTUs_detected = n_distinct(vOTU), .groups = "drop")

# ==== Count vOTUs per function and sample (present) ====
vOTUs_per_function <- df_amg_expanded %>%
  group_by(Sample, FUNCTION_ID) %>%
  summarise(n_vOTUs = n_distinct(vOTU), .groups = "drop")

# ==== Eligible functions: ≥5 unique vOTUs in ≥1 ecosystem (OR) ====
functions_with_min_vOTUs <- df_amg_expanded %>%
  distinct(vOTU, FUNCTION_ID, Group) %>%
  count(FUNCTION_ID, Group, name = "n_vOTUs") %>%
  pivot_wider(names_from = Group, values_from = n_vOTUs, values_fill = 0) %>%
  filter(Mar_Menor >= 5 | Medit_Sea >= 5) %>%
  pull(FUNCTION_ID)

# ==== Build complete grid with zeros for analysis (eligible funcs only) ====
samples_df <- total_detected_vOTUs %>%
  mutate(Group = ifelse(grepl("^[A-Z][a-z]{2}_\\d{4}$", Sample), "Mar_Menor", "Medit_Sea"))

base_grid <- tidyr::crossing(
  Sample = samples_df$Sample,
  FUNCTION_ID = functions_with_min_vOTUs
)

vOTUs_percent_filtered <- base_grid %>%
  left_join(vOTUs_per_function, by = c("Sample", "FUNCTION_ID")) %>%
  mutate(n_vOTUs = coalesce(n_vOTUs, 0L)) %>%
  left_join(total_detected_vOTUs, by = "Sample") %>%
  mutate(percent_vOTUs = (n_vOTUs / total_vOTUs_detected) * 100) %>%
  left_join(samples_df %>% select(Sample, Group), by = "Sample")

# ==== Median differences (pp) for all eligible functions ====
medians_diff_all <- vOTUs_percent_filtered %>%
  group_by(FUNCTION_ID, Group) %>%
  summarise(median_pp = median(percent_vOTUs, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = median_pp,
              names_prefix = "median_", values_fill = 0) %>%
  mutate(effect_pp = round(median_Mar_Menor - median_Medit_Sea, 3))

# ==== Wilcoxon (two-sided, exact=FALSE) + BH ====
wilcox_results <- vOTUs_percent_filtered %>%
  group_by(FUNCTION_ID) %>%
  summarise(
    test = list(tryCatch(
      stats::wilcox.test(percent_vOTUs ~ Group, data = cur_data(), exact = FALSE),
      error = function(e) NULL
    )),
    .groups = "drop"
  ) %>%
  filter(!sapply(test, is.null)) %>%
  mutate(
    p_value   = purrr::map_dbl(test, "p.value"),
    statistic = purrr::map_dbl(test, ~ unname(.$statistic))
  ) %>%
  select(-test) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))

# ==== Significant functions: BH p_adj < 0.05 & |effect_pp| ≥ threshold ====
sig_tbl <- wilcox_results %>%
  inner_join(medians_diff_all, by = "FUNCTION_ID") %>%
  filter(p_adj < 0.05, abs(effect_pp) >= EFFECT_THRESHOLD_PP) %>%
  arrange(p_adj)

sig_functions <- sig_tbl$FUNCTION_ID

# ==== # of exclusive vOTUs per group/function (for labeling) ====
n_vOTUs_per_group <- df_amg_expanded %>%
  filter(FUNCTION_ID %in% sig_functions) %>%
  distinct(vOTU, FUNCTION_ID, Group) %>%
  count(FUNCTION_ID, Group, name = "n_vOTUs") %>%
  pivot_wider(names_from = Group, values_from = n_vOTUs, values_fill = 0) %>%
  rename(`# vOTUs Mar_Menor` = Mar_Menor,
         `# vOTUs Medit_Sea` = Medit_Sea)

# ==== Boxplots (only if significant functions exist) ====
if (length(sig_functions) > 0) {
  p3 <- ggplot(
    vOTUs_percent_filtered %>% filter(FUNCTION_ID %in% sig_functions),
    aes(x = Group, y = percent_vOTUs, fill = Group)
  ) +
    geom_boxplot(outlier.shape = NA, alpha = 0.7, color = "black") +
    geom_point(
      aes(color = Group, size = Group),
      position = position_jitter(width = 0.2),
      alpha = 0.65,
      shape = 16
    ) +
    facet_wrap(~ FUNCTION_ID, scales = "free_y", ncol = 2) +
    scale_size_manual(values = c("Mar_Menor" = 0.8, "Medit_Sea" = 0.8)) +
    scale_fill_manual(values = c("Mar_Menor" = "lightgreen", "Medit_Sea" = "lightblue")) +
    scale_color_manual(values = c("Mar_Menor" = "darkgreen", "Medit_Sea" = "steelblue")) +
    labs(
      title = "Significant viral functions (≥5 exclusive vOTUs in ≥1 ecosystem; zeros included)",
      y = "Percent of vOTUs (%)", x = ""
    ) +
    theme_minimal(base_size = 12) +
    guides(color = "none")
  
  ggsave(file.path(output_dir, "boxplot_significant_functions_exclusive_filtered.png"),
         p3, width = 12, height = 8)
  print(p3)
} else {
  message("No significant viral functions after BH p_adj < 0.05 & |Δ| ≥ 0.1 pp.")
}

ggsave("boxplots_AMG_categories_w_significant_differences.svg", p3, width=14, height=20, units="cm", dpi=300)

# ==== Median-difference barplot (pp) — with zeros ====
if (length(sig_functions) > 0) {
  medians_diff_sig <- medians_diff_all %>%
    filter(FUNCTION_ID %in% sig_functions) %>%
    left_join(sig_tbl %>% select(FUNCTION_ID, p_adj), by = "FUNCTION_ID") %>%
    left_join(n_vOTUs_per_group, by = "FUNCTION_ID") %>%
    mutate(
      significance = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE ~ ""
      ),
      Label_txt = str_replace_all(FUNCTION_ID, "_", " "),
      Label_txt = glue("{Label_txt} ({`# vOTUs Medit_Sea`}/{`# vOTUs Mar_Menor`})"),
      Label     = fct_reorder(Label_txt, effect_pp),
      hjust_pos = if_else(effect_pp >= 0, -0.2, 1.2)
    )
  
  p_median <- ggplot(medians_diff_sig, aes(x = effect_pp, y = Label, fill = effect_pp)) +
    geom_col(width = 1, color = "black", linewidth = 0.3) +
    geom_text(aes(label = significance, hjust = hjust_pos), size = 5) +
    scale_fill_gradient2(low = "#9ecae1", mid = "white", high = "#a1d99b",
                         midpoint = 0, limits = range(medians_diff_sig$effect_pp),
                         name = "Δ median (pp)") +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "Δ median (pp) · Mar_Menor − Medit_Sea",
         y = "AMG functional category",
         title = "Significant viral functions (median difference; zeros included)") +
    theme_minimal(base_size = 13) +
    theme(
      axis.text.y = element_text(size = 12, color ="black"),
      axis.text.x = element_text(size = 12, color="black"),
      axis.title  = element_text(size = 12),
      plot.title  = element_text(size = 12),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
    )
  
  ggsave(file.path(output_dir, "Size_effect_plot_AMGs_Med-MM.svg"),
         p_median, width = 23, height = 8, units="cm", dpi=300)
  print(p_median)
}

# ==== Global tests (sample-level) ====
test_amg_fraction <- stats::wilcox.test(percent_with_AMG ~ Group, data = df_amg_fraction)
cat("Wilcoxon test (percent vOTUs with AMG, exclusive): p =", test_amg_fraction$p.value, "\n")
print(rstatix::wilcox_effsize(df_amg_fraction, percent_with_AMG ~ Group))

# ==== Tab for ALL functions (complete grid with zeros) ====
all_functions <- df_amg_expanded %>% distinct(FUNCTION_ID) %>% pull()

base_grid_all <- tidyr::crossing(
  Sample = samples_df$Sample,
  FUNCTION_ID = all_functions
)

all_functions_complete <- base_grid_all %>%
  left_join(vOTUs_per_function, by = c("Sample", "FUNCTION_ID")) %>%
  mutate(n_vOTUs = coalesce(n_vOTUs, 0L)) %>%
  left_join(total_detected_vOTUs, by = "Sample") %>%
  mutate(percent_vOTUs = (n_vOTUs / total_vOTUs_detected) * 100) %>%
  left_join(samples_df %>% select(Sample, Group), by = "Sample") %>%
  arrange(Group, Sample, desc(percent_vOTUs), FUNCTION_ID)

# ==== Export tables ====
wilcox_export_all <- medians_diff_all %>%
  left_join(wilcox_results, by = "FUNCTION_ID") %>%
  left_join(
    df_amg_expanded %>%
      distinct(vOTU, FUNCTION_ID, Group) %>%
      count(FUNCTION_ID, Group, name = "n_vOTUs") %>%
      pivot_wider(names_from = Group, values_from = n_vOTUs, values_fill = 0) %>%
      rename(`# vOTUs Mar_Menor` = Mar_Menor, `# vOTUs Medit_Sea` = Medit_Sea),
    by = "FUNCTION_ID"
  ) %>%
  select(FUNCTION_ID,
         median_Mar_Menor, median_Medit_Sea, effect_pp,
         p_value, p_adj, statistic,
         `# vOTUs Mar_Menor`, `# vOTUs Medit_Sea`) %>%
  arrange(p_adj)

wilcox_export_filtered <- wilcox_export_all %>%
  filter(p_adj < 0.05, abs(effect_pp) >= EFFECT_THRESHOLD_PP)

# ==== Parameters sheet ====
sheet_params <- tibble::tibble(
  Item = c(
    "Input data filter",
    "Function eligibility rule",
    "Univariate test",
    "Multiple testing correction",
    "Zeros included in test/effect",
    "Threshold |Δ median| (pp)",
    "Exported statistics",
    "I/O paths"
  ),
  Value = c(
    "Only vOTUs exclusive to one ecosystem (n_groups == 1)",
    "≥5 vOTUs in ≥1 ecosystem (OR)",
    "Wilcoxon rank-sum (two-sided, exact=FALSE) on percent_vOTUs",
    "BH (FDR) correction on per-function p-values",
    "Yes (complete grid with zeros)",
    EFFECT_THRESHOLD_PP,
    "Wilcoxon W, p-value, BH p_adj; Δ median (pp); #vOTUs per group",
    paste0("input: ", file.path(data_dir, "Med-MM_vOTU_abund_table_con_funciones.csv"),
           " | output: ", output_dir)
  )
)

# ==== Assemble workbook & write ====
sheets <- list(
  p1_percent_vOTUs_with_AMG     = df_amg_fraction,
  p3_sig_functions_input        = if (length(sig_functions) > 0)
    vOTUs_percent_filtered %>% filter(FUNCTION_ID %in% sig_functions) else NULL,
  p7_effect_median_diff_summary = if (length(sig_functions) > 0) medians_diff_sig else NULL,
  All_functions_counts_percents = all_functions_complete %>%
    select(Sample, Group, FUNCTION_ID, n_vOTUs, total_vOTUs_detected, percent_vOTUs),
  Wilcoxon_results_all          = wilcox_export_all,
  Wilcoxon_results_filtered     = wilcox_export_filtered,
  Params                        = sheet_params
) %>% purrr::compact()

out_xlsx_all <- file.path(output_dir, "AMGs_ALL_RESULTS.xlsx")

safe_write_xlsx <- function(sheets_list, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (file.exists(path)) {
    removed <- tryCatch({ file.remove(path) }, warning = function(e) FALSE, error = function(e) FALSE)
    if (!removed) {
      path <- file.path(
        dirname(path),
        paste0(tools::file_path_sans_ext(basename(path)),
               "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".xlsx")
      )
      message("Destination locked; using alternate filename: ", path)
    }
  }
  writexl::write_xlsx(sheets_list, path)
  message("Excel created: ", path)
}

safe_write_xlsx(sheets, out_xlsx_all)
