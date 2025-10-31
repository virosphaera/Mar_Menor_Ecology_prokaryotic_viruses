# ─────────────────────────────────────────────────────────
# Load libraries
# ─────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(ggplot2)
library(rstatix)
library(forcats)
library(purrr)

# ─────────────────────────────────────────────────────────
# Parameters (adjustable)
# ─────────────────────────────────────────────────────────
P_THRESHOLD <- 0.05
EFFECT_THRESHOLD <- 0.1
MIN_VOTUS_PER_ECOSYS <- 5

# ─────────────────────────────────────────────────────────
# Define paths
# ─────────────────────────────────────────────────────────
input_dir <- "HOSTS"
abundance_file <- file.path(input_dir, "Med-MM_vOTU_abund_table.csv")
out_dir <- "HOSTS"

# ─────────────────────────────────────────────────────────
# Read and process data
# ─────────────────────────────────────────────────────────
abund <- read_csv2(abundance_file)

# Convert decimal commas to dots if needed
abund[,-c(1:2)] <- abund[,-c(1:2)] %>%
  mutate(across(everything(), ~ as.numeric(str_replace(., ",", "."))))

combined <- abund %>%
  pivot_longer(cols = -c(vOTU, Lowest_taxon),
               names_to = "Sample", values_to = "Abundance") %>%
  mutate(
    Present = ifelse(Abundance > 0, 1, 0),
    Virus_abundance_table_ID = vOTU,
    Sample_group = case_when(
      str_detect(Sample, "^TARA_") | str_detect(Sample, "^Med-") ~ "Medit_Sea",
      str_detect(Sample, "^[A-Z][a-z]{2}_\\d{4}$") ~ "Mar_Menor",
      TRUE ~ NA_character_
    )
  )

# ─────────────────────────────────────────────────────────
# % of vOTUs per host taxon per sample (all vOTUs)
# ─────────────────────────────────────────────────────────
votus_total_sample <- combined %>%
  filter(Abundance > 0) %>%
  group_by(Sample) %>%
  summarise(Total_vOTUs = n_distinct(Virus_abundance_table_ID), .groups = "drop")

combined_split <- combined %>%
  filter(Abundance > 0, Lowest_taxon != "No_host") %>%
  separate_rows(Lowest_taxon, sep = " OR ") %>%
  distinct(Sample, Lowest_taxon, Virus_abundance_table_ID, Sample_group, .keep_all = TRUE)

votus_per_sample_taxa <- combined_split %>%
  group_by(Sample, Sample_group, Lowest_taxon) %>%
  summarise(Num_vOTUs = n(), .groups = "drop") %>%
  left_join(votus_total_sample, by = "Sample") %>%
  mutate(Percent = Num_vOTUs / Total_vOTUs * 100)

# Total vOTUs per host (all vOTUs)
total_votus_by_host <- combined_split %>%
  distinct(Lowest_taxon, Virus_abundance_table_ID) %>%
  count(Lowest_taxon, name = "Num_total_vOTUs_host")

write_csv2(total_votus_by_host, file.path(dirname(abundance_file), "votus_total_per_host.csv"))
write_csv2(votus_per_sample_taxa, file.path(dirname(abundance_file), "votus_per_sample_taxa.csv"))

# ─────────────────────────────────────────────────────────
# Fill zeros for all Sample × Lowest_taxon combinations
#    (include 0% for effect size and Wilcoxon)
# ─────────────────────────────────────────────────────────
samples_df <- combined %>%
  filter(!is.na(Sample_group)) %>%
  distinct(Sample, Sample_group)

all_taxa <- combined_split %>%
  distinct(Lowest_taxon)

base_grid <- tidyr::crossing(samples_df, all_taxa)

counts_detected <- combined_split %>%
  group_by(Sample, Sample_group, Lowest_taxon) %>%
  summarise(Num_vOTUs = n(), .groups = "drop")

votus_per_sample_taxa_complete <- base_grid %>%
  left_join(counts_detected,
            by = c("Sample", "Sample_group", "Lowest_taxon")) %>%
  mutate(Num_vOTUs = coalesce(Num_vOTUs, 0L)) %>%
  left_join(votus_total_sample, by = "Sample") %>%
  mutate(Percent = if_else(Total_vOTUs > 0, Num_vOTUs / Total_vOTUs * 100, 0))

# ─────────────────────────────────────────────────────────
# Eligible taxa (≥ MIN_VOTUS_PER_ECOSYS in ≥1 ecosystem)  ← OR condition
# ─────────────────────────────────────────────────────────
taxa_counts_per_ecosys <- combined_split %>%
  group_by(Lowest_taxon, Sample_group) %>%
  summarise(n_vOTUs = n_distinct(Virus_abundance_table_ID), .groups = "drop")

eligible_taxa <- taxa_counts_per_ecosys %>%
  tidyr::pivot_wider(names_from = Sample_group, values_from = n_vOTUs) %>%
  mutate(
    Mar_Menor = coalesce(Mar_Menor, 0L),
    Medit_Sea = coalesce(Medit_Sea, 0L)
  ) %>%
  filter(Mar_Menor >= MIN_VOTUS_PER_ECOSYS | Medit_Sea >= MIN_VOTUS_PER_ECOSYS) %>%
  pull(Lowest_taxon)

# ─────────────────────────────────────────────────────────
# Wilcoxon test + Effect size (includes 0%)
# ─────────────────────────────────────────────────────────
base_tbl <- votus_per_sample_taxa_complete %>%
  filter(Lowest_taxon %in% eligible_taxa)

valid_taxa <- base_tbl %>%
  group_by(Lowest_taxon) %>%
  summarise(n_groups = n_distinct(Sample_group), .groups = "drop") %>%
  filter(n_groups == 2) %>%
  pull(Lowest_taxon)

medians_by_taxon <- base_tbl %>%
  filter(Lowest_taxon %in% valid_taxa) %>%
  group_by(Lowest_taxon, Sample_group) %>%
  summarise(median_percent = median(Percent), .groups = "drop") %>%
  pivot_wider(names_from = Sample_group, values_from = median_percent, names_prefix = "median_") %>%
  rename_with(~ str_replace_all(., " ", "_"))

wilcox_summary_all <- base_tbl %>%
  filter(Lowest_taxon %in% valid_taxa) %>%
  group_by(Lowest_taxon) %>%
  summarise(data = list(pick(everything())), .groups = "drop") %>%
  mutate(
    test = map(data, ~ {
      tmp <- .x
      counts <- table(tmp$Sample_group)
      if (all(c("Mar_Menor", "Medit_Sea") %in% names(counts)) && all(counts >= 2)) {
        rstatix::wilcox_test(tmp, Percent ~ Sample_group, exact = FALSE)
      } else {
        tibble(statistic = NA, p = NA)
      }
    })
  ) %>%
  unnest(test) %>%
  select(-data) %>%
  left_join(medians_by_taxon, by = "Lowest_taxon") %>%
  mutate(
    effect = round(median_Mar_Menor - median_Medit_Sea, 3),
    significance = case_when(
      p <= 0.001 ~ "***",
      p <= 0.01  ~ "**",
      p <= 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# Save ALL Wilcoxon results (no final filtering)
wilcox_all_out <- wilcox_summary_all %>%
  left_join(
    taxa_counts_per_ecosys %>%
      tidyr::pivot_wider(names_from = Sample_group, values_from = n_vOTUs,
                         values_fill = 0, names_prefix = "n_"),
    by = "Lowest_taxon"
  ) %>%
  rename(N_Mar_Menor = n_Mar_Menor, N_Medit_Sea = n_Medit_Sea) %>%
  arrange(p)

write_csv2(wilcox_all_out, file.path(dirname(abundance_file), "wilcoxon_results_all.csv"))

# Final filter (p ≤ 0.05 and |effect| ≥ 0.1)
wilcox_summary <- wilcox_summary_all %>%
  filter(!is.na(p), p <= P_THRESHOLD, abs(effect) >= EFFECT_THRESHOLD)

# Save FILTERED Wilcoxon results (effect size + statistic)
wilcox_filtered_out <- wilcox_summary %>%
  left_join(
    taxa_counts_per_ecosys %>%
      tidyr::pivot_wider(names_from = Sample_group, values_from = n_vOTUs,
                         values_fill = 0, names_prefix = "n_"),
    by = "Lowest_taxon"
  ) %>%
  rename(N_Mar_Menor = n_Mar_Menor, N_Medit_Sea = n_Medit_Sea) %>%
  select(Lowest_taxon, N_Mar_Menor, N_Medit_Sea,
         median_Mar_Menor, median_Medit_Sea,
         effect, p, significance, statistic) %>%
  arrange(p)

write_csv2(wilcox_filtered_out, file.path(dirname(abundance_file), "wilcoxon_results_filtered.csv"))

# ─────────────────────────────────────────────────────────
# Effect size plot
# ─────────────────────────────────────────────────────────
plot_df <- wilcox_summary %>%
  mutate(
    Short_taxon = str_replace(Lowest_taxon, "^[a-z]__", ""),
    Short_taxon = fct_reorder(Short_taxon, effect)
  )

if (nrow(plot_df) > 0) {
  effect_size_plot <- ggplot(plot_df, aes(x = effect, y = Short_taxon, fill = effect)) +
    geom_col(width = 1, color = "black", size = 0.2) +
    geom_text(aes(label = significance),
              hjust = ifelse(plot_df$effect >= 0, -0.2, 1.2), size = 5) +
    scale_fill_gradient2(
      low = "#9ecae1", mid = "white", high = "#a1d99b",
      midpoint = 0,
      limits = c(min(plot_df$effect), max(plot_df$effect)),
      name = "Effect size"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10, color = "black"),
      axis.text.x = element_text(size = 10, color = "black"),
      axis.title = element_text(size = 12),
      plot.title = element_text(size = 14, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, size = 0.3)
    ) +
    labs(
      x = "Effect size (Mar Menor - Medit Sea)",
      y = "Host taxon",
      title = paste0("Taxa (all vOTUs; p ≤ ", P_THRESHOLD,
                     " and |effect| ≥ ", EFFECT_THRESHOLD,
                     "; ≥", MIN_VOTUS_PER_ECOSYS, " vOTUs in ≥1 ecosystem; including 0%)")
    )
  print(p)
} else {
  message("No taxa met the criteria for the plot.")
}

ggsave(file.path(out_dir, "Effect_size_plot_at_least_5_vOTUs_in_at_least_1_ecosystem_and_including_those_w_0percent_as_well.svg"),
       effect_size_plot, width = 22, height = 14, units = "cm", dpi = 300)

# ─────────────────────────────────────────────────────────
# Boxplot + save input (counts and percents)
# ─────────────────────────────────────────────────────────
sig_taxa <- wilcox_summary %>% pull(Lowest_taxon)

if (length(sig_taxa) > 0) {
  votu_counts <- combined_split %>%
    filter(Lowest_taxon %in% sig_taxa) %>%
    group_by(Lowest_taxon) %>%
    summarise(vOTU_count = n_distinct(Virus_abundance_table_ID), .groups = "drop")
  
  box_df <- base_tbl %>%
    filter(Lowest_taxon %in% sig_taxa) %>%
    left_join(votu_counts, by = "Lowest_taxon") %>%
    mutate(
      Short_taxon = str_replace(Lowest_taxon, "^[a-z]__", ""),
      Short_taxon = factor(
        Short_taxon,
        levels = votu_counts %>%
          mutate(Short_taxon = str_replace(Lowest_taxon, "^[a-z]__", "")) %>%
          arrange(desc(vOTU_count)) %>%
          pull(Short_taxon)
      )
    )
  
  # Save boxplot input (counts + percents per sample/taxon)
  boxplot_out <- box_df %>%
    select(Sample, Sample_group, Lowest_taxon, Short_taxon,
           Num_vOTUs, Total_vOTUs, Percent, vOTU_count)
  write_csv2(boxplot_out, file.path(dirname(abundance_file), "boxplot_input_counts_percents.csv"))
  
  p_box <- ggplot(box_df, aes(x = Short_taxon, y = Percent, fill = Sample_group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    geom_jitter(
      aes(color = Sample_group),
      position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75),
      size = 0.8, alpha = 0.6, show.legend = FALSE
    ) +
    scale_fill_manual(
      values = c("Mar_Menor" = "lightgreen", "Medit_Sea" = "skyblue"),
      labels = c("Mar_Menor" = "Mar Menor", "Medit_Sea" = "Medit Sea"),
      name = NULL
    ) +
    scale_color_manual(
      values = c("Mar_Menor" = "darkgreen", "Medit_Sea" = "steelblue"),
      labels = c("Mar_Menor" = "Mar Menor", "Medit_Sea" = "Medit Sea"),
      name = NULL
    ) +
    labs(
      title = paste0("% of vOTUs by host taxon (all vOTUs; p ≤ ",
                     P_THRESHOLD, " and |effect| ≥ ", EFFECT_THRESHOLD,
                     "; ≥", MIN_VOTUS_PER_ECOSYS, " vOTUs in ≥1 ecosystem; including 0%)"),
      y = "Percent of vOTUs (%)",
      x = "Host taxon"
    ) +
    theme_minimal(base_size = 15) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),
      axis.text.y = element_text(hjust = 1, color = "black"),
      legend.position = "top",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3)
    )
  
  print(p_box)
  
} else {
  message("No taxa met the criteria for the boxplot.")
}

ggsave(file.path(out_dir, "Box_plot_at_least_5_vOTUs_in_at_least_1_ecosystem_and_including_those_w_0percent_as_well.svg"),
       p_box, width = 19, height = 14, units = "cm", dpi = 300)

# ─────────────────────────────────────────────────────────
# Save both plot datasets + all-host data in one Excel file
# ─────────────────────────────────────────────────────────
library(writexl)

effect_export <- wilcox_summary %>%
  left_join(
    taxa_counts_per_ecosys %>%
      tidyr::pivot_wider(names_from = Sample_group, values_from = n_vOTUs,
                         values_fill = 0, names_prefix = "n_"),
    by = "Lowest_taxon"
  ) %>%
  rename(N_Mar_Menor = n_Mar_Menor, N_Medit_Sea = n_Medit_Sea) %>%
  mutate(Short_taxon = str_replace(Lowest_taxon, "^[a-z]__", "")) %>%
  select(Lowest_taxon, Short_taxon,
         N_Mar_Menor, N_Medit_Sea,
         median_Mar_Menor, median_Medit_Sea,
         effect, p, significance, statistic) %>%
  arrange(p)

if (exists("box_df")) {
  boxplot_export <- box_df %>%
    select(Sample, Sample_group, Lowest_taxon, Short_taxon,
           Num_vOTUs, Total_vOTUs, Percent, vOTU_count)
} else {
  boxplot_export <- tibble(
    Sample = character(), Sample_group = character(),
    Lowest_taxon = character(), Short_taxon = character(),
    Num_vOTUs = numeric(), Total_vOTUs = numeric(),
    Percent = numeric(), vOTU_count = numeric()
  )
}

all_hosts_export <- votus_per_sample_taxa_complete %>%
  select(Sample, Sample_group, Lowest_taxon, Num_vOTUs, Total_vOTUs, Percent) %>%
  arrange(Sample_group, Sample, desc(Percent), Lowest_taxon)

params_export <- tibble::tibble(
  Parameter = c("P_THRESHOLD", "EFFECT_THRESHOLD", "MIN_VOTUS_PER_ECOSYS",
                "Effect_definition", "Data_source", "Zeros_included"),
  Value     = c(P_THRESHOLD, EFFECT_THRESHOLD, MIN_VOTUS_PER_ECOSYS,
                "median(Mar_Menor) - median(Medit_Sea)",
                abundance_file, "Yes (for Wilcoxon & effect)")
)

output_dir <- "C:/Users/guillermo.dominguez/Documents/EPRA/VIRUS_MAR_MENOR/MM_prok_viruses_dynamics_and_comparison_w_Medit/OUTPUTS/H
