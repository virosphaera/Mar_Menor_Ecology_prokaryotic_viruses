# =========================
# 0) Load libraries
# =========================
library(tidyverse)
library(openxlsx)
library(forcats)
library(glue)

# =========================
# 1) Define paths
# =========================
input_file <- "C:/Users/guillermo.dominguez/Documents/EPRA/VIRUS_MAR_MENOR/MM_prok_viruses_dynamics_and_comparison_w_Medit/DATA/LYSOGENY/Med-MM_vOTU_abund_table_lysogeny_markers.csv"
output_dir <- "C:/Users/guillermo.dominguez/Documents/EPRA/VIRUS_MAR_MENOR/MM_prok_viruses_dynamics_and_comparison_w_Medit/OUTPUTS/LYSOGENY/"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

EFFECT_THRESHOLD_PP <- 0.1  # minimum |Δ median| (pp)

# =========================
# 2) Read and preprocess
# =========================
df <- read_delim(input_file, delim = ";", locale = locale(decimal_mark = ",")) %>%
  separate_rows(Lowest_host_taxon, sep = " OR ") %>%
  filter(Lowest_host_taxon != "No_host")

metadata_cols <- c("vOTU", "Lowest_host_taxon", "Lysogeny_markers")
sample_cols <- setdiff(colnames(df), metadata_cols)

get_group <- function(x) {
  case_when(
    grepl("^(Oct|Feb|Jul|Aug|Dec|Apr|Nov|May)_\\d{4}$", x) ~ "Mar_Menor",
    grepl("^(TARA|Med)", x) ~ "Medit_Sea",
    TRUE ~ NA_character_
  )
}

# =========================
# 3) Exclude shared vOTUs (keep exclusives)
# =========================
votu_group_lookup <- df %>%
  select(vOTU, all_of(sample_cols)) %>%
  pivot_longer(-vOTU, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Abundance = as.numeric(str_replace(Abundance, ",", "."))) %>%
  filter(Abundance > 0) %>%
  mutate(Group = get_group(Sample)) %>%
  filter(!is.na(Group)) %>%
  distinct(vOTU, Group)

shared_votus <- votu_group_lookup %>%
  count(vOTU) %>%
  filter(n == 2) %>%
  pull(vOTU)

df_exclusive <- df %>% filter(!vOTU %in% shared_votus)

# =========================
# 4) Global comparison (per sample) — always plots
# =========================
df_summary <- df_exclusive %>%
  select(vOTU, Lysogeny_markers, all_of(sample_cols)) %>%
  pivot_longer(all_of(sample_cols), names_to = "Sample", values_to = "Abundance") %>%
  mutate(Abundance = as.numeric(str_replace(Abundance, ",", "."))) %>%
  filter(Abundance > 0) %>%
  group_by(Sample) %>%
  summarise(
    Total_vOTUs = n_distinct(vOTU),
    Marker_vOTUs = n_distinct(vOTU[Lysogeny_markers == "Yes"]),
    Marker_percent = (Marker_vOTUs / Total_vOTUs) * 100,
    .groups = "drop"
  ) %>%
  mutate(Group = get_group(Sample)) %>%
  filter(Group %in% c("Mar_Menor", "Medit_Sea"))

global_wilcox <- stats::wilcox.test(Marker_percent ~ Group, data = df_summary)

global_wilcox_df <- tibble(
  test = "Wilcoxon_global_EXCLUSIVE",
  group1 = "Mar_Menor",
  group2 = "Medit_Sea",
  n1 = sum(df_summary$Group == "Mar_Menor"),
  n2 = sum(df_summary$Group == "Medit_Sea"),
  statistic = unname(global_wilcox$statistic),
  p_value = global_wilcox$p.value
) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "BH"),
    signif = case_when(
      p_adj <= 0.001 ~ "***",
      p_adj <= 0.01 ~ "**",
      p_adj <= 0.05 ~ "*",
      TRUE ~ "ns"
    )
  )

gg_global <- ggplot(df_summary, aes(x = Group, y = Marker_percent, fill = Group)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_point(
    aes(color = Group, size = Group),
    position = position_jitter(width = 0.2),
    alpha = 0.65,
    shape = 16
  ) +
  scale_size_manual(values = c("Mar_Menor" = 0.8, "Medit_Sea" = 0.8)) +
  scale_fill_manual(values = c("Mar_Menor" = "lightgreen", "Medit_Sea" = "skyblue"),
                    labels = c("Mar_Menor" = "Mar Menor", "Medit_Sea" = "Medit Sea")) +
  scale_color_manual(values = c("Mar_Menor" = "darkgreen", "Medit_Sea" = "steelblue"),
                     labels = c("Mar_Menor" = "Mar Menor", "Medit_Sea" = "Medit Sea")) +
  labs(y = "Percent of temperate vOTUs (%)") +
  theme_minimal(base_size = 12) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black")
  )

ggsave(".../boxplot_global_lysogeny_marker.svg",
       gg_global, width = 8, height = 7.8, units ="cm", dpi = 300)
print(gg_global)

# =========================
# 5) Host-level analysis with filters + zeros
# =========================
long_excl <- df_exclusive %>%
  pivot_longer(all_of(sample_cols), names_to = "Sample", values_to = "Abundance") %>%
  mutate(
    Abundance = as.numeric(str_replace(Abundance, ",", ".")),
    Group = get_group(Sample)
  ) %>%
  filter(Group %in% c("Mar_Menor", "Medit_Sea"))

# ---- Eligibility (≥5 unique vOTUs in ≥1 ecosystem) ----
eligible_hosts <- long_excl %>%
  filter(Abundance > 0) %>%
  distinct(Lowest_host_taxon, vOTU, Group) %>%
  count(Lowest_host_taxon, Group, name = "n_vOTUs") %>%
  pivot_wider(names_from = Group, values_from = n_vOTUs, values_fill = 0) %>%
  filter(Mar_Menor >= 5 | Medit_Sea >= 5) %>%
  pull(Lowest_host_taxon)

if (length(eligible_hosts) == 0) {
  message("No eligible hosts (≥5 vOTUs in ≥1 ecosystem). Continuing with global and export only.")
}

# ---- Counts per Sample × Host (present) ----
counts_present <- long_excl %>%
  filter(Abundance > 0) %>%
  group_by(Sample, Group, Lowest_host_taxon) %>%
  summarise(
    Total_vOTUs = n_distinct(vOTU),
    Temperate_vOTUs = n_distinct(vOTU[Lysogeny_markers == 'Yes']),
    .groups = "drop"
  )

# ---- Complete grid with zeros (Sample × eligible Host) ----
samples_df <- long_excl %>% distinct(Sample, Group)

df_host_complete <- tidyr::crossing(
  samples_df,
  tibble(Lowest_host_taxon = eligible_hosts)
) %>%
  left_join(counts_present, by = c("Sample", "Group", "Lowest_host_taxon")) %>%
  mutate(
    Total_vOTUs = coalesce(Total_vOTUs, 0L),
    Temperate_vOTUs = coalesce(Temperate_vOTUs, 0L),
    Temperate_percent = if_else(Total_vOTUs > 0, Temperate_vOTUs / Total_vOTUs * 100, 0)
  )

# ---- Medians and effect (pp) per host ----
medians_diff_all <- df_host_complete %>%
  group_by(Lowest_host_taxon, Group) %>%
  summarise(median_pp = median(Temperate_percent, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Group, values_from = median_pp,
              names_prefix = "median_", values_fill = 0) %>%
  mutate(effect_pp = round(median_Mar_Menor - median_Medit_Sea, 3))

# ---- Wilcoxon per host (includes 0%) + BH ----
wilcox_all <- df_host_complete %>%
  group_by(Lowest_host_taxon) %>%
  summarise(
    test = list(tryCatch(
      stats::wilcox.test(Temperate_percent ~ Group, data = cur_data(), exact = FALSE),
      error = function(e) NULL
    )),
    .groups = "drop"
  ) %>%
  filter(!sapply(test, is.null)) %>%
  mutate(
    p_value  = purrr::map_dbl(test, "p.value"),
    statistic = purrr::map_dbl(test, ~ unname(.$statistic))
  ) %>%
  select(-test) %>%
  left_join(medians_diff_all, by = "Lowest_host_taxon") %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  arrange(p_adj)

# ---- Final filtering: BH p_adj < 0.05 & |Δ| ≥ 0.1 pp ----
wilcox_filtered <- wilcox_all %>%
  filter(!is.na(p_adj), p_adj < 0.05, abs(effect_pp) >= EFFECT_THRESHOLD_PP)

sig_hosts <- wilcox_filtered$Lowest_host_taxon

# =========================
# 6) Host-level plots
# =========================
# --- A) Significant hosts (if any) ---
if (length(sig_hosts) > 0) {
  gg_sig_box <- ggplot(
    df_host_complete %>% filter(Lowest_host_taxon %in% sig_hosts),
    aes(x = Group, y = Temperate_percent, fill = Group)
  ) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, shape = 21, size = 1.5, alpha = 0.7, color = "black") +
    facet_wrap(~ Lowest_host_taxon, scales = "free_y") +
    labs(
      title = "Temperate vOTUs by host (EXCLUSIVE; zeros included; BH p_adj<0.05 & |Δ|≥0.1 pp)",
      x = "Group", y = "Percent temperate vOTUs (%)"
    ) +
    scale_fill_manual(values = c("Mar_Menor" = "lightgreen", "Medit_Sea" = "skyblue")) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
      legend.position = "top"
    )
  ggsave(file.path(output_dir, "faceted_boxplot_temperate_by_host_EXCLUSIVE_SIGNIFICANT_filtered.png"),
         gg_sig_box, width = 12, height = 8)
  print(gg_sig_box)
  
  medians_sig <- medians_diff_all %>%
    filter(Lowest_host_taxon %in% sig_hosts) %>%
    left_join(wilcox_filtered %>% select(Lowest_host_taxon, p_adj), by = "Lowest_host_taxon") %>%
    mutate(
      significance = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE ~ ""
      ),
      Label = fct_reorder(Lowest_host_taxon, effect_pp),
      hjust_pos = if_else(effect_pp >= 0, -0.2, 1.2)
    )
  
  p_eff <- ggplot(medians_sig, aes(x = effect_pp, y = Label, fill = effect_pp)) +
    geom_col(width = 1, color = "black", linewidth = 0.3) +
    geom_text(aes(label = significance, hjust = hjust_pos), size = 5) +
    scale_fill_gradient2(low = "#9ecae1", mid = "white", high = "#a1d99b",
                         midpoint = 0, limits = range(medians_sig$effect_pp),
                         name = "Δ median (pp)") +
    labs(x = "Δ median (pp) · Mar_Menor − Medit_Sea",
         y = "Host taxon",
         title = "Median difference in temperate % by host (zeros included)") +
    theme_minimal(base_size = 12) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3))
  ggsave(file.path(output_dir, "effect_median_diff_barplot_temperate_by_host.png"),
         p_eff, width = 10, height = 6)
  print(p_eff)
}

# --- B) Fallback: if NO significant hosts, plot ALL eligible ---
if (length(sig_hosts) == 0 && length(eligible_hosts) > 0) {
  message("No significant hosts; plotting all eligible hosts (diagnostic).")
  
  gg_all_box <- ggplot(
    df_host_complete, aes(x = Group, y = Temperate_percent, fill = Group)
  ) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.2, shape = 21, size = 1.3, alpha = 0.6, color = "black") +
    facet_wrap(~ Lowest_host_taxon, scales = "free_y") +
    labs(
      title = "Temperate vOTUs by host — ALL ELIGIBLE (zeros included)",
      x = "Group", y = "Percent temperate vOTUs (%)"
    ) +
    scale_fill_manual(values = c("Mar_Menor" = "lightgreen", "Medit_Sea" = "skyblue")) +
    theme_minimal(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold", size = 9),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
      legend.position = "top"
    )
  ggsave(file.path(output_dir, "faceted_boxplot_temperate_by_host_ELIGIBLE_all.png"),
         gg_all_box, width = 12, height = 8)
  print(gg_all_box)
  
  medians_elig <- medians_diff_all %>%
    left_join(wilcox_all %>% select(Lowest_host_taxon, p_adj), by = "Lowest_host_taxon") %>%
    mutate(
      significance = case_when(
        p_adj < 0.001 ~ "***",
        p_adj < 0.01  ~ "**",
        p_adj < 0.05  ~ "*",
        TRUE ~ ""
      ),
      Label = fct_reorder(Lowest_host_taxon, effect_pp),
      hjust_pos = if_else(effect_pp >= 0, -0.2, 1.2)
    )
  
  p_eff_all <- ggplot(medians_elig, aes(x = effect_pp, y = Label, fill = effect_pp)) +
    geom_col(width = 1, color = "black", linewidth = 0.3) +
    geom_text(aes(label = significance, hjust = hjust_pos), size = 4) +
    scale_fill_gradient2(low = "#9ecae1", mid = "white", high = "#a1d99b",
                         midpoint = 0, limits = range(medians_elig$effect_pp),
                         name = "Δ median (pp)") +
    labs(x = "Δ median (pp) · Mar_Menor − Medit_Sea",
         y = "Host taxon",
         title = "Median difference by host — ALL ELIGIBLE (zeros included)") +
    theme_minimal(base_size = 12) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3))
  ggsave(file.path(output_dir, "effect_median_diff_barplot_temperate_by_host_ELIGIBLE_all.png"),
         p_eff_all, width = 10, height = 6)
  print(p_eff_all)
}

# =========================
# 7) Save all results in Excel
# =========================
wb <- openxlsx::createWorkbook()

openxlsx::addWorksheet(wb, "Global_EXCLUSIVE")
openxlsx::writeData(wb, "Global_EXCLUSIVE", global_wilcox_df)

openxlsx::addWorksheet(wb, "By_host_ALL")
openxlsx::writeData(wb, "By_host_ALL", wilcox_all)

openxlsx::addWorksheet(wb, "By_host_FILTERED")
openxlsx::writeData(wb, "By_host_FILTERED", wilcox_filtered)

openxlsx::addWorksheet(wb, "Per_sample_GLOBAL_input")
openxlsx::writeData(
  wb, "Per_sample_GLOBAL_input",
  df_summary %>% dplyr::select(Sample, Group, Total_vOTUs, Marker_vOTUs, Marker_percent)
)

openxlsx::addWorksheet(wb, "Per_sample_HOST_input")
openxlsx::writeData(
  wb, "Per_sample_HOST_input",
  df_host_complete %>% select(Sample, Group, Lowest_host_taxon,
                              Total_vOTUs, Temperate_vOTUs, Temperate_percent)
)

openxlsx::saveWorkbook(wb, file.path(output_dir, "Lysogeny_tests_and_inputs_EXCLUSIVE.xlsx"), overwrite = TRUE)

cat("\nScript completed. Results saved in:", output_dir, "\n")
