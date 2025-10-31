# ===========================================================================
# Load libraries
# ===========================================================================
library(dplyr)
library(tidyr)
library(stringr)
library(readr)
library(pheatmap)
library(viridis)
library(tibble)
library(grid)
library(purrr)

# ===========================================================================
# Define paths
# ===========================================================================
input_file <- "MM_vOTU_abund_table.csv"
env_file <- "environmental_variables.csv"
out_dir <- "HOSTS"

# ===========================================================================
# Load and process abundance data
# ===========================================================================
df <- read_csv2(input_file)

# Convert comma decimals to numeric
abund_cols <- colnames(df)[!(colnames(df) %in% c("vOTU", "Host"))]
df[abund_cols] <- lapply(df[abund_cols], function(x) as.numeric(str_replace(x, ",", ".")))

sample_order <- abund_cols

df_long <- df %>%
  separate_rows(Host, sep = " OR ") %>%
  mutate(Lowest_taxon = ifelse(is.na(Host) | Host == "", "No_host", Host)) %>%
  select(-Host)

abund_long <- df_long %>%
  pivot_longer(cols = all_of(sample_order), names_to = "Sample", values_to = "Abundance") %>%
  filter(Abundance > 0) %>%
  group_by(Sample, vOTU) %>%
  mutate(nTaxa = n()) %>%
  ungroup() %>%
  mutate(Fraction_abundance = Abundance / nTaxa) %>%
  distinct(Sample, vOTU, Lowest_taxon, .keep_all = TRUE) %>%
  mutate(Sample = factor(Sample, levels = sample_order))

# ===========================================================================
# HEATMAP: Top 20 host taxa by number of vOTUs (excluding No_host)
# ===========================================================================
abund_long_hm <- abund_long %>% filter(Lowest_taxon != "No_host")

top_taxa <- abund_long_hm %>%
  group_by(Lowest_taxon) %>%
  summarise(Num_vOTUs = n_distinct(vOTU), .groups = "drop") %>%
  arrange(desc(Num_vOTUs)) %>%
  slice_head(n = 20) %>%
  pull(Lowest_taxon)

heat_data <- abund_long_hm %>%
  filter(Lowest_taxon %in% top_taxa) %>%
  group_by(Lowest_taxon, Sample) %>%
  summarise(Total_abund = sum(Fraction_abundance), .groups = "drop") %>%
  pivot_wider(names_from = Sample, values_from = Total_abund, values_fill = 0)

heat_matrix <- heat_data %>%
  column_to_rownames("Lowest_taxon") %>%
  as.matrix()

rownames(heat_matrix) <- gsub("^[a-z]__+", "", rownames(heat_matrix))
heat_matrix_log <- log10(heat_matrix + 1)

# ===========================================================================
# Plot heatmap
# ===========================================================================
plot_heatmap <- pheatmap(
  heat_matrix_log,
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  color         = viridis(100),
  fontsize      = 10,
  angle_col     = 45,
  main          = "Log(x+1) relative abundance (top 20 host taxa)",
  filename      = NA
)

print(plot_heatmap)

# Customize axis text formatting
library(pheatmap)
library(grid)
library(viridis)

plot_heatmap <- pheatmap(
  heat_matrix_log,
  cluster_rows  = TRUE,
  cluster_cols  = FALSE,
  color         = viridis(100),
  fontsize      = 10,
  fontsize_row  = 10,
  fontsize_col  = 10,
  angle_col     = 45,
  main          = "Log(x+1) relative abundance (top 20 host taxa)",
  filename      = NA,
  plot          = FALSE
)

# Identify available grobs
g_names <- vapply(plot_heatmap$gtable$grobs, function(g) g$name, character(1))

row_idx <- grep("row_names", g_names)
col_idx <- grep("col_names", g_names)

# Apply color (and size if desired) to all matches
if (length(row_idx)) {
  for (i in row_idx) plot_heatmap$gtable$grobs[[i]]$gp <- gpar(col = "black", fontsize = 9)
}
if (length(col_idx)) {
  for (i in col_idx) plot_heatmap$gtable$grobs[[i]]$gp <- gpar(col = "black", fontsize = 8)
}

grid.newpage()
grid.draw(plot_heatmap$gtable)

ggsave(file.path(out_dir, "Heatmap_top20_host_taxa_over_time.svg"), plot_heatmap, width = 17, height = 15, units = "cm", dpi = 300)

# ===========================================================================
# Spearman correlation with environmental variables (FDR corrected)
# ===========================================================================
# Load environmental data
env_df <- read_csv2(env_file)

# Convert to numeric (replace comma with dot)
env_df <- env_df %>%
  mutate(across(-Sample, ~ as.numeric(str_replace(., ",", "."))))

# Prepare abundance matrix: taxa x Sample
abundance_matrix <- heat_data %>%
  column_to_rownames("Lowest_taxon")

# Transpose to Sample x Taxon and align with env_df
abundance_for_corr <- as.data.frame(t(abundance_matrix)) %>%
  rownames_to_column("Sample")

# Join environmental and abundance data
corr_df <- inner_join(env_df, abundance_for_corr, by = "Sample")

# Get variable names
env_vars <- colnames(env_df)[-1]
taxon_vars <- colnames(abundance_for_corr)[-1]

# Run correlations
corr_results <- map_dfr(env_vars, function(env_var) {
  map_dfr(taxon_vars, function(taxon) {
    test <- cor.test(corr_df[[env_var]], corr_df[[taxon]], method = "spearman")
    tibble(
      Taxon = taxon,
      Env_var = env_var,
      Rho = test$estimate,
      P = test$p.value
    )
  })
}) %>%
  mutate(P_adjusted = p.adjust(P, method = "fdr")) %>%
  arrange(P_adjusted)

# Show top significant results
print(corr_results %>% filter(P_adjusted <= 0.05))

warnings()

# Save all results
write_csv2(corr_results,
           "spearman_all_results.csv")

# Save only significant results (FDR ≤ 0.05)
sig_corr <- corr_results %>% filter(P_adjusted <= 0.05)
write_csv2(sig_corr,
           "spearman_significant_results.csv")
