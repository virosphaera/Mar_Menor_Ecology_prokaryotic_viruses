# ==== LOAD LIBRARIES ====
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(ggplot2)
library(grid)
library(viridis)

# ==== PATHS ====
input_abund <- "MM_vOTU_abund_table_w_functions_wo_underscores.csv"
output_dir <- "OUTPUTS/AMGs/"

# ==== READ DATA ====
abund <- read.csv(input_abund, sep = ";", check.names = FALSE)

# ==== CONVERT COMMAS TO DOTS IN NUMERIC COLUMNS ====
abund[, 4:ncol(abund)] <- lapply(abund[, 4:ncol(abund)], function(x) as.numeric(str_replace(x, ",", ".")))

# ==== USE COLUMN 3 AS AMG FUNCTION ====
colnames(abund)[3] <- "AMG"

# ==== EXPAND MULTIPLE AMGs ("AND" or "OR") ====
if (!"AMG" %in% colnames(abund)) stop("Column 'AMG' not found")

abund_expanded <- abund %>%
  filter(!is.na(AMG) & AMG != "" & AMG != "No_AMG") %>%
  mutate(AMG = str_trim(AMG)) %>%
  separate_rows(AMG, sep = "\\s*(AND|OR)\\s*") %>%
  filter(!is.na(AMG) & AMG != "" & AMG != "No_AMG")

if (nrow(abund_expanded) == 0) stop("No valid AMGs after expansion.")

# ==== CALCULATE NUMBER OF vOTUs PER AMG ====
votu_counts <- abund_expanded %>%
  group_by(AMG) %>%
  summarise(n_vOTUs = n_distinct(vOTU)) %>%
  arrange(desc(n_vOTUs))

# ==== GROUP ABUNDANCES BY AMG ====
# Explicitly exclude the 'AMG_count' column from the summary
abund_grouped <- abund_expanded %>%
  select(-AMG_count) %>%
  group_by(AMG) %>%
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%
  ungroup() %>%
  left_join(votu_counts, by = "AMG") %>%
  arrange(desc(n_vOTUs))

# ==== LOG(x+1) TRANSFORMATION ====
sample_names <- setdiff(colnames(abund_grouped), c("AMG", "n_vOTUs"))
abund_log <- abund_grouped
abund_log[, sample_names] <- log1p(abund_grouped[, sample_names])

# ==== HEATMAP ====
mat <- as.matrix(abund_log[, sample_names])
rownames(mat) <- abund_log$AMG

# Row clustering
clust_rows <- hclust(dist(mat), method = "ward.D2")

# Create heatmap
heatmap_file <- file.path(output_dir, "heatmap_logx1_AMG_abundances_clustered.png")
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

heatmap_plot <- pheatmap(mat,
                         scale = "none",
                         cluster_rows = clust_rows,
                         cluster_cols = FALSE,
                         show_rownames = TRUE,
                         color = viridis(100),
                         angle_col = 45,
                         main = "log(x+1) AMG abundances clustered",
                         silent = TRUE)

# Display and save
grid::grid.newpage()
grid::grid.draw(heatmap_plot$gtable)

png(heatmap_file, width = 1000, height = 1200, res = 150)
grid::grid.newpage()
grid::grid.draw(heatmap_plot$gtable)
dev.off()

message("Heatmap saved to: ", heatmap_file)

ggsave("viral_abundance_per_AMG_functional_category_over_time.svg", 
       heatmap_plot, dpi=300, width=19, height=23, units="cm")


# ==== LOAD ENVIRONMENTAL DATA ====
env_path <- "environmental_variables.csv"

env_data <- read_delim(env_path, delim = ";", locale = locale(decimal_mark = ",")) %>%
  mutate(across(-Sample, ~ as.numeric(gsub(",", ".", as.character(.)))))  # convert comma decimals to numeric

# ==== PREPARE VIRAL FUNCTION ABUNDANCE MATRIX ====

# Filter only samples with significant functions
viral_sig <- vOTUs_percent %>%
  filter(FUNCTION_ID %in% sig_functions) %>%
  select(Sample, FUNCTION_ID, percent_vOTUs) %>%
  pivot_wider(names_from = FUNCTION_ID, values_from = percent_vOTUs)

# ==== JOIN ENVIRONMENTAL DATA WITH VIRAL DATA BY SAMPLE ====
cor_data <- env_data %>%
  inner_join(viral_sig, by = "Sample")

# ==== ENVIRONMENTAL VARIABLES AND FUNCTIONS ====
env_vars <- colnames(env_data)[-1]
func_vars <- colnames(viral_sig)[-1]

# ==== INITIALIZE RESULTS LIST ====
cor_results <- list()

# ==== CALCULATE SPEARMAN CORRELATIONS WITH FDR CORRECTION ====
for (func in func_vars) {
  for (env_var in env_vars) {
    x <- cor_data[[func]]
    y <- cor_data[[env_var]]
    
    # Remove NA pairs
    valid_idx <- !is.na(x) & !is.na(y)
    if (sum(valid_idx) > 3) { # minimum number of points to correlate
      cor_test <- cor.test(x[valid_idx], y[valid_idx], method = "spearman")
      
      cor_results <- append(cor_results, list(
        tibble(
          Function = func,
          Environmental_Variable = env_var,
          Spearman_rho = cor_test$estimate,
          p_value = cor_test$p.value
        )
      ))
    }
  }
}

cor_df <- bind_rows(cor_results)

# ==== FDR CORRECTION PER FUNCTION (BH) ====
cor_df <- cor_df %>%
  group_by(Function) %>%
  mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
  ungroup() %>%
  arrange(p_adj)

# ==== PRINT SIGNIFICANT CORRELATIONS ====
cat("\n=== Significant Spearman correlations (FDR < 0.05) ===\n\n")
sig_correlations <- cor_df %>% filter(p_adj < 0.05)
print(sig_correlations, n = Inf)

# ==== SAVE TO CSV ====
write_csv2(sig_correlations, file.path(output_dir, "spearman_correlation_significant_FDR.csv"))


# ==== GET TOP 1 DATE WITH HIGHEST ABUNDANCE FOR EACH AMG ====

# Raw abundances (before log(x+1)) for each AMG and sample
abund_raw <- abund_grouped %>%
  select(AMG, all_of(sample_names))

# Transform to long format: AMG - Sample - Abundance
abund_long <- abund_raw %>%
  pivot_longer(cols = -AMG, names_to = "Sample", values_to = "Abundance")

# For each AMG, get the top 1 sample with the highest abundance
top3_abund <- abund_long %>%
  group_by(AMG) %>%
  slice_max(order_by = Abundance, n = 1, with_ties = FALSE) %>%
  arrange(AMG, desc(Abundance)) %>%
  ungroup()

# Display result
cat("\nTop sample with the highest abundance per AMG:\n")
print(top3_abund)

# Optional: save to CSV
write_csv2(top3_abund, file.path(output_dir, "top3_abundances_per_AMG.csv"))
