# =========================
# 0) Load libraries
# =========================
library(tidyverse)
library(readr)
library(openxlsx)
library(Hmisc)

# =========================
# 1) Input file path
# =========================
input_file <- "MM_vOTU_abund_table_lysogeny_markers.csv"

# =========================
# 2) Read and clean data
# =========================
# - Expand Lowest_host_taxon when multiple hosts are separated by " OR "
# - Trim spaces and remove rows with "No_host"
df <- read_delim(input_file, delim = ";", locale = locale(decimal_mark = ",")) %>%
  separate_rows(Lowest_host_taxon, sep = " OR ") %>%
  mutate(Lowest_host_taxon = str_trim(Lowest_host_taxon)) %>%
  filter(Lowest_host_taxon != "No_host")

# =========================
# 3) Define temporal order of samples
# =========================
ordered_dates <- c("Oct_2019", "Feb_2020", "Jul_2020", "Oct_2020", "Dec_2020",
                   "Feb_2021", "Apr_2021", "Jul_2021", "Aug_2021", "Nov_2021",
                   "Feb_2022", "May_2022", "Jul_2022", "Oct_2022")

present_date_cols <- intersect(ordered_dates, colnames(df))
if (length(present_date_cols) == 0) {
  stop("No expected date columns found in the input file.")
}

# =========================
# 4) Function to summarize data by host taxon
# =========================
summarize_taxon <- function(taxon_name) {
  df %>%
    pivot_longer(cols = all_of(present_date_cols), names_to = "Sampling_date", values_to = "Abundance") %>%
    filter(Lowest_host_taxon == taxon_name, Abundance > 0) %>%
    group_by(Sampling_date) %>%
    summarise(
      total_vOTUs = n_distinct(vOTU),
      lysogenic_vOTUs = sum(Lysogeny_markers == "Yes"),
      percentage = 100 * lysogenic_vOTUs / total_vOTUs,
      label = paste0(lysogenic_vOTUs, " / ", total_vOTUs),
      Host = gsub("^g__", "", taxon_name),
      .groups = "drop"
    ) %>%
    mutate(Sampling_date = factor(Sampling_date, levels = present_date_cols))
}

# =========================
# 5) Hosts of interest
# =========================
taxa <- c("g__Pelagibacter", "g__Vibrio")

# Generate summaries per host if data are available
summaries <- lapply(taxa, summarize_taxon)
names(summaries) <- taxa
valid_summaries <- summaries[sapply(summaries, nrow) > 0]

# Combine results
combined_summary <- bind_rows(valid_summaries)
if (nrow(combined_summary) == 0) stop("No data available for the specified taxa.")

# Order host factor
present_hosts <- unique(combined_summary$Host)
combined_summary$Host <- factor(
  combined_summary$Host,
  levels = intersect(c("Pelagibacter", "Vibrio"), present_hosts)
)

# =========================
# 6) Plot lysogenic vOTUs (%)
# =========================
plot <- ggplot(combined_summary, aes(x = Sampling_date, y = percentage)) +
  geom_bar(stat = "identity", fill = "lightgreen", color = "black") +
  geom_text(aes(label = label), vjust = -0.3, size = 3) +
  facet_wrap(~ Host, ncol = 1, scales = "fixed") +
  coord_cartesian(ylim = c(0, 20)) +
  labs(title = NULL, x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_blank(),
    axis.text.x = element_text(size = 10, angle = 45, hjust = 1, color = "black"),
    axis.text.y = element_text(size = 10, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank()
  )

print(plot)

# Save figure
ggsave(
  filename = "SAR11_Vibrio_lysogeny_markers_over_time.svg",
  plot = plot, width = 6, height = 4, dpi = 300
)

# =========================
# 7) Save long summary
# =========================
write_delim(
  combined_summary,
  "combined_summary.csv",
  delim = ";"
)

# Print lysogeny percentages by host and date
combined_summary %>%
  select(Host, Sampling_date, percentage) %>%
  arrange(Host, Sampling_date) %>%
  pivot_wider(names_from = Sampling_date, values_from = percentage) %>%
  print(n = Inf, width = Inf)

# =========================
# 8) Correlation with environmental variables
# =========================
env_file <- "environmental_variables.csv"
lysogeny_file <- "combined_summary.csv"

env_df <- read_delim(env_file, delim = ";", locale = locale(decimal_mark = ","))
lysogeny_df <- read_delim(lysogeny_file, delim = ";") %>%
  filter(Host != "ARS21")  # Explicitly exclude ARS21

lysogeny_df <- lysogeny_df %>%
  mutate(Sampling_date = factor(Sampling_date, levels = unique(env_df$Sample)))

# Spearman correlation per host and environmental variable
taxon_cols <- unique(lysogeny_df$Host)

results <- lysogeny_df %>%
  select(Sampling_date, Host, percentage) %>%
  pivot_wider(names_from = Host, values_from = percentage) %>%
  inner_join(env_df, by = c("Sampling_date" = "Sample")) %>%
  pivot_longer(cols = all_of(taxon_cols), names_to = "Host", values_to = "Lysogeny") %>%
  pivot_longer(cols = -c(Sampling_date, Host, Lysogeny), names_to = "Variable", values_to = "Value") %>%
  group_by(Host, Variable) %>%
  summarise(
    rho = suppressWarnings(cor(Lysogeny, Value, method = "spearman", use = "complete.obs")),
    p_value = suppressWarnings(cor.test(Lysogeny, Value, method = "spearman")$p.value),
    .groups = "drop"
  ) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr")) %>%
  arrange(p_adj)

print(results)
cat("\nSignificant correlations (FDR < 0.05):\n")
print(results %>% filter(p_adj < 0.05))

write_csv(
  results,
  "lysogeny_environmental_spearman.csv"
)

# =========================
# 9) Export counts and percentages to Excel
# =========================
output_dir <- "LYSOGENY"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

long_by_host_date <- combined_summary %>%
  select(Host, Sampling_date, lysogenic_vOTUs, total_vOTUs, percentage)

perc_wide <- combined_summary %>%
  select(Host, Sampling_date, percentage) %>%
  pivot_wider(names_from = Sampling_date, values_from = percentage)

lysogenic_counts_wide <- combined_summary %>%
  select(Host, Sampling_date, lysogenic_vOTUs) %>%
  pivot_wider(names_from = Sampling_date, values_from = lysogenic_vOTUs)

total_counts_wide <- combined_summary %>%
  select(Host, Sampling_date, total_vOTUs) %>%
  pivot_wider(names_from = Sampling_date, values_from = total_vOTUs)

wb <- createWorkbook()
addWorksheet(wb, "Long_by_host_date")
writeData(wb, "Long_by_host_date", long_by_host_date)
addWorksheet(wb, "Perc_wide")
writeData(wb, "Perc_wide", perc_wide)
addWorksheet(wb, "Lysogenic_counts_wide")
writeData(wb, "Lysogenic_counts_wide", lysogenic_counts_wide)
addWorksheet(wb, "Total_counts_wide")
writeData(wb, "Total_counts_wide", total_counts_wide)

saveWorkbook(
  wb,
  file.path(output_dir, "combined_summary_counts_perc.xlsx"),
  overwrite = TRUE
)
cat("\nExcel file saved in:", file.path(output_dir, "combined_summary_counts_perc.xlsx"), "\n")

# =========================
# 10) Correlation among host lysogeny percentages
# =========================
host_matrix <- combined_summary %>%
  filter(Host != "ARS21") %>%
  select(Sampling_date, Host, percentage) %>%
  pivot_wider(names_from = Host, values_from = percentage)

cor_matrix <- cor(
  host_matrix %>% select(-Sampling_date),
  use = "pairwise.complete.obs",
  method = "spearman"
)

cor_test <- Hmisc::rcorr(
  as.matrix(host_matrix %>% select(-Sampling_date)),
  type = "spearman"
)

print(cor_matrix)
print(cor_test$P)

rho_df <- as.data.frame(cor_matrix)
rho_df <- tibble::rownames_to_column(rho_df, var = "Host")

pval_df <- as.data.frame(cor_test$P)
pval_df <- tibble::rownames_to_column(pval_df, var = "Host")

wb <- createWorkbook()
addWorksheet(wb, "Spearman_rho")
writeData(wb, "Spearman_rho", rho_df)
addWorksheet(wb, "Spearman_pvals")
writeData(wb, "Spearman_pvals", pval_df)
saveWorkbook(
  wb,
  file = "host_temperate_correlations.xlsx",
  overwrite = TRUE
)
cat("\nCorrelation matrices saved in Excel.\n")
