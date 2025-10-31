## Load required libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(dunn.test)  # For performing Dunn test with Bonferroni correction

# Define file paths
file_path <- "MM_vOTU_abund_table.csv"
output_dir <- "DYNAMIC_CLASSES"

# Output file names
output_pie_chart <- file.path(output_dir, "vOTU_classification_pie_chart.png")
output_boxplot <- file.path(output_dir, "vOTU_abundance_boxplot.png")
output_fraction_boxplot <- file.path(output_dir, "vOTU_fraction_boxplot.png")
output_fraction_barplot <- file.path(output_dir, "vOTU_fraction_barplot.png")

# Create output directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

###############################################################################
## 1. Read data and classify vOTUs
###############################################################################
data <- read.csv(file_path, sep = ";", dec = ",", header = TRUE)
num_samples <- ncol(data) - 1  # The first column is assumed to be vOTU ID

vOTU_classification <- data %>%
  mutate(detections = rowSums(.[, -1] > 0)) %>%
  mutate(classification = case_when(
    detections >= 1 & detections <= 3   ~ "Sporadic",
    detections >= 4 & detections <= 10  ~ "Intermittent",
    detections >= 11 & detections <= 14 ~ "Persistent"
  )) %>%
  mutate(classification = factor(classification,
                                 levels = c("Sporadic", "Intermittent", "Persistent"))) %>%
  select(1, detections, classification)

# Print classification table
cat("Classification table:\n")
print(table(vOTU_classification$classification))
cat("\n('vOTU_classification.csv' will not be saved.)\n")

###############################################################################
## 2. Define colors and sample order
###############################################################################
category_colors <- c(
  "Sporadic"     = "#D9F580",
  "Intermittent" = "lightgreen",
  "Persistent"   = "darkgreen"
)

sample_order <- colnames(data)[2:ncol(data)]

###############################################################################
## 3. Pie chart (global vOTU distribution by type)
###############################################################################
classification_counts <- vOTU_classification %>%
  group_by(classification) %>%
  summarise(count = n())

pie_chart <- ggplot(classification_counts, aes(x = "", y = count, fill = classification)) +
  geom_bar(width = 1, stat = "identity") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = category_colors) +
  theme_void() +
  ggtitle("Distribution of vOTUs by number of detections") +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 16)
  )

print(pie_chart)
ggsave(output_pie_chart, plot = pie_chart, width = 6, height = 6)
cat("\nPie chart saved as 'vOTU_classification_pie_chart.png'.\n")

###############################################################################
## 4. Prepare long-format data and compute log(abundance+1)
###############################################################################
data_long <- data %>%
  pivot_longer(cols = -1, names_to = "sample", values_to = "abundance") %>%
  left_join(vOTU_classification, by = colnames(data)[1]) %>%
  mutate(sample = factor(sample, levels = sample_order)) %>%
  filter(abundance > 0) %>%
  mutate(log_abundance = log(abundance + 1))

###############################################################################
## 5. Boxplot of log(abundance+1) by vOTU type
###############################################################################
boxplot_chart <- ggplot(data_long, aes(x = classification, y = log_abundance, fill = classification)) +
  geom_boxplot(
    outlier.shape = 21,
    outlier.fill = "white",
    outlier.colour = "black",
    outlier.size = 1,
    outlier.stroke = 0.5
  ) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") +
  scale_fill_manual(values = category_colors) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18),
    legend.position = "none"
  ) +
  labs(title = "Distribution of log(abundance+1) by vOTU type",
       x = "vOTU Type",
       y = "Log(Abundance + 1)")

print(boxplot_chart)
ggsave(output_boxplot, plot = boxplot_chart, width = 8, height = 6)
cat("\nAbundance boxplot saved as 'vOTU_abundance_boxplot.png'.\n")

###############################################################################
## 6. Boxplot of vOTU percentage per type in each sample
###############################################################################
df_total <- data_long %>%
  group_by(sample) %>%
  summarise(total_vOTUs = n(), .groups = "drop")

df_types <- data_long %>%
  group_by(sample, classification) %>%
  summarise(n_vOTUs_type = n(), .groups = "drop")

df_fraction <- df_types %>%
  left_join(df_total, by = "sample") %>%
  mutate(percentage = (n_vOTUs_type / total_vOTUs) * 100)

boxplot_fraction <- ggplot(df_fraction, aes(x = classification, y = percentage, fill = classification)) +
  geom_boxplot(
    outlier.shape = 21,
    outlier.fill = "white",
    outlier.colour = "black",
    outlier.size = 1,
    outlier.stroke = 0.5
  ) +
  stat_summary(fun = mean, geom = "point", shape = 20, size = 3, color = "red") +
  scale_fill_manual(values = category_colors) +
  theme_minimal() +
  theme(
    text = element_text(size = 18),
    plot.title = element_text(size = 20, face = "bold"),
    axis.title = element_text(size = 18),
    legend.position = "none"
  ) +
  labs(title = "Percentage of vOTUs per type in each sample",
       x = "vOTU Type",
       y = "Percent of vOTUs in sample")

print(boxplot_fraction)
ggsave(output_fraction_boxplot, plot = boxplot_fraction, width = 8, height = 6)
cat("\nvOTU percentage boxplot saved as 'vOTU_fraction_boxplot.png'.\n")

###############################################################################
## 7. Stacked barplot of vOTU percentages by sample
###############################################################################
df_fraction <- df_fraction %>%
  mutate(sample = factor(sample, levels = sample_order))

barplot_fraction <- ggplot(df_fraction, aes(x = sample, y = percentage, fill = classification)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = category_colors) +
  theme_minimal() +
  labs(title = "Percent of detected vOTUs per type in each sample",
       x = "Time point",
       y = "Percent of vOTUs (%)") +
  theme(
    text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  )

print(barplot_fraction)
ggsave(file.path(output_dir, "dynamic_classes_stackedbarplot.svg"),
       plot = barplot_fraction, width = 9, height = 6, units = "cm")

cat("\nStacked barplot saved as 'vOTU_fraction_barplot.png'.\n")

###############################################################################
## 8. Min/max percentages and top 5 samples per class
###############################################################################
range_percentages <- df_fraction %>%
  group_by(classification) %>%
  summarise(
    min_percentage = round(min(percentage, na.rm = TRUE), 1),
    max_percentage = round(max(percentage, na.rm = TRUE), 1),
    .groups = "drop"
  )
print(range_percentages)

top_bottom_5 <- df_fraction %>%
  group_by(classification) %>%
  arrange(desc(percentage)) %>%
  mutate(rank_desc = row_number()) %>%
  arrange(percentage) %>%
  mutate(rank_asc = row_number()) %>%
  ungroup() %>%
  filter(rank_desc <= 5 | rank_asc <= 5) %>%
  arrange(classification, rank_asc)

cat("\nTop 5 samples with HIGHEST and LOWEST vOTU percentages per type:\n")
print(top_bottom_5 %>% select(classification, sample, percentage, rank_asc, rank_desc))

###############################################################################
## 9. Barplot of vOTU detection frequency
###############################################################################
detection_frequency <- vOTU_classification %>%
  count(detections)

barplot_frequency <- ggplot(detection_frequency, aes(x = factor(detections), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue2", width = 0.9) +
  labs(
    title = "vOTU detection frequency",
    x = "Detections",
    y = "Number of vOTUs"
  ) +
  theme_minimal(base_size = 10) +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3))

print(barplot_frequency)
ggsave(file.path(output_dir, "vOTU_detection_frequency_barplot.svg"),
       plot = barplot_frequency, width = 7.5, height = 4)

cat("\nDetection frequency barplot saved as 'vOTU_detection_frequency_barplot.png'.\n")

###############################################################################
## 10. Summary statistics and Excel export
###############################################################################
total_vOTUs <- nrow(vOTU_classification)
type_counts <- table(vOTU_classification$classification)
type_percent <- round(100 * type_counts / total_vOTUs, 2)

cat("\nSummary of vOTUs by type:\n")
for (t in names(type_counts)) {
  cat(sprintf("- %s: %d vOTUs (%.2f%%)\n", t, type_counts[t], type_percent[t]))
}

# Save tables to Excel
id_col <- names(data)[1]
dynamic_table <- vOTU_classification %>%
  dplyr::rename(vOTU_ID = !!id_col)

stacked_table <- df_fraction %>%
  dplyr::select(sample, classification, percentage) %>%
  dplyr::mutate(sample = factor(sample, levels = sample_order)) %>%
  dplyr::arrange(sample, classification)

xlsx_path <- file.path(output_dir, "vOTU_dynamic_type_and_stackedbar.xlsx")

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "vOTU_dynamic_type")
openxlsx::writeData(wb, "vOTU_dynamic_type", dynamic_table)
openxlsx::addWorksheet(wb, "stacked_barplot_data")
openxlsx::writeData(wb, "stacked_barplot_data", stacked_table)
openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)
message("Tables saved in: ", xlsx_path)
