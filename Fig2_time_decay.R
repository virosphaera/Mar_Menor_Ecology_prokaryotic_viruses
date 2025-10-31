library(ggplot2)
library(reshape2)
library(dplyr)
library(vegan)
library(lubridate)
library(betapart)

# 1. Load data
input_dir <- "C:/Users/guillermo.dominguez/Documents/EPRA/VIRUS_MAR_MENOR/MM_prok_viruses_dynamics_and_comparison_w_Medit/DATA/ABUNDANCE/"
file_path <- paste0(input_dir, "MM_vOTU_abund_table.csv")
data <- read.csv(file_path, sep=";", dec=",", row.names=1, check.names=FALSE)
output_dir <- "C:/Users/guillermo.dominguez/Documents/EPRA/VIRUS_MAR_MENOR/MM_prok_viruses_dynamics_and_comparison_w_Medit/OUTPUTS/TIME_DECAY"

# 2. Initial formatting
data$vOTU <- rownames(data)
rownames(data) <- NULL  
data_long <- melt(data, id.vars = "vOTU", variable.name = "Sample", value.name = "Abundance")
data_long$LogAbundance <- log(data_long$Abundance + 1)
data_long <- subset(data_long, LogAbundance > 0)

# 3. Dates
data_long$Month <- as.numeric(match(substr(data_long$Sample, 1, 3), month.abb))
data_long$Year <- as.numeric(substr(data_long$Sample, 5, 8))
data_long$Date <- as.Date(paste0(data_long$Year, "-", data_long$Month, "-01"))
min_date <- min(data_long$Date, na.rm = TRUE)
data_long$MonthsSinceStart <- as.numeric(difftime(data_long$Date, min_date, units = "days")) / 30  

# 4. Abundance matrix
abundance_matrix <- dcast(data_long, Sample ~ vOTU, value.var = "LogAbundance", fun.aggregate = mean)
abundance_matrix[is.na(abundance_matrix)] <- 0  
sample_labels <- abundance_matrix$Sample
abundance_matrix <- abundance_matrix[,-1]  
rownames(abundance_matrix) <- sample_labels

# 5. Bray–Curtis partitioning
bray_parts <- bray.part(abundance_matrix)
bc_total <- as.matrix(bray_parts$bray)
bc_turnover <- as.matrix(bray_parts$bray.bal)
bc_grad <- as.matrix(bray_parts$bray.gra)

# 6. Sample pairs and temporal distances
sample_pairs <- expand.grid(Sample1 = sample_labels, Sample2 = sample_labels)
sample_pairs <- sample_pairs[sample_pairs$Sample1 != sample_pairs$Sample2, ]
sample_pairs$TimeDifference <- abs(data_long$MonthsSinceStart[match(sample_pairs$Sample1, data_long$Sample)] - 
                                     data_long$MonthsSinceStart[match(sample_pairs$Sample2, data_long$Sample)])
sample_pairs$TimeDifference <- round(sample_pairs$TimeDifference)

# 7. Function to extract values from symmetric matrices
get_lower_tri_values <- function(mat, labels) {
  vals <- NULL
  for (i in 1:(nrow(mat)-1)) {
    for (j in (i+1):ncol(mat)) {
      vals <- rbind(vals, data.frame(Sample1 = labels[i], Sample2 = labels[j], Value = mat[i, j]))
    }
  }
  return(vals)
}

# 8. Extract dissimilarities for each component
dist_total <- get_lower_tri_values(bc_total, sample_labels)
dist_turn <- get_lower_tri_values(bc_turnover, sample_labels)
dist_grad <- get_lower_tri_values(bc_grad, sample_labels)

# 9. Combine with time differences
sample_pairs_key <- sample_pairs %>%
  mutate(pair_id = paste(pmin(Sample1, Sample2), pmax(Sample1, Sample2), sep = "_"))

dist_total <- dist_total %>%
  mutate(pair_id = paste(pmin(Sample1, Sample2), pmax(Sample1, Sample2), sep = "_"),
         Type = "Total") %>%
  rename(Distance = Value)

dist_turn <- dist_turn %>%
  mutate(pair_id = paste(pmin(Sample1, Sample2), pmax(Sample1, Sample2), sep = "_"),
         Type = "Turnover") %>%
  rename(Distance = Value)

dist_grad <- dist_grad %>%
  mutate(pair_id = paste(pmin(Sample1, Sample2), pmax(Sample1, Sample2), sep = "_"),
         Type = "Abundance") %>%
  rename(Distance = Value)

all_distances <- bind_rows(dist_total, dist_turn, dist_grad)

# Add TimeDifference
all_distances <- left_join(all_distances, sample_pairs_key[, c("pair_id", "TimeDifference")], by = "pair_id")

# 10. Correlations by component
for (type in unique(all_distances$Type)) {
  cat(paste0("\n--- Spearman correlation for ", type, " ---\n"))
  d <- subset(all_distances, Type == type)
  print(cor.test(d$TimeDifference, d$Distance, method = "spearman"))
}

# 11. Comparative plot
summary_plot_data <- all_distances %>%
  group_by(Type, TimeDifference) %>%
  summarise(MeanDistance = mean(Distance, na.rm = TRUE), .groups = "drop")

ggplot(all_distances, aes(x = TimeDifference, y = Distance)) +
  geom_point(alpha = 0.3, color = "gray") +
  geom_line(data = summary_plot_data, aes(x = TimeDifference, y = MeanDistance, color = Type), linewidth = 1.2) +
  geom_point(data = summary_plot_data, aes(x = TimeDifference, y = MeanDistance, fill = Type),
             size = 3, shape = 21, color = "black") +
  facet_wrap(~Type, scales = "free_y") +
  theme_minimal() +
  labs(title = "Temporal decay of community similarity components",
       x = "Time Difference (months)", y = "Dissimilarity") +
  theme(strip.text = element_text(face = "bold"), axis.text = element_text(size = 10))

# Create a list to store plots by type
plots_by_type <- list()

# Custom colors by type (optional)
type_colors <- c("Total" = "plum", "Turnover" = "lightgreen", "Abundance" = "blue2")

# Generate individual plots for each component
for (t in unique(all_distances$Type)) {
  dist_data <- subset(all_distances, Type == t)
  summary_data <- subset(summary_plot_data, Type == t)
  
  p <- ggplot(dist_data, aes(x = TimeDifference, y = Distance)) +
    geom_point(alpha = 0.3, color = "gray") +
    geom_line(data = summary_data, aes(x = TimeDifference, y = MeanDistance), 
              color = type_colors[t], linewidth = 1.2) +
    geom_point(data = summary_data, aes(x = TimeDifference, y = MeanDistance), 
               size = 1.5, shape = 21, fill = type_colors[t], color = "black") +
    theme_minimal() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      panel.grid = element_blank(),
      axis.text = element_text(size = 8, color = "black"),
      axis.title = element_text(size = 9),
      plot.title = element_text(face = "bold", hjust = 0.5)
    ) +
    labs(title = paste("Temporal decay -", t),
         x = "Time Difference (months)", y = "Dissimilarity")
  
  # Store each plot
  plots_by_type[[t]] <- p
}

# Display each plot separately
print(plots_by_type$"Abundance")
print(plots_by_type$"Total")
print(plots_by_type$"Turnover")

# Save each plot as SVG
ggsave(file.path(output_dir, "decay_Total.svg"),
       plot = plots_by_type[["Total"]], width = 7.5, height = 4, units = "cm", dpi=300)

ggsave(file.path(output_dir, "decay_Turnover.svg"),
       plot = plots_by_type[["Turnover"]], width = 7.5, height = 4, units = "cm", dpi=300)

ggsave(file.path(output_dir, "decay_Abundance.svg"),
       plot = plots_by_type[["Abundance"]], width = 7.67, height = 4, units = "cm", dpi=300)


# ============================
# Mantel tests by component
# ============================

library(vegan)

# 1. Temporal distance matrix
sample_dates <- data_long %>%
  select(Sample, Date) %>%
  distinct() %>%
  arrange(Sample)
sample_dates <- sample_dates[match(rownames(abundance_matrix), sample_dates$Sample), ]
time_vector <- as.numeric(sample_dates$Date)
time_matrix <- dist(time_vector) / 30  # convert days to months

# 2. Function to run Mantel test on a dissimilarity matrix
run_mantel <- function(dist_matrix, time_matrix, method = "spearman", permutations = 9999, label = "") {
  cat(paste0("\n========================\n"))
  cat(paste0("Mantel Test: ", label, "\n"))
  result <- mantel(dist_matrix, time_matrix, method = method, permutations = permutations)
  print(result)
  return(result)
}

# 3. Run the three tests
mantel_total <- run_mantel(bc_total, time_matrix, label = "Bray-Curtis Total")
mantel_turnover <- run_mantel(bc_turnover, time_matrix, label = "Turnover (bray.bal)")
mantel_abundance <- run_mantel(bc_grad, time_matrix, label = "Abundance Gradient (bray.gra)")


# ==============================================
# Save tables with results from decay plots and Mantel tests
# ==============================================

# 1. Save decay table: Mean dissimilarity vs TimeDifference per component
decay_table_path <- file.path(output_dir, "temporal_decay_summary_by_component.csv")
write.csv(summary_plot_data, decay_table_path, row.names = FALSE)
cat(paste0("Saved decay summary table to:\n", decay_table_path, "\n"))

# 2. Create data frame with Mantel test results
mantel_results_df <- data.frame(
  Component = c("Total", "Turnover", "Abundance"),
  Mantel_r = c(
    mantel_total$statistic,
    mantel_turnover$statistic,
    mantel_abundance$statistic
  ),
  P_value = c(
    mantel_total$signif,
    mantel_turnover$signif,
    mantel_abundance$signif
  ),
  Permutations = c(
    mantel_total$permutations,
    mantel_turnover$permutations,
    mantel_abundance$permutations
  )
)

# Save Mantel test results
mantel_table_path <- file.path(output_dir, "mantel_test_results_by_component.csv")
write.csv(mantel_results_df, mantel_table_path, row.names = FALSE)
cat(paste0("Saved Mantel test results to:\n", mantel_table_path, "\n"))
