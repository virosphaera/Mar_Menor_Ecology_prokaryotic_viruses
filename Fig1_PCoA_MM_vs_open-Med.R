# ============================================================================================================
# Principal Coordinate and Silhouette Analyses for k-means Clustering: Mar Menor vs. Open Mediterranean Sea
# ============================================================================================================

# --- Install (if needed) ---
# install.packages("remotes")
# remotes::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")


# --- Packages ---
library(compositions)
library(vegan)
library(ggplot2)
library(stringr)
library(dplyr)
library(ggrepel)
library(cluster)
library(tidyr)
library(openxlsx)
library(pairwiseAdonis)


# --- Parameters ---
point_size <- 5
label_text_size <- 4
axis_text_size <- 11
axis_title_size <- 13
title_size <- 13

# --- Load data ---
file_path <- "Med-MM_vOTU_abund_table.csv"
vOTU <- read.table(file_path, sep = ";", dec = ",", header = TRUE, row.names = 1)
trans_df <- t(vOTU)

# --- CLR + Aitchison distance ---
df_clr <- clr(trans_df + 0.001)
dist_aitchison <- dist(df_clr, method = "euclidean")

# --- 1) Silhouette vs k ---
ks <- 2:10
set.seed(123)
sil_means <- sapply(ks, function(k) {
  km <- kmeans(df_clr, centers = k, nstart = 50)
  mean(silhouette(km$cluster, dist_aitchison)[, 3])
})
metrics_df <- data.frame(k = ks, silhouette_mean = sil_means)

p_sil_Med_MM <- ggplot(metrics_df, aes(k, silhouette_mean)) +
  # Line connecting the points
  geom_line(color = "black", linewidth = 0.9) +
  # Points
  geom_point(color = "darkgrey", size = 4) +
  theme_minimal(base_size = 14) +
  theme(
    # Plot border
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    panel.background = element_blank(),
    # Axis text
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    # Axis titles
    axis.title.x = element_text(size = 14, color = "black", face = "plain"),
    axis.title.y = element_text(size = 14, color = "black", face = "plain")
  ) +
  labs(
    x = "Number of clusters (k)",
    y = "Average silhouette width"
  )

print(p_sil_Med_MM)

ggplot2::ggsave(
  "Silhouette_vs_k_MM-Med.svg",
  p_sil_Med_MM, width = 8, height = 8, units = "cm", dpi = 300
)

# --- 2) PCoA for k = 4 ---
set.seed(123)
k <- 4
km <- kmeans(df_clr, centers = k, nstart = 50)
cluster_labels <- factor(km$cluster)
perc_SS <- round(100 * km$betweenss / km$totss, 2)

sample_names <- rownames(trans_df)
ecosystem <- ifelse(str_detect(sample_names, "^[A-Z][a-z]{2}_[0-9]{4}"),
                    "Mar_Menor", "Medit_Sea")

custom_labels <- case_when(
  str_detect(sample_names, "^TARA_007") ~ "Algerian",
  str_detect(sample_names, "^TARA_009") ~ "Men-Sar",
  str_detect(sample_names, "^TARA_018") ~ "Sicily",
  str_detect(sample_names, "^TARA_023") ~ "Adriatic",
  str_detect(sample_names, "^TARA_025") ~ "Ionian",
  str_detect(sample_names, "^TARA_030") ~ "Levantine",
  str_detect(sample_names, "^Med") ~ "near-MM",
  TRUE ~ sample_names
)

pcoa <- cmdscale(dist_aitchison, k = 2, eig = TRUE)
pcoa_df <- as.data.frame(pcoa$points)
colnames(pcoa_df) <- c("PCoA1", "PCoA2")
pcoa_df$Sample <- sample_names
pcoa_df$Ecosystem <- factor(ecosystem, levels = c("Mar_Menor", "Medit_Sea"))
pcoa_df$Label <- custom_labels
pcoa_df$Cluster <- cluster_labels

eig_values <- pcoa$eig
var_explained <- round(100 * eig_values[1:2] / sum(eig_values[eig_values > 0]), 2)

shape_vals <- c(21, 22, 24, 23)  # circle, square, triangle, diamond
fill_vals  <- c("Mar_Menor" = "lightgreen", "Medit_Sea" = "skyblue")

p_pcoa <- ggplot(pcoa_df, aes(PCoA1, PCoA2, fill = Ecosystem, shape = Cluster)) +
  geom_point(size = point_size, stroke = 0.6, color = "black") +
  geom_text_repel(aes(label = Label), size = label_text_size, max.overlaps = Inf) +
  scale_fill_manual(values = fill_vals) +
  scale_shape_manual(values = shape_vals) +
  theme_minimal(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    axis.text = element_text(size = axis_text_size),
    axis.title = element_text(size = axis_title_size),
    plot.title = element_text(size = title_size, face = "bold")
  ) +
  labs(
    title = paste0("PCoA (Aitchison) | Ecosystem (fill) & k-means (shape, k = ", k, ") | %SS = ", perc_SS),
    x = paste0("PCoA1 (", var_explained[1], "%)"),
    y = paste0("PCoA2 (", var_explained[2], "%)")
  )
print(p_pcoa)

# Perform PERMANOVA for the clusters
adonis2(dist_aitchison ~ Cluster, data = pcoa_df)

# Post-hoc PERMANOVA (pairwise)
pairwise.adonis(dist_aitchison, pcoa_df$Cluster, p.adjust.m = "fdr")

# --- 3) Export data from both plots ---
silhouette_table <- metrics_df
pcoa_table <- pcoa_df %>% select(Sample, Ecosystem, Cluster, PCoA1, PCoA2)

wb <- createWorkbook()
addWorksheet(wb, "Silhouette_vs_k")
writeData(wb, "Silhouette_vs_k", silhouette_table)
addWorksheet(wb, "PCoA_k4")
writeData(wb, "PCoA_k4", pcoa_table)
saveWorkbook(wb, "silhouette_and_pcoa_data.xlsx", overwrite = TRUE)


# --- Permutation test for several k values ---
perm_test_kmeans <- function(df_clr, dist_mat, k, n_perm = 199, nstart = 10) {
  # Observed silhouette
  km_obs <- kmeans(df_clr, centers = k, nstart = nstart)
  sil_obs <- mean(silhouette(km_obs$cluster, dist_mat)[, 3])
  
  # Permutations (only shuffling sample order)
  sil_perm <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    df_perm <- df_clr[sample(nrow(df_clr)), , drop = FALSE]
    km_perm <- kmeans(df_perm, centers = k, nstart = nstart)
    sil_perm[i] <- mean(silhouette(km_perm$cluster, dist_mat)[, 3])
  }
  
  p_val <- (sum(sil_perm >= sil_obs) + 1) / (n_perm + 1)
  return(list(silhouette_obs = sil_obs, p_value = p_val))
}

# Example: test for several k values and save results
ks <- 2:6
n_perm <- 199  # increase to 999 for final version if runtime allows

perm_results <- lapply(ks, function(k) {
  res <- perm_test_kmeans(df_clr, dist_aitchison, k, n_perm = n_perm)
  data.frame(k = k, silhouette_mean = res$silhouette_obs, p_value_perm = res$p_value)
})
perm_results_df <- do.call(rbind, perm_results)

print(perm_results_df)


# --- Tables to export ---
silhouette_table <- metrics_df
pcoa_table <- pcoa_df %>% dplyr::select(Sample, Ecosystem, Cluster, PCoA1, PCoA2)

# Silhouette per sample for k = 4
sil_k4 <- cluster::silhouette(km$cluster, dist_aitchison)
sil_k4_df <- data.frame(
  Sample    = rownames(df_clr),
  Cluster   = factor(km$cluster),
  sil_width = sil_k4[, 3],
  row.names = NULL
)

# K-means summary for k = 4
kmeans_summary <- data.frame(
  k            = 4,
  totss        = km$totss,
  betweenss    = km$betweenss,
  tot_withinss = km$tot.withinss,
  perc_SS      = round(100 * km$betweenss / km$totss, 2),
  iter         = km$iter
)
kmeans_sizes <- data.frame(Cluster = seq_along(km$size), Size = km$size)

# PCoA eigenvalues
eig_df <- data.frame(
  PC            = seq_along(pcoa$eig),
  eigenvalue    = pcoa$eig,
  var_explained = 100 * pcoa$eig / sum(pcoa$eig[pcoa$eig > 0])
)

# --- Export to Excel (one sheet per table + extras) ---
output_dir <- "PARTITION_PLOTS"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
xlsx_path <- file.path(output_dir, "silhouette_and_pcoa_data_ecosystems.xlsx")

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, "Silhouette_vs_k")
openxlsx::writeData(wb, "Silhouette_vs_k", silhouette_table)

openxlsx::addWorksheet(wb, "PCoA_k4")
openxlsx::writeData(wb, "PCoA_k4", pcoa_table)

openxlsx::addWorksheet(wb, "Silhouette_per_sample_k4")
openxlsx::writeData(wb, "Silhouette_per_sample_k4", sil_k4_df)

openxlsx::addWorksheet(wb, "kmeans_k4_summary")
openxlsx::writeData(wb, "kmeans_k4_summary", kmeans_summary)
openxlsx::writeData(wb, "kmeans_k4_summary", kmeans_sizes, startRow = nrow(kmeans_summary) + 3)

openxlsx::addWorksheet(wb, "PCoA_eigenvalues")
openxlsx::writeData(wb, "PCoA_eigenvalues", eig_df)

openxlsx::saveWorkbook(wb, xlsx_path, overwrite = TRUE)

# (optional) save figures in PNG format as well
ggplot2::ggsave(file.path(output_dir, "PCoA_k4.png"), plot = p_pcoa, width = 8, height = 6, dpi = 300)

message("Tables saved to: ", xlsx_path)
