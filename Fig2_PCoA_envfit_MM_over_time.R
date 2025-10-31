# --- Libraries ---
library(compositions)
library(vegan)
library(ggplot2)
library(ggrepel)
library(stringr)
library(cluster)
library(dplyr)
library(grid)        # for unit()
library(openxlsx)    # for exporting tables
library(pairwiseAdonis)

# --- Parameters ---
point_size <- 5
label_text_size <- 3
axis_text_size <- 10
axis_title_size <- 11
title_size <- 16
k_fix <- 3  # fixed k for PCoA

# --- Paths ---
abundance_file <- "MM_vOTU_abund_TABLE.csv"
env_file       <- "environmental_variables.csv"

# --- Load data ---
vOTU     <- read.table(abundance_file, sep=";", dec=",", header=TRUE, row.names=1)
metadata <- read.table(env_file, sep=";", dec=",", header=TRUE, row.names=1)

# --- CLR + Aitchison ---
df_clr <- clr(t(vOTU) + 0.001)

# Align samples (intersection)
common <- intersect(rownames(df_clr), rownames(metadata))
stopifnot(length(common) >= 4)
df_clr   <- df_clr[common, , drop=FALSE]
metadata <- metadata[common, , drop=FALSE]

# Aitchison distances
dist_aitchison <- dist(df_clr, method="euclidean")

# ========= 1) Silhouette vs k (2–10) =========
ks <- 2:min(10, nrow(df_clr)-1)
set.seed(123)
sil_means <- sapply(ks, function(k){
  km <- kmeans(df_clr, centers=k, nstart=50)
  mean(cluster::silhouette(km$cluster, dist_aitchison)[,3])
})
metrics_df <- data.frame(k=ks, silhouette_mean=sil_means)

p_sil_MM_over_time <- ggplot(metrics_df, aes(k, silhouette_mean)) +
  geom_line(color="black", linewidth=0.9) +
  geom_point(color="darkgrey", size=4) +
  theme_minimal(base_size=14) +
  theme(
    panel.border = element_rect(color="black", fill=NA, linewidth=0.3),
    panel.background = element_blank(),
    axis.text.x = element_text(size=12, color="black"),
    axis.text.y = element_text(size=12, color="black"),
    axis.title.x = element_text(size=14, color="black"),
    axis.title.y = element_text(size=14, color="black")) +
  labs(x="Number of clusters (k)", y="Average silhouette width")
print(p_sil)

ggplot2::ggsave("Silhouette_vs_k_MM_over_time.svg",
                p_sil_MM_over_time, width = 8, height = 8, units ="cm", dpi = 300)

# ========= 2) PCoA with k-means k=3 + envfit (no labels) =========
# k-means on CLR (k=3)
set.seed(123)
km3 <- kmeans(df_clr, centers=k_fix, nstart=50)
clusters_k3 <- factor(km3$cluster)
perc_SS_k3  <- round(100 * km3$betweenss/km3$totss, 1)

# PCoA on Aitchison distance
pcoa <- cmdscale(dist_aitchison, k=2, eig=TRUE)
pcoa_df <- as.data.frame(pcoa$points)
colnames(pcoa_df) <- c("PCoA1","PCoA2")
pcoa_df$Sample  <- rownames(pcoa_df)
pcoa_df$Cluster <- clusters_k3[match(pcoa_df$Sample, rownames(df_clr))]

# Assign seasons based on "Mon_YYYY"
pcoa_df <- pcoa_df %>%
  mutate(
    Mon3   = str_to_title(str_extract(Sample, "^[A-Za-z]{3}")),
    MonNum = match(Mon3, month.abb),
    Season = case_when(
      MonNum %in% c(3,4,5)   ~ "spring",
      MonNum %in% c(6,7,8)   ~ "summer",
      MonNum %in% c(9,10,11) ~ "autumn",
      MonNum %in% c(12,1,2)  ~ "winter",
      TRUE ~ "other"
    )
  )
season_cols <- c(spring="green3", summer="gold", autumn="orange", winter="dodgerblue", other="grey70")
shape_vals  <- c(21,22,24)  # circle, square, triangle
ve <- round(100 * pcoa$eig[1:2] / sum(pcoa$eig[pcoa$eig > 0]), 1)

# envfit on PCoA (scaled environmental variables)
env_scaled <- scale(metadata)
fit_pcoa  <- envfit(pcoa_df[,c("PCoA1","PCoA2")], env_scaled, permutations=999)
vecs      <- scores(fit_pcoa, display="vectors")
pvals     <- fit_pcoa$vectors$pvals
sig_idx   <- which(pvals < 0.05)

# Data frame of arrows (no labels in plot)
sig_df <- if (length(sig_idx)) {
  rng <- max(diff(range(pcoa_df$PCoA1)), diff(range(pcoa_df$PCoA2)))
  data.frame(
    Variable = rownames(vecs)[sig_idx],
    x = 0, y = 0,
    xend = vecs[sig_idx,1] * 0.35 * rng,
    yend = vecs[sig_idx,2] * 0.35 * rng,
    row.names = NULL
  )
} else data.frame(Variable=character(0), x=numeric(0), y=numeric(0), xend=numeric(0), yend=numeric(0))

# PCoA plot (dark grey arrows, NO variable labels)
p_pcoa <- ggplot(pcoa_df, aes(PCoA1, PCoA2, fill=Season, shape=Cluster)) +
  geom_point(size=point_size, stroke=0.8, color="black") +
  ggrepel::geom_text_repel(aes(label=Sample), size=label_text_size, max.overlaps=Inf,
                           box.padding=0.25, point.padding=0.15, min.segment.length=0) +
  geom_segment(data=sig_df, aes(x=x, y=y, xend=xend, yend=yend),
               inherit.aes=FALSE, arrow=arrow(length=unit(0.28,"cm")),
               color="grey20", alpha=0.95, linewidth=1.2) +
  scale_fill_manual(values=season_cols, name="Season") +
  scale_shape_manual(values=shape_vals, name="Cluster (k = 3)") +
  theme_minimal(base_size=14) +
  theme(
    panel.border = element_rect(color="black", fill=NA, linewidth=0.3),
    axis.text = element_text(color="black", size=axis_text_size),
    axis.title = element_text(size=axis_title_size),
    plot.title = element_text(size=title_size, face="bold"),
    legend.position = "right"
  ) +
  labs(
    title = paste0("PCoA (Aitchison) · k-means k = ", k_fix, "  |  %SS = ", perc_SS_k3),
    x = paste0("PCoA1 (", ve[1], "%)"),
    y = paste0("PCoA2 (", ve[2], "%)")
  )
print(p_pcoa)

# Save as SVG
ggplot2::ggsave(file.path(output_dir, "PCoA_k4_envfit.svg"),
                plot = p_pcoa, width = 14, height = 10, units = "cm")

# Run PERMANOVA for the previously assigned clusters
adonis2(dist_aitchison ~ Cluster, data = pcoa_df)

# Post-hoc test for PERMANOVA
pairwise.adonis(dist_aitchison, pcoa_df$Cluster, p.adjust.m = "fdr")

# ========= 3) Tables with plotted values =========
# (i) Silhouette vs k
silhouette_table <- metrics_df

# (ii) PCoA points with cluster and season
pcoa_table <- pcoa_df %>%
  select(Sample, Cluster, Season, PCoA1, PCoA2)

# (iii) Significant envfit vectors (for reference; not labeled in plot)
envfit_table <- if (nrow(sig_df)) {
  tibble(
    Variable = rownames(vecs)[sig_idx],
    r        = as.numeric(fit_pcoa$vectors$r[sig_idx]),
    p        = as.numeric(pvals[sig_idx]),
    x        = 0, y = 0,
    xend     = sig_df$xend,
    yend     = sig_df$yend
  )
} else {
  tibble(Variable=character(0), r=numeric(0), p=numeric(0),
         x=numeric(0), y=numeric(0), xend=numeric(0), yend=numeric(0))
}

# --- Export to Excel ---
wb <- createWorkbook()
addWorksheet(wb, "Silhouette_vs_k")
writeData(wb, "Silhouette_vs_k", silhouette_table)

addWorksheet(wb, "PCoA_points")
writeData(wb, "PCoA_points", pcoa_table)

addWorksheet(wb, "Envfit_vectors_sig")
writeData(wb, "Envfit_vectors_sig", envfit_table)

# --- Save tables to target path ---
output_dir <- "PARTITION_PLOTS"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

xlsx_path <- file.path(output_dir, "MM_temporal_silhouette_pcoa_envfit.xlsx")
saveWorkbook(wb, xlsx_path, overwrite = TRUE)

message("Tables saved to: ", xlsx_path)
