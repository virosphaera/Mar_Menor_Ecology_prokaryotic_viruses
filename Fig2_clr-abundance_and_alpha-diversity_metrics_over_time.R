# --------------------------------------------------------
# CLR ANALYSIS SCRIPT FOR vOTUs IN PICOPLANKTON METAGENOMES
# CREATES AND SAVES: CLR MATRIX + CLR BOXPLOT + BARPLOTS (SVG)
# --------------------------------------------------------

# Load required libraries
suppressPackageStartupMessages({
  library(nortest)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(compositions)
  library(rstatix)
  library(tibble)
  library(vegan)
  library(Hmisc)
  library(pheatmap)
})

# --------------------------------------------------------
# TEXT PARAMETERS AND COMMON THEME
# --------------------------------------------------------
axis_text_size    <- 9    # axis tick labels (X and Y)
y_axis_title_size <- 9    # Y-axis title
x_angle <- 45
x_size  <- 7

my_theme <- function(show_x = TRUE) {
  theme_minimal(base_size = 9) +
    theme(
      panel.border   = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      legend.position = "none",
      axis.title.x   = element_blank(),
      axis.title.y   = element_text(size = y_axis_title_size, colour = "black"),
      axis.text.y    = element_text(size = axis_text_size, colour = "black"),
      axis.text.x    = if (show_x) element_text(angle = x_angle, hjust = 1, size = x_size, colour = "black") else element_blank(),
      axis.ticks.x   = if (show_x) element_line() else element_blank(),
      plot.margin    = margin(t = 0, r = 5.5, b = 0, l = 5.5)
    )
}

# --------------------------------------------------------
# PATHS
# --------------------------------------------------------
input_votu <- "vOTU_relative_abundances_MM_time_series.csv"
ts_file    <- "total_nts_in_reads_per_metagenome_MM_time_series.csv"
out_dir    <- "clr-abundance_alpha-diversity_over_time"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# --------------------------------------------------------
# DATA LOADING AND CLR TRANSFORMATION
# --------------------------------------------------------
vOTU <- read.table(input_votu, sep = ";", dec = ",", header = TRUE, row.names = 1)

# Clean potential whitespace in sample names
colnames(vOTU) <- trimws(colnames(vOTU))

cat("✅ Total sum per sample:\n"); print(colSums(vOTU))
cat("\n🔍 Maximum value:\n"); print(max(vOTU))

# CLR transformation with pseudo-count
df_pseudo <- vOTU + 0.001
df_clr    <- as.data.frame(clr(df_pseudo))

# Ensure clean and aligned column names
colnames(df_clr) <- trimws(colnames(df_clr))

clr_medians <- apply(df_clr, 2, median)

# Save CLR matrix (vOTU ID as first column)
df_clr_out <- df_clr %>% rownames_to_column("vOTU")
write.table(df_clr_out,
            file = file.path(out_dir, "vOTU_CLR_matrix.csv"),
            sep = ";", dec = ",", row.names = FALSE)

# Convert to long format for boxplot (without forcing incorrect factor levels)
df_long <- df_clr_out %>%
  pivot_longer(cols = -vOTU, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Sample = trimws(Sample))

# --------------------------------------------------------
# CLR ABUNDANCE BOXPLOT (object + export)
# --------------------------------------------------------
# Fixed sample order
lvl <- colnames(df_clr)

# Median per sample
meds <- df_long |>
  dplyr::mutate(Sample = factor(Sample, levels = lvl)) |>
  dplyr::group_by(Sample) |>
  dplyr::summarise(y = median(Abundance, na.rm = TRUE), .groups = "drop") |>
  dplyr::mutate(x = as.numeric(Sample))

p_clr <- ggplot(df_long, aes(x = factor(Sample, levels = lvl), y = Abundance)) +
  geom_boxplot(
    fill = "lightgreen", color = "black",
    linewidth = 0.3,
    fatten = 0,           # hide default median line
    outlier.shape = NA,
    na.rm = TRUE
  ) +
  # Median as a thin filled rectangle (in data units)
  geom_rect(
    data = meds,
    aes(xmin = x - 0.35, xmax = x + 0.35,
        ymin = y - 0.1, ymax = y + 0.1),   # controls “thickness” visually
    inherit.aes = FALSE,
    fill = "black", color = NA
  ) +
  labs(x = "Sampling date", y = "CLR-Abundance") +
  scale_x_discrete(limits = lvl) +
  coord_cartesian(ylim = c(-7, 10)) +
  my_theme(show_x = TRUE)

p_clr

# --------------------------------------------------------
# ECOLOGICAL METRICS AND SUMMARY TABLE
# --------------------------------------------------------
shannon_index   <- apply(vOTU, 2, function(x) diversity(x, index = "shannon"))
richness_obs    <- apply(vOTU, 2, function(x) sum(x > 0))
pielou_evenness <- shannon_index / log(richness_obs)

ts_data <- read.table(ts_file, sep = ";", dec = ",", header = TRUE, row.names = 1)
ts_data <- ts_data[colnames(vOTU), , drop = FALSE]
vOTUs_per_Gb <- (richness_obs / ts_data[,1]) * 1e9

div_table <- data.frame(
  Sample = factor(colnames(vOTU), levels = colnames(vOTU)),
  CLR_Median = clr_medians,
  Observed_Richness = richness_obs,
  vOTUs_per_Gb = vOTUs_per_Gb,
  Shannon = shannon_index,
  Pielou = pielou_evenness
)

write.table(div_table, file = file.path(out_dir, "diversity_stats_CLR.csv"),
            sep = ";", dec = ",", row.names = FALSE)

# --------------------------------------------------------
# LIMITS AND BARPLOTS AS SEPARATE OBJECTS
# --------------------------------------------------------
y_limits <- list(
  Observed_Richness = c(2000, 8500),
  vOTUs_per_Gb      = c(300, 750),
  Shannon           = c(7, 8.5),
  Pielou            = c(0.90, 0.95)
)

p_richness <- ggplot(div_table, aes(x = Sample, y = Observed_Richness)) +
  geom_col(fill = "lightgreen", color = "black", linewidth = 0.3) +
  coord_cartesian(ylim = y_limits[["Observed_Richness"]]) +
  labs(title = NULL, x = NULL, y = "Richness") +
  my_theme(show_x = TRUE)

p_richGb <- ggplot(div_table, aes(x = Sample, y = vOTUs_per_Gb)) +
  geom_col(fill = "lightgreen", color = "black", linewidth = 0.3) +
  coord_cartesian(ylim = y_limits[["vOTUs_per_Gb"]]) +
  labs(title = NULL, x = NULL, y = "Richness/Gb") +
  my_theme(show_x = TRUE)

p_shannon <- ggplot(div_table, aes(x = Sample, y = Shannon)) +
  geom_col(fill = "lightgreen", color = "black", linewidth = 0.3) +
  coord_cartesian(ylim = y_limits[["Shannon"]]) +
  labs(title = NULL, x = NULL, y = "Shannon index") +
  my_theme(show_x = TRUE)

p_pielou <- ggplot(div_table, aes(x = Sample, y = Pielou)) +
  geom_col(fill = "lightgreen", color = "black", linewidth = 0.3) +
  coord_cartesian(ylim = y_limits[["Pielou"]]) +
  labs(title = NULL, x = NULL, y = "Evenness") +
  my_theme(show_x = TRUE)

# (Optional) view in RStudio
print(p_clr); print(p_richness); print(p_richGb); print(p_shannon); print(p_pielou)
