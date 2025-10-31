# ==========================================================================================
# Fig1 -- mean-abundance and alpha-diversity metrics -Mar Menor versus open Mediterranean Sea
# ==========================================================================================

# ---- LIBRARIES ----
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(vegan)
  library(tibble)
  library(cowplot)  # align_plots, plot_grid
  library(grid)     # unit(), grid.newpage, grid.draw
  library(gtable)
})

# ---- INPUT AND OUTPUT PATHS ----
input_combined <- "Med-MM_vOTU_abund_table.csv"
total_nts_file_medit <- "total_nts_in_reads_per_metagenome_Medit.csv"
total_nts_file_mm <- "total_nts_in_reads_per_metagenome_MM_time_series.csv"

output_dir <- "Mean-abundance_Alpha-diversity_metrics_MM_vs_open-Med"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ---- TEXT PARAMETERS ----
axis_text_size    <- 9   # axis tick labels
y_axis_title_size <- 9   # Y-axis title
x_angle <- 45
x_size  <- 9

# ---- PANEL PARAMETERS (PLOTTING AREA) ----
# Use these two values to unify the size of the “black square”
panel_w_cm <- 9.0
panel_h_cm <- 3.2   # reduce if you want a more compact column

# ---- THEME (no plot title) with minimal vertical margins ----
my_theme <- function(show_x = FALSE) {
  theme_minimal(base_size = 9) +
    theme(
      panel.border   = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      legend.position = "none",
      axis.title.x   = element_blank(),
      axis.title.y   = element_text(size = y_axis_title_size, colour = "black"),
      axis.text.y    = element_text(size = axis_text_size,  colour = "black"),
      axis.text.x    = if (show_x) element_text(angle = x_angle, hjust = 1, size = x_size, colour = "black") else element_blank(),
      axis.ticks.x   = if (show_x) element_line() else element_blank(),
      # zero outer vertical margins to allow tight stacking of panels
      plot.margin    = margin(t = 0, r = 5.5, b = 0, l = 5.5)
    )
}

# ---- BOXPLOT (MM vs MED) ----
plot_box <- function(metric, ylab, show_x = TRUE, diversity_metrics) {
  ggplot(diversity_metrics, aes(x = Group, y = .data[[metric]], fill = Group)) +
    geom_boxplot(width = 0.45, outlier.shape = NA) +
    geom_point(
      aes(color = Group, size = Group),   # different point sizes by group
      position = position_jitter(width = 0.15),
      alpha = 0.65,
      shape = 16                          # solid circle
    ) +
    scale_size_manual(values = c("Mar_Menor" = 0.8, "Medit_Sea" = 0.8)) +
    scale_fill_manual(values = c("Mar_Menor" = "lightgreen", "Medit_Sea" = "lightblue")) +
    scale_color_manual(values = c("Mar_Menor" = "darkgreen", "Medit_Sea" = "steelblue")) +
    scale_x_discrete(limits = c("Mar_Menor", "Medit_Sea"),
                     expand = expansion(mult = c(0.4, 0.4))) +
    labs(y = ylab, x = NULL) +
    my_theme(show_x)
}

# ---- DATA LOADING ----
abund_all <- read_delim(input_combined, delim = ";", locale = locale(decimal_mark = ","))

total_nts_medit <- read_delim(total_nts_file_medit, delim = ";", locale = locale(decimal_mark = ",")) %>%
  rename(Total_nts = `Total number of metagenome read nts`)

total_nts_mm <- read_delim(total_nts_file_mm, delim = ";", locale = locale(decimal_mark = ",")) %>%
  rename(Total_nts = `Total number of metagenome read nts`)

total_nts <- bind_rows(total_nts_medit, total_nts_mm)

# ---- LONG/WIDE FORMAT AND METRIC CALCULATIONS ----
sample_order_all <- colnames(abund_all)[-1]

abund_long <- abund_all %>%
  pivot_longer(-1, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Sample = factor(Sample, levels = sample_order_all)) %>%
  mutate(Group = case_when(
    str_starts(Sample, "TARA") ~ "Medit_Sea",
    str_starts(Sample, "Med")  ~ "Medit_Sea",
    str_detect(Sample, "^[A-Za-z]{3}_\\d{4}$") ~ "Mar_Menor",
    TRUE ~ "Unknown"
  )) %>%
  filter(Group %in% c("Mar_Menor", "Medit_Sea"))

abund_wide <- abund_long %>%
  select(-Group) %>%
  pivot_wider(names_from = Sample, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames(var = names(abund_all)[1]) %>%
  t() %>%
  as.data.frame()

diversity_metrics <- data.frame(
  Sample = rownames(abund_wide),
  Richness = rowSums(abund_wide > 0),
  Shannon = diversity(abund_wide, index = "shannon"),
  Mean_abundance_all = rowMeans(abund_wide),
  Pielou = NA_real_
)

diversity_metrics$Pielou <- with(diversity_metrics,
                                 ifelse(Richness > 1, Shannon / log(Richness), NA))

diversity_metrics <- diversity_metrics %>%
  left_join(total_nts, by = "Sample") %>%
  mutate(
    Group = case_when(
      str_starts(Sample, "TARA") ~ "Medit_Sea",
      str_starts(Sample, "Med")  ~ "Medit_Sea",
      str_detect(Sample, "^[A-Za-z]{3}_\\d{4}$") ~ "Mar_Menor",
      TRUE ~ "Unknown"
    ),
    Group = factor(Group, levels = c("Mar_Menor", "Medit_Sea")),
    Richness_per_Gb = (Richness * 1e9) / Total_nts
  ) %>%
  filter(Group %in% c("Mar_Menor", "Medit_Sea"))

# ---- CREATE THE 5 BOXPLOTS (ggplot) ----
p_mean   <- plot_box("Mean_abundance_all", "Mean abundance",  show_x = TRUE, diversity_metrics)
p_rich   <- plot_box("Richness",           "Richness",         show_x = TRUE, diversity_metrics)
p_richGb <- plot_box("Richness_per_Gb",    "Richness/Gb",      show_x = TRUE, diversity_metrics)
p_shan   <- plot_box("Shannon",            "Shannon",          show_x = TRUE, diversity_metrics)
p_even   <- plot_box("Pielou",             "Evenness",         show_x = TRUE, diversity_metrics)

p_mean
p_rich
p_richGb
p_shan
p_even
