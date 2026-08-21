#-------------------------------------------------------------------#
#--- Stage 19: Figures                                             ---#
#-------------------------------------------------------------------#
#
# Purpose: Render every publication figure from the processed data
#          produced by the earlier pipeline stages.
#
# Depends on: 13_phyloseq_objects.R, 14_ordination.R, 17_ddPCR_analysis.R
#

# Output: outputs/figures/LOD_data_ddPCR.tiff / .pdf
#         outputs/figures/Figure3.tiff / .pdf
#         outputs/figures/ordination_overview.tiff / .pdf
#         outputs/figures/<marker>_ordination_hull.tiff / .pdf
#         outputs/figures/Laminariaceae_comparison.tiff / .pdf
#         outputs/figures/<marker>_indicator_species_venn.tiff / .pdf
#         outputs/figures/Figure 6 <marker>.tiff / .pdf
#

library(ggplot2)
library(dplyr)
library(patchwork)
library(phyloseq)
library(ggVennDiagram)

dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

# Saves every figure as BOTH a .tiff and a .pdf.
save_figure <- function(filename, plot, ...) {
  ggsave(filename, plot, ...)

  pdf_filename <- sub("\\.tiff$", ".pdf", filename, ignore.case = TRUE)
  if (identical(pdf_filename, filename)) {
    pdf_filename <- paste0(filename, ".pdf")
  }
  ggsave(pdf_filename, plot, ...)
}


recode_site <- function(x) {
  lookup <- c(
    "Tory" = "TC", "TC" = "TC",
    "Pelorus" = "PS", "PS" = "PS",
    "InnerQC" = "Inner QCS", "Inner QC" = "Inner QCS",
    "OuterQC" = "Outer QCS", "Outer QC" = "Outer QCS"
  )
  x_chr <- as.character(x)
  out <- unname(lookup[x_chr])
  
  out[is.na(out)] <- x_chr[is.na(out)]
  out
}

#=====================================================================
# Figure 2 -- ddPCR LOD/LOQ
#=====================================================================

lod_loq_summary <- readRDS("data/processed/lod_loq_summary.rds")
base_y <- 1e-3

LOD_data <- ggplot(lod_loq_summary, aes(x = Inputconc.)) +
  geom_rect(aes(
    xmin = as.numeric(factor(Inputconc.)) - 0.4,
    xmax = as.numeric(factor(Inputconc.)) + 0.4,
    ymin = base_y,
    ymax = Mean1),
    fill = "steelblue") +
  scale_y_log10(limits = c(base_y, NA)) +
  labs(y = "Copies/μL", x = "Input concentration") +
  theme_minimal() +
  theme(
    text = element_text(family = "sans", color = "black", size = 8),
    axis.title.y = element_text(size = 8),
    axis.text = element_text(size = 8, colour = "black")
  ) +
  geom_errorbar(aes(
    x = as.numeric(factor(Inputconc.)),
    ymin = Mean1 - SE,
    ymax = Mean1 + SE),
    width = 0.2,
    color = "black")

save_figure("outputs/figures/LOD_data_ddPCR.tiff", plot = LOD_data)

#=====================================================================
# Figure 3 -- ddPCR field detections
#=====================================================================

positive_only <- readRDS("data/processed/positive_detections.rds")

plot_theme <- theme_minimal() +
  theme(
    text = element_text(family = "sans", colour = "black", size = 12),
    axis.title.y = element_text(size = 12),
    axis.text = element_text(size = 12, colour = "black")
  )

colour_palette <- c("#201321ff", "#3485a5ff", "#81d8b0ff")

DISTANCE_LEVELS <- c("0", "50", "100")
Y_MAX <- max(positive_only$Conc.copies.µL., na.rm = TRUE) * 1.05

create_location_plot <- function(df, location_name, title) {
  plot_data <- df %>%
    filter(Location == location_name) %>%
    mutate(Distance = factor(Distance, levels = DISTANCE_LEVELS))

  ggplot(plot_data, aes(x = Distance, y = Conc.copies.µL.)) +
    geom_jitter(aes(colour = Transect, shape = Replicate), width = 0.2, alpha = 0.6, size = 2) +
    geom_hline(yintercept = 2.6, linetype = "dotted") +
    scale_x_discrete(limits = DISTANCE_LEVELS) +
    scale_color_manual(values = colour_palette) +
    coord_cartesian(ylim = c(0, Y_MAX)) +
    labs(x = NULL, y = "copies/µL", title = title) +
    plot_theme +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
}


tc_plot       <- create_location_plot(positive_only, "TC", recode_site("TC"))
outer_qc_plot <- create_location_plot(positive_only, "Outer QC", recode_site("Outer QC"))
inner_qc_plot <- create_location_plot(positive_only, "Inner QC", recode_site("Inner QC"))
pelorus_plot  <- create_location_plot(positive_only, "PS", recode_site("PS"))

fig3 <- (tc_plot | outer_qc_plot) / (inner_qc_plot | pelorus_plot)

save_figure("outputs/figures/Figure3.tiff", fig3)

#=====================================================================
# Ordination figures -- v1: combined multi-panel
#=====================================================================

ordination_data <- readRDS("data/processed/ordination_data.rds")

ordination_data <- lapply(ordination_data, function(od) {
  od$plot_data$site <- recode_site(od$plot_data$site)
  od$hulls$site <- recode_site(od$hulls$site)
  od
})

site_cols <- c(
  "Inner QCS" = "#3c3162ff",
  "Outer QCS" = "#3670a0ff",
  "PS" = "#3cb2adff",
  "TC" = "#b6e5c4ff"
)

marker_titles <- c(
  Bacterial = "Bacterial", COI = "COI", rbcL = "rbcL",
  V7_18S = "V7 18S", V9_18S = "V9 18S"
)

make_hull_plot_v1 <- function(od, title) {
  ggplot(od$plot_data, aes(Axis.1, Axis.2)) +
    geom_polygon(data = od$hulls, aes(fill = site, group = site), alpha = 0.2, colour = NA) +
    geom_point(aes(colour = site, shape = Distance), size = 3, stroke = 0.6) +
    scale_shape_manual(values = c("0" = 16, "50" = 17, "100" = 15)) +
    scale_fill_manual(values = site_cols) +
    scale_colour_manual(values = site_cols) +
    labs(title = title, x = "Axis 1", y = "Axis 2") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 12),
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = 12),
      axis.title = element_text(size = 12)
    )
}

hull_plots_v1 <- Map(make_hull_plot_v1, ordination_data, marker_titles[names(ordination_data)])

ord_plot <- (hull_plots_v1$Bacterial | hull_plots_v1$V7_18S | hull_plots_v1$V9_18S) /
  (hull_plots_v1$COI | hull_plots_v1$rbcL)

save_figure("outputs/figures/ordination_overview.tiff", ord_plot)

#=====================================================================
# Ordination figures -- v2
#=====================================================================

ordination_theme <- theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 12),
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = 12),
    axis.title = element_text(size = 12),
    strip.text = element_text(colour = "black", size = 12)
  )

make_hull_plot_v2 <- function(od, title) {
  ggplot(od$plot_data, aes(Axis.1, Axis.2)) +
    geom_polygon(data = od$hulls, aes(fill = site, group = site), alpha = 0.2) +
    geom_point(aes(color = site, shape = Distance), size = 3) +
    scale_fill_manual(values = site_cols) +
    scale_color_manual(values = site_cols) +
    labs(title = title, x = "Axis 1", y = "Axis 2") +
    ordination_theme
}

hull_plots_v2 <- Map(make_hull_plot_v2, ordination_data, marker_titles[names(ordination_data)])

save_figure("outputs/figures/bacterial_ordination_hull.tiff", hull_plots_v2$Bacterial)
save_figure("outputs/figures/V7_ordination_hull.tiff", hull_plots_v2$V7_18S)
save_figure("outputs/figures/V9_ordination_hull.tiff", hull_plots_v2$V9_18S)
save_figure("outputs/figures/COI_ordination_hull.tiff", hull_plots_v2$COI)
save_figure("outputs/figures/rbcL_ordination_hull.tiff", hull_plots_v2$rbcL)


#=====================================================================
# Laminariaceae metabarcoding-vs-ddPCR comparison
#=====================================================================

laminariaceae_comparison <- readRDS("data/processed/laminariaceae_comparison.rds")

# Standardize the `Location` labels used as facet titles below (see
# recode_site() near the top of this script).
laminariaceae_comparison$Location <- recode_site(laminariaceae_comparison$Location)

LAMINARIACEAE_GENE_COLOURS <- setNames(colour_palette, c("COI", "18S-V7", "18S-V9"))

laminariaceae_plot <- ggplot(
  laminariaceae_comparison,
  aes(x = Gene, y = Number_Laminariaceae, fill = Gene)
) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_y_log10(labels = scales::scientific) +
  scale_fill_manual(values = LAMINARIACEAE_GENE_COLOURS) +
  facet_wrap(~Location) +
  labs(title = NULL, x = NULL, y = "Abundance of Laminariaceae OTUs") +
  theme_minimal(base_size = 12) +
  theme(
    text = element_text(color = "black"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 12),
    axis.text.y = element_text(angle = 0, hjust = 1, colour = "black", size = 12),
    strip.text = element_text(colour = "black", size = 12),
    axis.title.y = element_text(size = 12),
    legend.position = "none"
  )

save_figure("outputs/figures/Laminariaceae_comparison.tiff", laminariaceae_plot)

#=====================================================================
# Figure 6 -- indicator species class composition
#=====================================================================
# (No location/site labels appear in this figure -- panels are split by
# marker, not by site -- so recode_site() doesn't apply here.)
FIG6_TITLES <- c(
  Bacterial = "Bacterial", V7_18S = "18S-V7", V9_18S = "18S-V9",
  rbcL = "rbcL", COI = "COI"
)
FIG6_ITALIC <- c(
  Bacterial = FALSE, V7_18S = FALSE, V9_18S = FALSE, rbcL = TRUE, COI = FALSE
)
FIG6_FILENAMES <- c(
  Bacterial = "bacterial", V7_18S = "V7", V9_18S = "V9", rbcL = "rbcL", COI = "COI"
)

make_indicator_class_barplot <- function(overlap_df, title, italic_title = FALSE) {
  filtered <- overlap_df %>% filter(stat > 0.8)

  plot_data <- filtered %>%
    count(class) %>%
    mutate(prop = n / sum(n))

  # Relabel NA class as "Unknown"
  class_chr <- as.character(plot_data$class)
  class_chr[is.na(class_chr)] <- "Unknown"
  plot_data$class <- factor(class_chr)

  ggplot(plot_data, aes(x = "", y = prop, fill = class)) +
    geom_bar(stat = "identity", width = 1) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(
      x = NULL, y = paste0("Proportion of taxa (n=", nrow(filtered), ")"),
      title = title, fill = "Class"
    ) +
    scale_fill_viridis_d(option = "mako") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      text = element_text(color = "black"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      plot.title = element_text(color = "black", hjust = 0.5, face = if (italic_title) "italic" else "plain"),
      legend.title = element_text(color = "black"),
      legend.text = element_text(color = "black")
    )
}

for (dataset_name in names(indicator_results)) {
  overlap <- indicator_results[[dataset_name]]$tory_outerqc_overlap

  if (is.null(overlap) || !"stat" %in% colnames(overlap) ||
      nrow(overlap %>% filter(stat > 0.8)) == 0) {
    message(
      dataset_name,
      ": no ASVs with indicator value (stat) > 0.8 in the Tory/Outer QC ",
      "overlap -- skipping Figure 6 panel."
    )
    next
  }

  fig6_plot <- make_indicator_class_barplot(
    overlap, FIG6_TITLES[[dataset_name]], FIG6_ITALIC[[dataset_name]]
  )

  save_figure(
    paste0("outputs/figures/Figure 6 ", FIG6_FILENAMES[[dataset_name]], ".tiff"),
    fig6_plot
  )
}
