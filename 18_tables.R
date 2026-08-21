#-------------------------------------------------------------------#
#--- Stage 18: Export result tables                                ---#
#-------------------------------------------------------------------#
#
# Purpose: Write the final CSV tables reported alongside the
#          manuscript, gathering the results computed by
#          15_permanova_permdisp.R and 17_ddPCR_analysis.R in one
#          place. 
#
# Depends on: 15_permanova_permdisp.R, 17_ddPCR_analysis.R
#
# Output: outputs/tables/positive_detections.csv
#         outputs/tables/lod_loq_summary.csv
#         outputs/tables/pairwise_PERMANOVA_results_all.csv
#         outputs/tables/pairwise_PERMANOVA_results_summary.csv
#         outputs/tables/pairwise_PERMDISP_results.csv     
#           15_permanova_permdisp.R 
#         outputs/tables/<marker>_indic_species_kelp_regions.csv
#         

library(dplyr)

dir.create("outputs/tables", recursive = TRUE, showWarnings = FALSE)

# ---- ddPCR tables (from 17_ddPCR_analysis.R) ----
positive_only <- readRDS("data/processed/positive_detections.rds")
write.csv(positive_only, "outputs/tables/positive_detections.csv", row.names = FALSE)

lod_loq_summary <- readRDS("data/processed/lod_loq_summary.rds")
write.csv(lod_loq_summary, "outputs/tables/lod_loq_summary.csv", row.names = FALSE)

# ---- PERMANOVA / PERMDISP tables (from 15_permanova_permdisp.R) ----
combined_pairwise_results <- readRDS("data/processed/permanova_results.rds")
write.csv(
  combined_pairwise_results,
  "outputs/tables/pairwise_PERMANOVA_results_all.csv",
  row.names = FALSE
)

summary_table <- combined_pairwise_results %>%
  group_by(Dataset) %>%
  summarise(
    mean_R2 = round(mean(R2, na.rm = TRUE), 3),
    min_p = signif(min(p, na.rm = TRUE), 3),
    sig_pairs = sum(p < 0.05, na.rm = TRUE),
    total_pairs = n(),
    perc_sig = round(100 * sig_pairs / total_pairs, 1)
  )

write.csv(
  summary_table,
  "outputs/tables/pairwise_PERMANOVA_results_summary.csv",
  row.names = FALSE
)

permdisp_results <- readRDS("data/processed/permdisp_results.rds")
write.csv(
  permdisp_results,
  "outputs/tables/pairwise_PERMDISP_results.csv",
  row.names = FALSE
)

# ---- Indicator species tables (from 16_indicator_species.R) ----
# One CSV per marker, matching your original script's per-marker
# write.csv() calls (rbcL_indic_species_kelp_regions.csv etc.).
indicator_results_path <- "data/processed/indicator_species_results.rds"
if (file.exists(indicator_results_path)) {
  indicator_results <- readRDS(indicator_results_path)

  file_prefixes <- c(
    rbcL = "rbcL", COI = "COI", V7_18S = "v7", V9_18S = "v9", Bacterial = "bacterial"
  )

  for (dataset_name in names(indicator_results)) {
    overlap <- indicator_results[[dataset_name]]$tory_outerqc_overlap
    if (is.null(overlap) || nrow(overlap) == 0) next

    prefix <- file_prefixes[[dataset_name]]
    if (is.null(prefix)) prefix <- tolower(dataset_name)

    write.csv(
      overlap,
      paste0("outputs/tables/", prefix, "_indic_species_kelp_regions.csv"),
      row.names = FALSE
    )
  }
}
