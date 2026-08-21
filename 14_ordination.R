#-------------------------------------------------------------------#
#--- Stage 14: Beta-diversity ordination (PCoA)                   ---#
#-------------------------------------------------------------------#
#
# Purpose: Run PCoA ordination (Bray-Curtis) for each metabarcoding
#          marker and compute the convex hulls used to group samples
#          by site in the ordination figures. Produces the
#          plotting-ready data consumed by 19_figures.R 
#
# Depends on: 13_phyloseq_objects.R (data/processed/physeq_list.rds)
#
# Output: data/processed/ordination_data.rds
#
library(phyloseq)
library(dplyr)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

physeq_list <- readRDS("data/processed/physeq_list.rds")


compute_ordination <- function(ps) {
  ord <- ordinate(ps, method = "PCoA", distance = "bray")

  # Stable extraction method (no vegan::scores)
  plot_data <- plot_ordination(ps, ord, type = "samples")$data
  plot_data$Distance <- factor(plot_data$Distance, levels = c(0, 50, 100))

  if (!"site" %in% colnames(plot_data)) {
    warning("`site` column not found in ordination plot data -- check 12_metadata_processing.R output.")
  }

  hulls <- plot_data %>%
    dplyr::group_by(site) %>%
    dplyr::slice(chull(Axis.1, Axis.2))

  list(ordination = ord, plot_data = plot_data, hulls = hulls)
}

ordination_data <- lapply(physeq_list, compute_ordination)

saveRDS(ordination_data, "data/processed/ordination_data.rds")
