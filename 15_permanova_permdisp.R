#-------------------------------------------------------------------#
#--- Stage 15: PERMANOVA and PERMDISP                              ---#
#-------------------------------------------------------------------#
#
# Purpose: Pairwise PERMANOVA (adonis2) between sites for each marker
#
# Depends on: 13_phyloseq_objects.R (data/processed/physeq_list.rds)
#
# Output: data/processed/permanova_results.rds
#         data/processed/permdisp_results.rds   


library(phyloseq)
library(vegan)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

physeq_list <- readRDS("data/processed/physeq_list.rds")

combined_pairwise_results <- data.frame()

for (dataset_name in names(physeq_list)) {

  phy_obj <- physeq_list[[dataset_name]]
  dist_mat <- phyloseq::distance(phy_obj, method = "bray")
  meta <- data.frame(sample_data(phy_obj))
  if (!"site" %in% colnames(meta)) {
    warning(paste("Skipping", dataset_name, "- 'site' column not found"))
    next
  }

  sites <- unique(meta$site)

  for (i in 1:(length(sites) - 1)) {
    for (j in (i + 1):length(sites)) {

      subset_meta <- meta[meta$site %in% c(sites[i], sites[j]), ]

      subset_dist <- as.dist(
        as.matrix(dist_mat)[rownames(subset_meta), rownames(subset_meta)]
      )

      test <- adonis2(subset_dist ~ site, data = subset_meta)

      combined_pairwise_results <- rbind(
        combined_pairwise_results,
        data.frame(
          Dataset = dataset_name,
          Group1 = sites[i],
          Group2 = sites[j],
          F = test$F[1],
          R2 = test$R2[1],
          p = test$`Pr(>F)`[1]
        )
      )
    }
  }
}

combined_pairwise_results$R2 <- round(combined_pairwise_results$R2, 3)
combined_pairwise_results$p  <- signif(combined_pairwise_results$p, 3)

saveRDS(combined_pairwise_results, "data/processed/permanova_results.rds")

# ---- PERMDISP 
combined_permdisp_results <- data.frame()

for (dataset_name in names(physeq_list)) {

  phy_obj <- physeq_list[[dataset_name]]
  dist_mat <- phyloseq::distance(phy_obj, method = "bray")
  meta <- data.frame(sample_data(phy_obj))
  if (!"site" %in% colnames(meta)) next

  sites <- unique(meta$site)

  for (i in 1:(length(sites) - 1)) {
    for (j in (i + 1):length(sites)) {

      subset_meta <- meta[meta$site %in% c(sites[i], sites[j]), ]
      subset_dist <- as.dist(
        as.matrix(dist_mat)[rownames(subset_meta), rownames(subset_meta)]
      )

      disp <- betadisper(subset_dist, subset_meta$site)
      disp_test <- permutest(disp)

      combined_permdisp_results <- rbind(
        combined_permdisp_results,
        data.frame(
          Dataset = dataset_name,
          Group1 = sites[i],
          Group2 = sites[j],
          F = disp_test$tab$F[1],
          p = disp_test$tab$`Pr(>F)`[1]
        )
      )
    }
  }
}

combined_permdisp_results$p <- signif(combined_permdisp_results$p, 3)

saveRDS(combined_permdisp_results, "data/processed/permdisp_results.rds")
