#-------------------------------------------------------------------#
#--- Stage 16: Indicator species analysis                         ---#
#-------------------------------------------------------------------#
#
# Purpose: For each marker, runs indicspecies::multipatt()
#          across all four sites, then uses a four-set Venn diagram to
#          pull out the ASVs shared between Tory Channel and Outer QC
#          (the two M. pyrifera-present sites). These are the
#          "indicator species associated with kelp regions" results
#          reported in the manuscript. Plotting (the Venn diagrams
#          themselves) happens in 19_figures.R; CSV export happens in
#          18_tables.R -- this script only computes and saves results,
#          matching the pattern used by 14/15/17.
#
# Depends on: 13_phyloseq_objects.R (data/processed/physeq_list.rds)
#
# Output: data/processed/indicator_species_results.rds

library(phyloseq)
library(indicspecies)
library(vegan)
library(dplyr)
library(tibble)
library(ggVennDiagram)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

physeq_list <- readRDS("data/processed/physeq_list.rds")

# Markers whose ASV/OTU table should have Laminariaceae-family taxa
# removed before the indicator species analysis 
LAMINARIACEAE_EXCLUDE <- c("COI", "V7_18S", "V9_18S")

SITE_ORDER <- c("Tory", "Pelorus", "InnerQC", "OuterQC")
SITE_DISPLAY_NAMES <- c("Tory Channel", "Pelorus", "Inner QC", "Outer QC")

run_indicator_species <- function(dataset_name, ps) {

  if (dataset_name %in% LAMINARIACEAE_EXCLUDE) {
    ps <- subset_taxa(ps, family != "Laminariaceae")
  }

  otu <- as.data.frame(otu_table(ps))
  if (taxa_are_rows(ps)) otu <- as.data.frame(t(otu))
  group <- sample_data(ps)$site

  res <- multipatt(otu, group, func = "IndVal", duleg = FALSE, control = how(nperm = 999))

  signif_ids <- rownames(res$sign[res$sign$p.value <= 0.05, ])

  tax_df <- as.data.frame(tax_table(ps))
  sig_results <- cbind(tax_df[signif_ids, , drop = FALSE], res$sign[signif_ids, , drop = FALSE])

  if ("superkingdom" %in% colnames(sig_results)) {
    sig_results <- sig_results %>% filter(!is.na(superkingdom))  # see FLAG #4 above
  }

  # ---- Per-site significant-ASV lists (input to the Venn diagram) ----
  venn_lists <- lapply(SITE_ORDER, function(s) {
    col <- paste0("s.", s)
    if (!col %in% colnames(sig_results)) {
      warning(
        dataset_name, ": expected multipatt column '", col, "' not found -- ",
        "check that SITE_ORDER above matches the actual `site` values for ",
        "this marker. Available columns: ",
        paste(grep("^s\\.", colnames(sig_results), value = TRUE), collapse = ", ")
      )
      return(character(0))
    }
    sig_results %>%
      filter(.data[[col]] == 1, index != 0) %>%
      pull(asv)
  })
  names(venn_lists) <- SITE_DISPLAY_NAMES

  # ---- Tory vs Outer QC overlap  ----
  overlap <- NULL
  if (all(lengths(venn_lists) > 0 | seq_along(venn_lists) %in% c(1, 4))) {
    venn_obj <- ggVennDiagram::Venn(venn_lists)
    vd <- ggVennDiagram::process_data(venn_obj)
    region <- vd$regionData %>% filter(id == "1/4")  # set 1 = Tory Channel, set 4 = Outer QC
    asvs <- unlist(region$item)
    asvs <- asvs[asvs %in% taxa_names(ps)]

    res_df <- as.data.frame(res$sign) %>% tibble::rownames_to_column("ASV")

    overlap <- as.data.frame(tax_table(ps))[asvs, , drop = FALSE] %>%
      tibble::rownames_to_column("ASV") %>%
      left_join(res_df, by = "ASV")
  }

  list(
    multipatt = res,
    significant = sig_results,
    venn_lists = venn_lists,
    tory_outerqc_overlap = overlap
  )
}

indicator_results <- Map(run_indicator_species, names(physeq_list), physeq_list)

saveRDS(indicator_results, "data/processed/indicator_species_results.rds")

message(
  "16_indicator_species.R: done. Significant-taxon counts (p<=0.05) per marker: ",
  paste(names(indicator_results), vapply(indicator_results, function(x) nrow(x$significant), integer(1)), sep = "=", collapse = ", ")
)
