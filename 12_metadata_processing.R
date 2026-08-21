#-------------------------------------------------------------------#
#--- Stage 12: Metadata processing                                ---#
#-------------------------------------------------------------------#
#
# Purpose: Load the sample metadata key and produce a single, cleaned
#          metadata table used consistently by every downstream script
#          (13_phyloseq_objects.R onward).
#
# Input:   data/metadata_key.csv (expected columns: a sample ID
#          column, Distance, Replicate, Transect, and Site/site)
# Output:  data/processed/metadata_clean.rds

library(dplyr)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)


SAMPLE_ID_COLUMN <- "SampleID"

metadata_key <- read.csv("data/metadata_key.csv", stringsAsFactors = FALSE)

if (!SAMPLE_ID_COLUMN %in% colnames(metadata_key) && "SAMPLE" %in% colnames(metadata_key)) {
  metadata_key <- metadata_key %>% rename(SampleID = SAMPLE)
  SAMPLE_ID_COLUMN <- "SampleID"
}

if (!SAMPLE_ID_COLUMN %in% colnames(metadata_key)) {
  stop(
    "Could not find a sample ID column ('", SAMPLE_ID_COLUMN, "' or 'SAMPLE') ",
    "in data/metadata_key.csv. Update SAMPLE_ID_COLUMN at the top of this ",
    "script to match your file."
  )
}

if (!"site" %in% colnames(metadata_key) && "Site" %in% colnames(metadata_key)) {
  metadata_key$site <- metadata_key$Site
}

# Standardize Distance factor levels used throughout the pipeline
if ("Distance" %in% colnames(metadata_key)) {
  metadata_key$Distance <- factor(metadata_key$Distance, levels = c("0", "50", "100"))
}

# Basic column sanity check -- edit this list if your key has more/fewer fields
expected_cols <- c("SampleID", "Distance", "Replicate", "Transect", "site")
missing_cols <- setdiff(expected_cols, colnames(metadata_key))
if (length(missing_cols) > 0) {
  warning(
    "metadata_key.csv is missing expected column(s): ",
    paste(missing_cols, collapse = ", "),
    ". Downstream scripts that rely on these columns may fail."
  )
}

metadata_clean <- metadata_key

message("metadata_clean columns: ", paste(colnames(metadata_clean), collapse = ", "))

saveRDS(metadata_clean, "data/processed/metadata_clean.rds")
