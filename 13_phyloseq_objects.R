#-------------------------------------------------------------------#
#--- Stage 13: Assemble analysis-ready phyloseq objects           ---#
#-------------------------------------------------------------------#
#
# Purpose: Load the per-marker phyloseq objects (16S/bacterial, COI,
#          rbcL, 18S-V7, 18S-V9 -- output of the upstream DADA2 /
#          decontamination / read-depth filtering pipeline described
#          in the manuscript Methods, translate their lab/sequencing
#          sample IDs into field IDs, and attach the cleaned sample
#          metadata from 12_metadata_processing.R. Produces the single
#          list of analysis-ready phyloseq objects used by every
#          downstream analysis script (14, 15, 16). 

# Depends on: 12_metadata_processing.R (data/processed/metadata_clean.rds)
#
# Input:   data/bacterial.rds, data/COI.rds, data/rbcL.rds,
#          data/V7_18S.rds, data/V9_18S.rds
#          data/processed/metadata_clean.rds
#          data/sequencing submission forms
#       
# Output:  data/processed/physeq_list.rds


library(phyloseq)
library(dplyr)
library(readxl)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

metadata_clean <- readRDS("data/processed/metadata_clean.rds")

raw_files <- c(
  Bacterial = "data/bacterial.rds",
  COI       = "data/COI.rds",
  rbcL      = "data/rbcL.rds",
  V7_18S    = "data/V7_18S.rds",
  V9_18S    = "data/V9_18S.rds"
)

#=====================================================================
# Sample ID crosswalk (lab/sequencing ID -> field ID)
#=====================================================================

CROSSWALK_SPEC <- list(
  Bacterial = list(file = "data/sample submission.xlsx",  plates = 1:2, suffix = "_16S"),
  COI       = list(file = "data/sample submission.xlsx",  plates = 3:4, suffix = "_CO1"),
  V9_18S    = list(file = "data/sample submission.xlsx", plates = 1:2, suffix = "_v9"),
  V7_18S    = list(file = "data/sample submission.xlsx", plates = 3:4, suffix = "_v7"),
  rbcL      = list(file = "data/sample submission.xlsx", plates = 5:6, suffix = "_rbcL")
)

build_crosswalk <- function(file, plates, suffix) {
  plate_sheets <- grep("^Plate", readxl::excel_sheets(file), value = TRUE)
  if (max(plates) > length(plate_sheets)) {
    stop(
      "build_crosswalk(): ", file, " only has ", length(plate_sheets),
      " 'Plate...' sheet(s), but plate(s) ", paste(plates, collapse = ", "),
      " were requested. The submission form layout may have changed --",
      " update CROSSWALK_SPEC in this script."
    )
  }

  rows <- lapply(plate_sheets[plates], function(sheet) {
    df <- readxl::read_excel(file, sheet = sheet, col_names = FALSE, skip = 3)
    data.frame(lab_id = df[[2]], sample_name = df[[4]], stringsAsFactors = FALSE)
  })

  crosswalk <- do.call(rbind, rows)
  crosswalk <- crosswalk[!is.na(crosswalk$lab_id), ]
  crosswalk$field_id <- sub(paste0(suffix, "$"), "", crosswalk$sample_name)
  crosswalk
}

crosswalks <- lapply(CROSSWALK_SPEC, function(spec) {
  build_crosswalk(spec$file, spec$plates, spec$suffix)
})


KNOWN_EXCLUSIONS <- list(
  V7_18S = c("CAW-25-03-P4-51")  # duplicate submission of field ID I31003
)

# Translate a phyloseq object's lab/sequencing sample names to field
# IDs using its marker's crosswalk, preserving the original lab ID in
# sample_data (column `lab_id`) for traceability. 

translate_lab_ids <- function(ps, crosswalk) {
  lab_ids <- sample_names(ps)
  lookup <- setNames(crosswalk$field_id, crosswalk$lab_id)
  field_ids <- unname(lookup[lab_ids])
  field_ids[is.na(field_ids)] <- lab_ids[is.na(field_ids)]

  dup <- duplicated(field_ids) & field_ids %in% crosswalk$field_id
  if (any(dup)) {
    message(
      "translate_lab_ids(): dropping duplicate submission(s) -- lab ID(s) [",
      paste(lab_ids[dup], collapse = ", "), "] resolve to field ID(s) [",
      paste(field_ids[dup], collapse = ", "), "] already seen earlier in the ",
      "submission form. Keeping only the first occurrence (matches the ",
      "confirmed precedent for V7_18S -- see FLAG at top of script)."
    )
    ps <- prune_samples(lab_ids[!dup], ps)
    field_ids <- field_ids[!dup]
    lab_ids <- lab_ids[!dup]
  }

  sample_data(ps)$lab_id <- lab_ids
  sample_names(ps) <- field_ids
  ps
}

#=====================================================================
# Attach metadata
#=====================================================================

attach_metadata <- function(ps, metadata) {
  meta <- data.frame(sample_data(ps)) %>%
    tibble::rownames_to_column("SampleID")

  overlapping <- intersect(setdiff(colnames(metadata), "SampleID"), colnames(meta))
  if (length(overlapping) > 0) {
    message(
      "attach_metadata(): dropping pre-existing column(s) [",
      paste(overlapping, collapse = ", "),
      "] from the raw phyloseq sample_data before merging metadata_clean."
    )
    meta <- meta %>% select(-all_of(overlapping))
  }

  meta <- meta %>%
    left_join(metadata, by = "SampleID") %>%
    tibble::column_to_rownames("SampleID")

  sample_data(ps) <- meta
  ps
}


set_distance_levels <- function(ps) {
  sample_data(ps)$Distance <- factor(sample_data(ps)$Distance, levels = c("0", "50", "100"))
  ps
}

physeq_list <- lapply(names(raw_files), function(dataset_name) {
  ps <- readRDS(raw_files[[dataset_name]])

  exclude <- KNOWN_EXCLUSIONS[[dataset_name]]
  if (!is.null(exclude)) {
    present <- intersect(exclude, sample_names(ps))
    if (length(present) > 0) {
      message(
        dataset_name, ": applying known exclusion (ported from the ",
        "original author's script) -- dropping lab ID(s): ",
        paste(present, collapse = ", ")
      )
      ps <- prune_samples(setdiff(sample_names(ps), present), ps)
    }
  }

  ps <- translate_lab_ids(ps, crosswalks[[dataset_name]])
  ps <- attach_metadata(ps, metadata_clean)
  ps <- set_distance_levels(ps)
  ps
})
names(physeq_list) <- names(raw_files)


UNMATCHED_FRACTION_LIMIT <- 0.5  # stop if more than this fraction fails to match

for (dataset_name in names(physeq_list)) {
  ps <- physeq_list[[dataset_name]]
  meta <- data.frame(sample_data(ps))
  no_match <- rownames(meta)[is.na(meta$site)]

  if (length(no_match) / nsamples(ps) > UNMATCHED_FRACTION_LIMIT) {
    stop(
      dataset_name, ": ", length(no_match), " of ", nsamples(ps),
      " samples (", round(100 * length(no_match) / nsamples(ps)), "%) have no ",
      "match in metadata_key.csv after ID translation -- this looks like a ",
      "deeper ID problem, not a few stray controls, so refusing to continue.\n",
      "  Example unmatched (post-translation) sample IDs: ",
      paste(head(no_match, 5), collapse = ", "), "\n",
      "  Example metadata_key.csv SampleIDs:                ",
      paste(head(metadata_clean$SampleID, 5), collapse = ", "), "\n",
      "Check CROSSWALK_SPEC for '", dataset_name, "' against the actual ",
      "submission form layout."
    )
  }

  if (length(no_match) > 0) {
    message(
      dataset_name, ": dropping ", length(no_match),
      " sample(s) with no metadata_key.csv match (translated IDs): ",
      paste(no_match, collapse = ", ")
    )
    physeq_list[[dataset_name]] <- prune_samples(setdiff(sample_names(ps), no_match), ps)
  }
}


for (dataset_name in names(physeq_list)) {
  cols <- colnames(data.frame(sample_data(physeq_list[[dataset_name]])))
  if (!"site" %in% cols) {
    warning(
      "`site` column missing from ", dataset_name, " after attach_metadata(). ",
      "Columns present: ", paste(cols, collapse = ", ")
    )
  }
}

saveRDS(physeq_list, "data/processed/physeq_list.rds")
