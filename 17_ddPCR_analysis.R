#-------------------------------------------------------------------#
#--- Stage 17: ddPCR data analysis (LOD/LOQ + field detections)   ---#
#-------------------------------------------------------------------#
#
# Purpose: Process the ddPCR-related datasets used in the manuscript --
#          the dilution-series LOD/LOQ characterization , the field sediment-sample
#          detections, and the metabarcoding-vs-ddPCR Laminariaceae comparison 
#          into the summary objects consumed by 18_tables.R and 19_figures.R. 
#
# Input:   data/Probe LOQ LOD R.csv
#          data/Sample analysis.csv
#          data/Comparison_metabarcoding_ddPCR.csv
# Output:  data/processed/lod_loq_summary.rds
#          data/processed/positive_detections.rds
#          data/processed/laminariaceae_comparison.rds
#

library(dplyr)

dir.create("data/processed", recursive = TRUE, showWarnings = FALSE)

# ---- LOD/LOQ dilution series  ----
lod_loq_raw <- read.csv("data/Probe LOQ LOD R.csv")

lod_loq_summary <- lod_loq_raw %>%
  group_by(Inputconc.) %>%
  summarise(
    Mean = mean(copiesper),
    SD = sd(copiesper),
    SE = SD / sqrt(n()) 
  )

desired_order <- c(
  "14ng", "1.4ng", "0.14ng",
  "14pg", "1.4pg", "0.14pg",
  "14fg", "1.4fg", "0.14fg",
  "14ag", "1.4ag", "0.14ag",
  "NTC"
)

lod_loq_summary$Inputconc. <- factor(lod_loq_summary$Inputconc., levels = desired_order)

lod_loq_summary$Mean1 <- lod_loq_summary$Mean
lod_loq_summary$Mean1[lod_loq_summary$Mean1 == 0] <- NA

saveRDS(lod_loq_summary, "data/processed/lod_loq_summary.rds")

# ---- Field ddPCR detections ----
field_data <- read.csv("data/Sample analysis.csv")

field_data$Location <- factor(
  field_data$Location,
  levels = c("Pelorus", "Tory", "OuterQC", "InnerQC", "Pos", "NTC")
)


new_na <- sum(is.na(field_data$Location))
if (new_na > 0) {
  stop(
    "field_data$Location: ", new_na, " row(s) did not match any of the ",
    "expected Location values and became NA. Check the raw values with ",
    "table(read.csv('data/Sample analysis.csv')$Location) and update the ",
    "levels= list above to match."
  )
}

positive_only <- field_data %>%
  filter(Conc.copies.µL. >= 0.2) %>%
  filter(!Location %in% c("Pos", "NTC")) %>%
  mutate(Location = factor(
    Location,
    levels = c("Pelorus", "Tory", "OuterQC", "InnerQC"),
    labels = c("PS", "TC", "Outer QC", "Inner QC")
  ))

saveRDS(positive_only, "data/processed/positive_detections.rds")

# ---- Metabarcoding-vs-ddPCR Laminariaceae comparison  ----
# Compares the number of Laminariaceae OTUs detected by metabarcoding
# against the ddPCR results, by gene and by location.
laminariaceae_comparison <- read.csv("data/Comparison_metabarcoding_ddPCR.csv")

laminariaceae_comparison$Location <- factor(
  laminariaceae_comparison$Location,
  levels = c("Pelorus", "Tory", "InnerQC", "OuterQC")
)

saveRDS(laminariaceae_comparison, "data/processed/laminariaceae_comparison.rds")
