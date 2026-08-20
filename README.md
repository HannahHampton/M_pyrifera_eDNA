[README.md](https://github.com/user-attachments/files/31284946/README.md)
 Leveraging environmental DNA to unveil the role of giant kelp habitats as a blue carbon sink

Code accompanying the manuscript:

> Hampton HG, Crossett D, Giles EC, Pochon X, Laroche O, Scriver M, Kaur G, Geraldi NR, Jedrecka T, Zhang XB, Zaiko A. *Leveraging environmental DNA to unveil the role of giant kelp habitats as a blue carbon sink.* (manuscript in preparation)

Affiliations: Cawthron Institute (Nelson, NZ); University of Auckland Institute of Marine Science; Sequench Ltd (Nelson, NZ); Naturemetrics (Guildford, UK); Kelp Forest Foundation (Netherlands).

## Abstract

Giant kelp forests formed by *Macrocystis pyrifera* are foundational marine ecosystems that sustain biodiversity and generate organic carbon, with a portion transported and potentially sequestered in sediments beyond the habitats where the kelp grows. Widespread tracking of *M. pyrifera*-derived carbon is a challenge but could be aided by molecular tools. Here, we developed a digital droplet PCR (ddPCR) assay targeting the internal transcribed spacer 2 (ITS2) region, identified using comparative genomics and phylogenetic analyses. Application of the assay to sediment samples from sites with varying *M. pyrifera* presence revealed a gradient of *M. pyrifera* environmental DNA (eDNA) detection. Occupancy modelling indicated that *M. pyrifera* presence was the primary covariate influencing *M. pyrifera* eDNA presence in sediments. The estimated probability of detection was 0.6, comparable to detectability estimates reported for other marine species using eDNA assays. Complementary eDNA metabarcoding showed distinct community structures among locations and identified several taxa strongly associated with *M. pyrifera* habitats. However, metabarcoding failed to identify *M. pyrifera* at genus or species level, likely due to gaps in current reference databases, highlighting the need for more species-specific diagnostic tools such as qPCR or ddPCR. Overall, this study represents the first species-specific assay validated on marine sediments for the reliable detection of *M. pyrifera*, demonstrating its potential for monitoring *M. pyrifera* distributions and contributions to carbon sequestration in marine sediments.

## Repository contents


.
├── scripts/          Numbered R pipeline, run in order 
├── data/             Expected input data structure (see data/README.md) — raw data files are not included
│   └── processed/    Intermediate .rds files passed between pipeline stages — not tracked in git
└── outputs/          Final figures/ and tables/, created by the scripts at runtime — not tracked in git


## Scripts

The pipeline is organized as a sequence of numbered stages, each reading the previous stage's saved output from `data/processed/` and writing its own. Run them in order from the repository root (so relative paths resolve correctly).

| Script | Stage | Description |
|---|---|---|
| `scripts/12_metadata_processing.R` | Metadata | Loads and standardizes `data/metadata_key.csv` (sample ID column, `Distance` factor levels, `site` column) into a single cleaned metadata table used by every later stage. |
| `scripts/13_phyloseq_objects.R` | Phyloseq assembly | Loads the per-marker phyloseq objects (16S/bacterial, COI, rbcL, 18S-V7, 18S-V9), translates their lab/sequencing sample IDs to field IDs using a crosswalk built from the Sequench submission forms, and attaches the cleaned metadata, producing one analysis-ready `physeq_list`. |
| `scripts/14_ordination.R` | Ordination | Runs PCoA (Bray-Curtis) per marker and computes convex hulls by site — analysis only, no plotting. |
| `scripts/15_permanova_permdisp.R` | Stats | Pairwise PERMANOVA (`adonis2`) between sites per marker, plus a PERMDISP (`betadisper`) test. |
| `scripts/16_indicator_species.R` | Indicator species | Indicator species analysis (`indicspecies::multipatt`) comparing *M. pyrifera*-present vs. absent/low sites. |
| `scripts/17_ddPCR_analysis.R` | ddPCR | Processes the ddPCR LOD/LOQ dilution series, field sediment-sample detection results, and the metabarcoding-vs-ddPCR Laminariaceae comparison. |
| `scripts/18_tables.R` | Tables | Exports all result tables (PERMANOVA/PERMDISP, ddPCR summaries, indicator species) to `outputs/tables/`. |
| `scripts/19_figures.R` | Figures | Renders all publication figures (ddPCR LOD/LOQ, field detections, four ordination-figure variants, Laminariaceae comparison, indicator species Venn diagrams, Figure 6 class-composition bar charts) to `outputs/figures/`. |


## Software

Analyses were performed in R. 

## Data availability

Raw metabarcoding sequence data are deposited in the NCBI Sequence Read Archive under BioProject **PRJNA1402725**. Whole-genome resequencing data used for primer/probe design are available under the NCBI accessions and Aotearoa Genomics Data Repository ID listed in the manuscript. Processed data files referenced by these scripts (ddPCR results, cleaned phyloseq objects, sample metadata) are not included in this repository — see `data/README.md` for the expected file structure, and contact the corresponding author for access.

## License

Code in this repository is released under the MIT License — see [LICENSE](LICENSE).

## Contact

Hannah Hampton, Cawthron Institute — hannah.hampton@cawthron.org.nz
