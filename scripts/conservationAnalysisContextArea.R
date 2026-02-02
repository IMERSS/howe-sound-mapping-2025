# =========================================
# Compute ecological indices for Howe Sound vegetation communities
# (with effort-based rarefaction)
# =========================================

# Set working directory to project root (parent of current script)
if (!isTRUE(getOption('knitr.in.progress'))) {
  setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/.."))
}

# =========================================
# Load libraries
# =========================================

library(dplyr)
library(htmlwidgets)
library(jsonlite)
library(leaflet)
library(tidyverse)
library(raster)
library(reshape2)
library(scales)
library(sf)
library(vegan)
library(viridis)

# Source any custom functions
source("scripts/utils.R")

# =========================================
# Load vascular plant occurrence summary
# =========================================

BEC.x.plants <- read.csv("tabular_data/vascular-plants_FullStudyRegion.csv")

# Limit analyses to native species # RESTORE THIS FOR ANALYSIS LATER
# BEC.x.plants.native <- BEC.x.plants %>%
#  filter(establishmentMeans == 'native')

# =========================================
# RAREFACTION TO CONTROL FOR SAMPLING EFFORT
# =========================================

## NOTE DUE TO PROBLEMS W DATASET WE CANNOT PROCEED WITH RAREFACTION FOR NOW (TOO FEW SPECIES IN SOME ZONES)

# Quantify sampling effort per BEC zone
effort_by_zone <- BEC.x.plants %>%
  count(MAP_LABEL, name = "n_obs")

# Inspect effort distribution (recommended to run interactively)
print(summary(effort_by_zone$n_obs))

# Choose rarefaction depth (minimum effort across zones)
min_effort <- min(effort_by_zone$n_obs)

set.seed(42)  # Reproducibility

# Rarefy: equal number of observations per zone
BEC.x.plants.rarefied <- BEC.x.plants %>%
  group_by(MAP_LABEL) %>%
  slice_sample(n = min_effort) %>%
  ungroup()

# =========================================
# Create site-species presence/absence matrix for NMDS
# =========================================

site_species_matrix <- BEC.x.plants %>%
  dplyr::select(MAP_LABEL, species) %>% # NOTE: THIS NEW DATASET IS NOT FORMATTED CORRECTLY SO I AM RESORTING TO 'species' FOR NOW
  distinct() %>%
  mutate(presence = 1) %>%
  tidyr::pivot_wider(
    names_from = species,
    values_from = presence,
    values_fill = 0
  ) %>%
  column_to_rownames(var = "MAP_LABEL")

set.seed(42)  # For reproducibility of NMDS

# Run NMDS using Bray-Curtis distance
nmds_result <- metaMDS(
  site_species_matrix,
  distance = "bray",
  k = 2,
  trymax = 100
)

# Extract site scores (coordinates)
nmds_sites <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_sites$MAP_LABEL <- rownames(nmds_sites)

# Extract species scores
nmds_species <- as.data.frame(scores(nmds_result, display = "species"))
nmds_species$species <- rownames(nmds_species)

# Identify common vs. less common species
species_freq <- colSums(site_species_matrix)
common_species <- names(species_freq[species_freq >= 3])
rare_species   <- names(species_freq[species_freq < 3])

nmds_species$group <- ifelse(
  nmds_species$species %in% common_species,
  "Common (≥3 zones)",
  "Less common (<3 zones)"
)

# =========================================
# Compute rarity-weighted richness (RWR)
# =========================================

# Step 1: species rarity = inverse of number of zones occupied
species_rarity <- BEC.x.plants %>%
  distinct(species, MAP_LABEL) %>%
  group_by(species) %>%
  summarize(
    n_BECs = n_distinct(MAP_LABEL),
    rarity_weight = 1 / n_BECs,
    .groups = "drop"
  )

# Step 2: sum rarity weights per zone
RWR_BEC <- BEC.x.plants %>%
  distinct(species, MAP_LABEL) %>%
  left_join(species_rarity, by = "species") %>%
  group_by(MAP_LABEL) %>%
  summarize(
    n_species = n(),
    RWR = sum(rarity_weight, na.rm = TRUE),
    .groups = "drop"
  )

# =========================================
# Identify zone endemics
# =========================================

zone_endemics <- BEC.x.plants %>%
  distinct(MAP_LABEL, species) %>%
  group_by(species) %>%
  summarize(
    n_zones = n(),
    .groups = "drop"
  ) %>%
  right_join(
    BEC.x.plants %>% distinct(MAP_LABEL, species),
    by = "species"
  ) %>%
  group_by(MAP_LABEL) %>%
  summarize(
    n_species = n(),
    n_endemics = sum(n_zones == 1),
    endemic_prop = n_endemics / n_species,
    .groups = "drop"
  )

# =========================================
# Compute Jaccard dissimilarity between zones
# =========================================

beta_jaccard <- vegan::vegdist(
  site_species_matrix,
  method = "jaccard"
)

beta_matrix <- as.matrix(beta_jaccard)

zone_beta_index <- data.frame(
  MAP_LABEL = rownames(beta_matrix),
  mean_dissimilarity = rowMeans(beta_matrix)
)

# =========================================
# Merge all metrics and compute composite index
# =========================================

zone_metrics <- RWR_BEC %>%
  left_join(
    zone_endemics %>% dplyr::select(MAP_LABEL, n_endemics, endemic_prop),
    by = "MAP_LABEL"
  ) %>%
  left_join(zone_beta_index, by = "MAP_LABEL") %>%
  mutate(
    RWR_norm = RWR / max(RWR),
    composite_index = RWR_norm * mean_dissimilarity
  ) %>%
  dplyr::select(
    MAP_LABEL,
    n_species,
    RWR,
    RWR_norm,
    n_endemics,
    endemic_prop,
    mean_dissimilarity,
    composite_index
  )

print(zone_metrics, n = Inf)

zone_annotations <- c(
  # 1 BAFAun
  "Low overall biodiversity value; species-poor zone with low rarity-weighted richness and a single species confined to this zone within the study area, though composition is moderately distinct from other zones.",
  
  # 2 BGxh3
  "Very low biodiversity value; few species, no zone-restricted taxa, and low rarity-weighted richness, despite moderate compositional distinctness.",
  
  # 3 CDFmm
  "Very high biodiversity value; exceptionally species-rich with high rarity-weighted richness and multiple species confined to this zone within the study area, representing a floristically distinctive and regionally important community.",
  
  # 4 CMAunp
  "Moderate-to-high biodiversity value; relatively rich community with several zone-restricted species and moderate rarity-weighted richness, though species composition overlaps substantially with other zones.",
  
  # 5 CWHdm
  "Highest combined biodiversity value; the most rarity-rich zone with the greatest number of species confined to this zone within the study area and strong compositional distinctness, indicating exceptional conservation importance.",
  
  # 6 CWHds1
  "Moderate biodiversity value; intermediate rarity-weighted richness with few zone-restricted species and slightly reduced compositional distinctness relative to other coastal forest zones.",
  
  # 7 CWHms1
  "Moderate biodiversity value; moderately high rarity-weighted richness driven primarily by overall richness, but very few zone-restricted species and only moderate compositional uniqueness.",
  
  # 8 CWHvm1
  "Moderate biodiversity value; similar to other CWH variants with moderate rarity-weighted richness, few zone-restricted species, and average compositional differentiation.",
  
  # 9 CWHvm2
  "Moderate-to-high biodiversity value; elevated rarity-weighted richness with several species confined to this zone within the study area and moderately strong compositional distinctness.",
  
  # 10 CWHxm1
  "Very high biodiversity value; high rarity-weighted richness and numerous species confined to this zone within the study area combined with strong compositional uniqueness, ranking among the most important zones regionally.",
  
  # 11 ESSFdh2
  "Low biodiversity value; limited species richness, no zone-restricted taxa, and low rarity-weighted richness, despite moderate differentiation from other zones.",
  
  # 12 ESSFdv1
  "Low biodiversity value; species-poor zone with minimal rarity-weighted richness and no species confined to this zone within the study area.",
  
  # 13 ESSFdv2
  "Low-to-moderate biodiversity value; modest species richness but low rarity-weighted richness and no zone-restricted species.",
  
  # 14 ESSFdvp
  "Low-to-moderate biodiversity value; limited rarity and no zone-restricted species, with only moderate compositional distinctness.",
  
  # 15 ESSFdvw
  "Low biodiversity value; relatively species-poor with no zone-restricted taxa and reduced compositional differentiation.",
  
  # 16 ESSFmw2
  "Moderate biodiversity value; intermediate rarity-weighted richness with a small number of species confined to this zone within the study area, though compositional uniqueness is modest.",
  
  # 17 ESSFmwp
  "Low-to-moderate biodiversity value; limited richness and rarity, with a single species confined to this zone within the study area contributing marginally to overall biodiversity value.",
  
  # 18 ESSFmww
  "Low biodiversity value; species-poor zone lacking zone-restricted species and exhibiting low rarity-weighted richness.",
  
  # 19 ESSFxc3
  "Low-to-moderate biodiversity value; few species overall but containing a single species confined to this zone within the study area, yielding modest rarity-weighted richness.",
  
  # 20 ESSFxcp
  "Very low biodiversity value; minimal species richness, no zone-restricted taxa, and very low rarity-weighted richness.",
  
  # 21 ESSFxcw
  "Very low biodiversity value; extremely species-poor with no species confined to this zone within the study area and negligible rarity-weighted richness.",
  
  # 22 ESSFxv2
  "Low biodiversity value; limited richness and no zone-restricted species, despite moderate compositional differentiation.",
  
  # 23 ESSFxvp
  "Very low biodiversity value; few species and no zone-restricted taxa, contributing little to regional biodiversity patterns.",
  
  # 24 ESSFxvw
  "Extremely low biodiversity value; exceptionally species-poor zone with negligible rarity-weighted richness and no species confined to this zone within the study area.",
  
  # 25 IDFdc
  "Low-to-moderate biodiversity value; moderate species richness but low rarity-weighted richness and no zone-restricted species.",
  
  # 26 IDFdk1
  "Low biodiversity value; species-poor zone with minimal rarity-weighted richness and no zone-restricted taxa.",
  
  # 27 IDFdk2
  "Extremely low biodiversity value; represented by a single species with negligible contribution to regional rarity or compositional distinctness.",
  
  # 28 IDFdk3
  "Very low biodiversity value; few species and no zone-restricted taxa, yielding minimal rarity-weighted richness.",
  
  # 29 IDFww
  "Moderate biodiversity value; relatively species-rich interior forest zone, though rarity-weighted richness is low and no species are confined to this zone within the study area.",
  
  # 30 IDFww1
  "Low-to-moderate biodiversity value; limited rarity and no zone-restricted species despite moderate richness.",
  
  # 31 IDFxc
  "Moderate biodiversity value; moderate rarity-weighted richness with a single species confined to this zone within the study area, but limited compositional uniqueness.",
  
  # 32 IDFxh2
  "Low biodiversity value; species-poor zone with low rarity-weighted richness and no zone-restricted taxa.",
  
  # 33 IDFxw
  "Low biodiversity value; limited species richness and no zone-restricted taxa, contributing little to regional compositional distinctness.",
  
  # 34 IMAun
  "Low-to-moderate biodiversity value; modest richness but low rarity-weighted richness and no species confined to this zone within the study area.",
  
  # 35 IMAunp
  "Low-to-moderate biodiversity value; moderate richness with a single species confined to this zone within the study area, though overall rarity and compositional distinctness are limited.",
  
  # 36 MHmm1
  "Moderate biodiversity value; relatively rich montane zone but dominated by widespread species, resulting in few species restricted to this zone within the study area.",
  
  # 37 MHmm2
  "Moderate biodiversity value; moderate richness with a small number of species confined to this zone within the study area, but limited compositional differentiation.",
  
  # 38 MSdc1
  "Low biodiversity value; modest richness but low rarity-weighted richness and no zone-restricted species.",
  
  # 39 MSdc3
  "Low biodiversity value; similar to MSdc1 with limited rarity and no zone-restricted taxa.",
  
  # 40 MSxk3
  "Very low biodiversity value; species-poor zone with negligible rarity-weighted richness and no species confined to this zone within the study area.",
  
  # 41 PPxh2
  "Moderate biodiversity value; moderate rarity-weighted richness and a small number of species confined to this zone within the study area, contributing modestly to regional uniqueness."
)


zone_metrics$annotation <- zone_annotations

# Save for downstream analyses
write.csv(zone_metrics, "outputs/Context_Area_BEC_zone_metrics.csv", row.names = FALSE)


# Downscale biodiversity indices to map with VRI mapping data

# Load VRI and BEC veg indices

VRI <- read.csv('tabular_data/AHSBR_BEC-site-series.csv')

BEC <- read.csv("outputs/AHSBR_BEC_zone_metrics_normalized_sampling_effort.csv")

VRI <- VRI %>%
  left_join(
    BEC,
    by = c("BEC_LABEL" = "MAP_LABEL")
  )

# Apply penalty for structural stage
VRI <- VRI %>%
  mutate(
    stage_penalty = case_when(
      TEM_StrClass == "Sparse/Bryoid (Rock, Ice, Moss)" ~ 0,
      is.na(Age_Class) ~ NA_real_,
      TRUE ~ (9 - Age_Class) / 9
    )
  )

# Assign weightings for CDC Ranked ecological communities

# Helper to convert BC list to numeric weight
bc_weight <- function(x) {
  case_when(
    x == "Red"  ~ 2.0,
    x == "Blue" ~ 1.5,
    TRUE        ~ 1.0   # Yellow, NA, "", anything else
  )
}

VRI <- VRI %>%
  mutate(
    conservation_multiplier =
      S1_DECn * bc_weight(S1_BCList) +
      S2_DECn * bc_weight(S2_BCList) +
      S3_DECn * bc_weight(S3_BCList)
  )

VRI <- VRI %>%
  mutate(
    stage_modifier = 1 - 0.5 * stage_penalty,  
    final_biodiversity_index =
      composite_index *
      stage_modifier *
      conservation_multiplier
  )

# Normalize
VRI <- VRI %>%
  mutate(
    final_biodiversity_index_norm =
      (final_biodiversity_index - min(final_biodiversity_index, na.rm = TRUE)) /
      (max(final_biodiversity_index, na.rm = TRUE) - min(final_biodiversity_index, na.rm = TRUE))
  )


# Assess influence of modifiers on final biodiversity index
ggplot(VRI, aes(x = composite_index,
               y = final_biodiversity_index,
               color = stage_penalty)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_viridis_c(option = "C") +
  coord_equal() +
  labs(
    x = "Composite biodiversity index (unmodified)",
    y = "Final biodiversity index (with modifiers)",
    color = "Stage penalty",
    title = "Effect of stage and conservation modifiers on biodiversity index"
  ) +
  theme_minimal()

VRI %>%
  mutate(conservation_effect =
           final_biodiversity_index /
           (composite_index * stage_modifier)) %>%
  ggplot(aes(x = conservation_effect)) +
  geom_histogram(bins = 30, alpha = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  labs(
    x = "Conservation multiplier effect (conditional on stage)",
    y = "Count",
    title = "Distribution of conservation effects after accounting for stage"
  ) +
  theme_minimal()


# Write output
write.csv(VRI, "outputs/AHSBR_VRI_vascular_plant_biodiversity_indices.csv", row.names = FALSE)




# Call spatial data

# Layer 1: BEC zones
BEC <- mx_read("spatial_data/vectors/BEC")

# VRI temporarily rendered as tabular data for analysis until we can pull SHP from Google

# Layer 2: VRI
# VRI <- mx_read("spatial_data/vectors/VRI")


# Manipulate VRI as CSV / table for now until we can deal with the hefty dataset

# VRI <- st_drop_geometry(VRI)
# write.csv(VRI, "tabular_data/VRI.csv")













# Layer 2: hillshade raster
hillshade <- raster("spatial_data/rasters/Hillshade_80m.tif")

# Layer 3: coastline
coastline <- mx_read("spatial_data/vectors/Islands_and_Mainland")

# Layer 4: watershed boundary
watershed.boundary <- mx_read("spatial_data/vectors/Howe_Sound")

# Sumamrize taxa by BEC Zone

BEC.summary <- BEC.x.plants %>%
  filter(!is.na(MAP_LABEL), !is.na(species)) %>%
  group_by(MAP_LABEL) %>%
  summarise(
    n_species = n_distinct(species),
    species_list = paste(sort(unique(species)), collapse = ", ")
  ) %>%
  arrange(desc(n_species))





# Create map labels

becLabels <- list(CMAunp = "Coastal Mountain Heather Alpine Zone",
                  CWHdm  = "Dry Maritime Coastal Western Hemlock Zone",
                  CWHds1 = "Southern Dry Submaritime Coastal Western Hemlock Zone",
                  CWHms1 = "Southern Moist Submaritime Coastal Western Hemlock Zone",
                  CWHvm1 = "Submontane Very Wet Maritime Coastal Western Hemlock Zone",
                  CWHvm2 = "Montane Very Wet Maritime Coastal Western Hemlock Zone",
                  CWHxm1 = "Eastern Very Dry Maritime Coastal Western Hemlock Zone",
                  MHmm1 = "Windward Moist Maritime Mountain Hemlock Zone",
                  MHmm2 = "Leeward Moist Maritime Mountain Hemlock Zone")

# Create color palette for BEC Zones

# Following rough elevational gradient:  
# CDFmm, CWHxm1, CWHdm, CWHvm1, CWHvm2, CWHds1, CWHms1, MHmm1, MHmm2, ESSFmw2, CMAunp

palette = data.frame(
  cat = c("CWHxm1","CWHdm","CWHvm2","CWHvm1","CWHds1","CWHms1","MHmm1","MHmm2","CMAunp"), 
  col = c("#440154FF", "#472D7BFF","#3B528BFF","#2C728EFF","#21908CFF","#27AD81FF","#5DC863FF","#AADC32FF", "#FDE725FF")
)

# Reduce plants.x.BEC dataframe for visualization

BEC.x.plants$MAP_LABEL <- as.factor(BEC.x.plants$MAP_LABEL)

taxa <- BEC.x.plants %>%
  group_by(MAP_LABEL) %>%
  summarize(taxa = paste(species, collapse = ", "), .groups = "drop")

# Generate visualization data and aesthetics

paletteHash <- split(x = palette$col, f = palette$cat)
taxaHash <- split(x = taxa$taxa, f = taxa$MAP_LABEL)

vascularData <- list(palette = paletteHash, taxa = taxaHash, mapTitle = "Map 1. Ecological communities", regionLabels = becLabels)

BEC <- base::merge(BEC, palette, by.x ="MAP_LABEL", by.y="cat")

# Write plant x BEC Zone summary to JSON file for viz

write_json(vascularData, "viz_data/Vascular-plotData.json")

# Plot map

speciesMap <- leaflet(options=list(mx_mapId="Vascular")) %>%
  setView(-123.2194, 49.66076, zoom = 8.5) %>%
  addTiles(options = providerTileOptions(opacity = 0.5)) %>%
  addRasterImage(hillshade, opacity = 0.8) %>%
  addPolygons(data = coastline, color = "black", weight = 1.5, fill = FALSE) %>%
  addPolygons(data = BEC, fillColor = BEC$col, fillOpacity = 0.6, weight = 2, stroke = FALSE, options = labelToOption(BEC$MAP_LABEL)) %>%
  addPolygons(data = watershed.boundary, color = "black", weight = 4, fill = FALSE)

#Note that this statement is only effective in standalone R
print(speciesMap)
