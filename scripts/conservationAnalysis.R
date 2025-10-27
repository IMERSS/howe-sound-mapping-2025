# =========================================
# Visualize vascular plant diversity in Átl’ka7tsem by BEC unit
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
# Source hefty data from Google Drive
# =========================================

# downloadGdrive("", )

# =========================================
# Load vascular plant occurrence summary
# =========================================
BEC.x.plants <- read.csv("tabular_data/BEC_x_vascular_plants_summary_2024.csv")

# Limit analyses to native species
BEC.x.plants.native <- BEC.x.plants %>% filter(establishmentMeans == 'native')

# =========================================
# Create site-species presence/absence matrix for NMDS
# =========================================
site_species_matrix <- BEC.x.plants.native %>%
  dplyr::select(MAP_LABEL, scientificName) %>%
  distinct() %>%
  mutate(presence = 1) %>%
  tidyr::pivot_wider(
    names_from = scientificName,
    values_from = presence,
    values_fill = 0
  ) %>%
  column_to_rownames(var = "MAP_LABEL")

set.seed(42)  # For reproducibility of NMDS

# Run NMDS using Bray-Curtis distance
nmds_result <- metaMDS(site_species_matrix, distance = "bray", k = 2, trymax = 100)

# Check NMDS stress (goodness-of-fit)
nmds_result$stress

# Extract site scores (coordinates) for plotting
nmds_sites <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_sites$MAP_LABEL <- rownames(nmds_sites)

# Extract species scores
nmds_species <- as.data.frame(scores(nmds_result, display = "species"))
nmds_species$species <- rownames(nmds_species)

# Identify common species (occur in ≥3 zones)
species_freq <- colSums(site_species_matrix)
common_species <- names(species_freq[species_freq >= 3])
nmds_species_sub <- nmds_species %>% filter(species %in% common_species)

# =========================================
# Plot NMDS with selected species overlay
# =========================================
ggplot() +
  geom_point(data = nmds_sites, aes(x = NMDS1, y = NMDS2), size = 3, color = "blue") +
  geom_text(data = nmds_sites, aes(x = NMDS1, y = NMDS2, label = MAP_LABEL), vjust = -0.5) +
  geom_point(data = nmds_species_sub, aes(x = NMDS1, y = NMDS2), color = "red", alpha = 0.5, size = 1) +
  geom_text(data = nmds_species_sub, aes(x = NMDS1, y = NMDS2, label = species), color = "red", alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(title = "NMDS of Vegetation by BEC Zone with Species Overlay")

# =========================================
# Compute rarity-weighted richness (RWR)
# =========================================
# Step 1: compute rarity of each species = 1 / (# of zones it occurs in)
species_rarity <- BEC.x.plants.native %>%
  distinct(scientificName, MAP_LABEL) %>%
  group_by(scientificName) %>%
  summarize(
    n_BECs = n_distinct(MAP_LABEL),
    rarity_weight = 1 / n_BECs
  )

# Step 2: compute RWR per zone by summing rarity weights of species present
RWR_BEC <- BEC.x.plants.native %>%
  distinct(scientificName, MAP_LABEL) %>%
  left_join(species_rarity, by = "scientificName") %>%
  group_by(MAP_LABEL) %>%
  summarize(
    n_species = n(),
    RWR = sum(rarity_weight, na.rm = TRUE)
  )

# =========================================
# Identify zone endemics
# =========================================
zone_endemics <- BEC.x.plants.native %>%
  distinct(MAP_LABEL, scientificName) %>%
  group_by(scientificName) %>%
  summarize(n_zones = n(), .groups = "drop") %>%
  right_join(BEC.x.plants.native %>% distinct(MAP_LABEL, scientificName),
             by = "scientificName") %>%
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
beta_jaccard <- vegan::vegdist(site_species_matrix, method = "jaccard")

# Convert to matrix
beta_matrix <- as.matrix(beta_jaccard)

# Compute mean dissimilarity per zone
zone_beta_index <- data.frame(
  MAP_LABEL = rownames(beta_matrix),
  mean_dissimilarity = rowMeans(beta_matrix)
)

# =========================================
# Merge all metrics and normalize composite index
# =========================================
zone_metrics <- RWR_BEC %>%
  left_join(zone_endemics %>% dplyr::select(MAP_LABEL, n_endemics, endemic_prop), by = "MAP_LABEL") %>%
  left_join(zone_beta_index, by = "MAP_LABEL") %>%
  # Normalize RWR to 0–1 by max value
  mutate(
    RWR_norm = RWR / max(RWR),
    # Compute normalized composite index: (normalized RWR) × (mean dissimilarity)
    composite_index = RWR_norm * mean_dissimilarity
  ) %>%
  # Reorder columns for clarity
  dplyr::select(MAP_LABEL, n_species, RWR, RWR_norm, n_endemics, endemic_prop, mean_dissimilarity, composite_index)

# Inspect the merged dataframe
print(zone_metrics)

# Create annotation vector
zone_annotations <- c(
  "Fairly high due to moderate RWR despite fewer unique species",
  "Highest combined biodiversity value: many rare/native species + moderate uniqueness relative to other zones",
  "Low-to-moderate; few rare species and moderate uniqueness",
  "Low; few rare species and high overlap with other zones",
  "Lowest; very low RWR, low number of endemics, though beta uniqueness is highest",
  "Moderate value, primarily due to RWR rather than compositional uniqueness",
  "High value, driven by moderate RWR and relatively high beta diversity",
  "Lower value, few endemics and lower RWR relative to the max",
  "Moderate biodiversity value, RWR and uniqueness contributing roughly equally"
)

# Add to dataframe
zone_metrics$annotation <- zone_annotations

# Inspect
zone_metrics


# Save for downstream analyses
# write.csv(zone_metrics, "outputs/AHSBR_BEC_zone_metrics_normalized.csv", row.names = FALSE)


# Call spatial data

# Layer 1: BEC zones
BEC <- mx_read("spatial_data/vectors/BEC")

# Layer 2: VRI
# VRI <- mx_read("spatial_data/vectors/VRI")

# VRI temporarily rendered as tabular data for analysis until we can pull SHP from Google

VRI <- read.csv('tabular_data/VRI.csv')





# Layer 2: hillshade raster
hillshade <- raster("spatial_data/rasters/Hillshade_80m.tif")

# Layer 3: coastline
coastline <- mx_read("spatial_data/vectors/Islands_and_Mainland")

# Layer 4: watershed boundary
watershed.boundary <- mx_read("spatial_data/vectors/Howe_Sound")

# Sumamrize taxa by BEC Zone

BEC.summary <- BEC.x.plants %>%
  filter(!is.na(MAP_LABEL), !is.na(scientificName)) %>%
  group_by(MAP_LABEL) %>%
  summarise(
    n_species = n_distinct(scientificName),
    species_list = paste(sort(unique(scientificName)), collapse = ", ")
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
  summarize(taxa = paste(scientificName, collapse = ", "), .groups = "drop")

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
