# Visualize vascular plant diversity in Átl’ka7tsem by BEC unit

# Set relative paths (https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location-in-rstudio)
if (!isTRUE(getOption('knitr.in.progress'))) {
  setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/.."))
}

# Load libraries

library(dplyr)
library(leaflet)
library(raster)
library(reshape2)
library(scales)
library(sf)
library(jsonlite)
library(viridis)
library(htmlwidgets)

# Source dependencies

source("scripts/utils.R")

# Call tabular data

BEC.x.plants <- read.csv("tabular_data/BEC_x_vascular_plants_summary_2024.csv")

# Call spatial data

# Layer 1: BEC zones
BEC <- mx_read("spatial_data/vectors/BEC")

# Layer 2: hillshade raster
hillshade <- raster("spatial_data/rasters/Hillshade_80m.tif")

# Layer 3: coastline
coastline <- mx_read("spatial_data/vectors/Islands_and_Mainland")

# Layer 4: watershed boundary
watershed.boundary <- mx_read("spatial_data/vectors/Howe_Sound")

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
