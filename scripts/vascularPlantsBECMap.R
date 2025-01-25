# Visualize vascular plant diversity in Átl’ka7tsem by BEC unit

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

# Call map layers

# Layer 1: BEC zone x vascular plants
BEC <- mx_read("spatial_data/vectors/BEC")

# Layer 2: BEC zones
plants.x.BEC <- mx_read("spatial_data/vectors/BEC_x_plants_2024")

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

plants.x.BEC$MAP_LAB <- as.factor(plants.x.BEC$MAP_LAB)

palette = data.frame(
  cat = c("CWHxm1","CWHdm","CWHvm2","CWHvm1","CWHds1","CWHms1","MHmm1","MHmm2","CMAunp"), 
  col = c("#440154FF", "#472D7BFF","#3B528BFF","#2C728EFF","#21908CFF","#27AD81FF","#5DC863FF","#AADC32FF", "#FDE725FF")
)

# Reduce plants.x.BEC dataframe for visualization
taxa <- plants.x.BEC %>%
  group_by(MAP_LAB) %>%
  summarize(taxa = paste(scntfcN, collapse = ", "), .groups = "drop")

paletteHash <- split(x = palette$col, f = palette$cat)
taxaHash <- split(x = taxa$taxa, f = taxa$MAP_LAB)

vascularData <- list(palette = paletteHash, taxa = taxaHash, mapTitle = "Map 1. Ecological communities", regionLabels = becLabels)

taxa <- base::merge(taxa, palette, by.x ="MAP_LAB", by.y="cat")

# Write summarised plants to JSON file for viz

write(jsonlite::toJSON(vascularData, auto_unbox = TRUE, pretty = TRUE), "viz_data/Vascular-plotData.json")

# Plot map

speciesMap <- leaflet(options=list(mx_mapId="Vascular")) %>%
  setView(-123.2194, 49.66076, zoom = 8.5) %>%
  addTiles(options = providerTileOptions(opacity = 0.5)) %>%
  addRasterImage(hillshade, opacity = 0.8) %>%
  addPolygons(data = coastline, color = "black", weight = 1.5, fill = FALSE) %>%
  addPolygons(data = plants.x.BEC, fillColor = plants.x.BEC$col, fillOpacity = 0.6, weight = 2, stroke = FALSE, options = labelToOption(BEC$MAP_LABEL)) %>%
  addPolygons(data = watershed.boundary, color = "black", weight = 4, fill = FALSE)

#Note that this statement is only effective in standalone R
print(speciesMap)
