# Map Átl’ka7tsem's vascular plant diversity in relation to protected areas

# Set relative paths (https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location-in-rstudio)
if (!isTRUE(getOption('knitr.in.progress'))) {
  setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/.."))
}

# Load libraries

library(Dict)
library(dplyr)
library(leaflet)
library(plotly)
library(raster)
library(reshape2)
library(sf)
library(stringr)
library(jsonlite)
library(viridis)

# Source dependencies

source("scripts/utils.R")

# Read occurrence data (plants x protected area)

plants.x.protected.area <- read.csv("tabular_data/plants_x_protected_areas_2024.csv")

plants.x.protected.area$X <- NULL

# Load Protected Areas Shape

protected.areas <- mx_read("spatial_data/vectors/Protected_Areas")

# Create labels

protected.areas$prtct__[is.na(protected.areas$prtct__)] <- 0

protected.areas$label <- paste(protected.areas$prtctdA, ":", protected.areas$prtct__, "species", sep = " ")

# Load additional map layers

# Layer 1: hillshade raster
hillshade <- raster("spatial_data/rasters/Hillshade_80m.tif")

# Layer 2: coastline
coastline <- mx_read("spatial_data/vectors/Islands_and_Mainland")

# Layer 3: watershed boundary
watershed.boundary <- mx_read("spatial_data/vectors/Howe_Sound")

# Note: for the following map, species can be mapped to protected areas based on the 
# catalog called as 'plants.x.protected.area' 

# Create color palette

palette = data.frame(
  types = unique(protected.areas$prtctAT),
  colors = c('#8b4513', '#008000', '#4682b4', '#4b0082', '#ff0000', '#00ff00', '#00ffff', '#0000ff', '#ffff54', '#ff69b4', '#ffe4c4')
)

protected.areas$colors <- palette$colors[match(unlist(protected.areas$prtctAT), palette$types)]

protectedData <- list(mapTitle = "Map 4. Protected Areas")
write(jsonlite::toJSON(protectedData, auto_unbox = TRUE, pretty = TRUE), "viz_data/Protected-plotData.json")

# Plot map

protectedAreaMap <- leaflet(options=list(mx_mapId="Protected")) %>%
  setView(-123.2194, 49.66076, zoom = 8.5) %>%
  addTiles(options = providerTileOptions(opacity = 0.5)) %>%
  addRasterImage(hillshade, opacity = 0.8) %>%
  addPolygons(data = coastline, color = "black", weight = 1.5, fillOpacity = 0, fillColor = NA) %>%
  addPolygons(data = watershed.boundary, color = "black", weight = 4, fillOpacity = 0) %>% 
  addPolygons(data = protected.areas, fillColor = protected.areas$colors, fillOpacity = 0.8, weight = 0,
            label = paste(protected.areas$prtctdA, ":", protected.areas$prtct__, "species", sep = " "))

#Note that this statement is only effective in standalone R
print(protectedAreaMap)

# Create dataframe summarizing plant diversity by protected area type

types <- as.factor(unique(protected.areas$prtctAT))
count <- vector(mode="numeric", length=length(types))

protected.area.summary <- data.frame(types,count)

protected.area.summary$count <- protected.areas$prtctd_r__[match(unlist(protected.area.summary$types), protected.areas$prtctAT)]

protected.area.summary <- protected.area.summary[order(protected.area.summary$types),]

protected.area.summary$types <- factor(protected.area.summary$types, levels = unique(protected.area.summary$types)[order(protected.area.summary$count, decreasing = TRUE)])


# Create Plotly bar plot showing species diversity represented within protected area types

# First add color palette matching with map

protected.area.summary$colors <- palette$colors[match(unlist(protected.area.summary$types), palette$types)]

colormap <- setNames(object = protected.area.summary$colors,
                     nm = protected.area.summary$types)

# Plot

protected.area.plot <- plot_ly(
                          data = protected.area.summary,
                          y = ~types,
                          x = ~count,
                          color = ~types,
                          colors = colormap,
                          opacity = 0.8,
                          type = "bar"
                            ) %>% layout(xaxis = list(categoryorder = "category ascending", title = "Species Reported by Protected Area")) %>%
                                  layout(yaxis = list(title = "", width=1024)) %>%
                                  layout(meta = list(mx_widgetId = "protectedAreas"))

protected.area.plot 

# Prepare CSVs to export catalog data for each Protected Area Type

OECM.plants <- plants.x.protected.area %>% filter(protectedAreaType == 'Other Effective Area-Based Conservation Measure')
A.Park.plants <- plants.x.protected.area %>% filter(protectedAreaType == 'A - Park')
Private.Conservation.Area.plants <- plants.x.protected.area %>% filter(protectedAreaType == 'Privately Owned Conservation Area')
WHA.plants <- plants.x.protected.area %>% filter(protectedAreaType == 'Wildlife Habitat Areas')
PA.plants <- plants.x.protected.area %>% filter(protectedAreaType == 'Protected Area')
S2S.plants <- plants.x.protected.area %>% filter(protectedAreaType == 'Sea To Sky Wildland Zones')
MBS.plants <- plants.x.protected.area %>% filter(protectedAreaType == 'Migratory Bird Sanctuary')
Conservancy.plants <- plants.x.protected.area %>% filter(protectedAreaType == 'Conservancy')
ER.plants <- plants.x.protected.area %>% filter(protectedAreaType == 'Ecological Reserve')
OGMA.plants <- plants.x.protected.area %>% filter(protectedAreaType == 'Old Growth Management Areas (Mapped Legal)')
WMA.plants <- plants.x.protected.area %>% filter(protectedAreaType == 'Wildlife Management Area')

timedWrite(OECM.plants, "outputs/AHSBR_OECM_vascular_plants.csv")
timedWrite(A.Park.plants, "outputs/AHSBR_Class_A_Parks_vascular_plants.csv")
timedWrite(Private.Conservation.Area.plants, "outputs/AHSBR_Private_Conservation_Area_vascular_plants.csv")
timedWrite(WHA.plants, "outputs/AHSBR_Wildlife_Habitat_Area_vascular_plants.csv")
timedWrite(PA.plants, "outputs/AHSBR_Protected_Area_vascular_plants.csv")
timedWrite(S2S.plants, "outputs/AHSBR_Sea_to_Sky_Wildland_Zones_vascular_plants.csv")
timedWrite(MBS.plants, "outputs/AHSBR_Migratory_Bird_Sanctuary_vascular_plants.csv")
timedWrite(Conservancy.plants, "outputs/AHSBR_Conservancy_vascular_plants.csv")
timedWrite(ER.plants, "outputs/AHSBR_Ecological_Reserve_vascular_plants.csv")
timedWrite(OGMA.plants, "outputs/AHSBR_Old_Growth_Management_Areas_vascular_plants.csv")
timedWrite(WMA.plants, "outputs/AHSBR_Wildlife_Management_Area_vascular_plants.csv")


## Prepare summary of native plant diversity protected in the area

plant.summary <- timedFread("tabular_data/vascular_plant_summary_resynthesized_2024-11-14.csv");

native.plants <- plant.summary %>% filter(establishmentMeans == 'native')

protected.plants <- unique(plants.x.protected.area$scientificName)

protected.plants <- protected.plants %>% paste(collapse = '|')
  
protected.plants <- plant.summary %>% filter(str_detect(scientificName, protected.plants))

timedWrite(native.plants, "outputs/AHSBR_native_vascular_plant_species.csv")

timedWrite(protected.plants, "outputs/AHSBR_protected_native_vascular_plant_species.csv")