# Map history of vascular plant surveys in Átl’ka7tsem

# Set relative paths (https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location-in-rstudio)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

# Load libraries

library(dplyr)
library(gapminder)
library(gganimate)
library(ggplot2)
library(ggthemes)
library(gifski)
library(hrbrthemes)
library(jsonlite)
library(leaflet)
library(plotly)
library(raster)
library(sf)
library(tidyr)
library(viridis)

# Set relative paths (https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location-in-rstudio)
if (!isTRUE(getOption('knitr.in.progress'))) {
  setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/.."))
}

# Source dependencies

source("scripts/utils.R")
source("scripts/geomUtils.R")

cumulative.history <- read.csv("tabular_data/history_1897-2024_cumulative.csv")

year <- cumulative.history$year

# Plot analysis of historical collection activities

speciesPlot <- plot_ly(width = 700, height = 420)

speciesPlot <- speciesPlot %>%
      layout(title='Vascular plant species recorded in Átl’ḵa7tsem/Howe Sound 1890-2024', showlegend = FALSE,
      xaxis = list(title="Year", range = c(1900, 2024)),
      yaxis = list(title='Reported Species', range=c(0, max(cumulative.history$cum.spp)))
      )

steps <- list()

# General method cribbed from https://plotly.com/r/sliders/#sine-wave-slider

for (i in 1:length(year)) {
    args <- list('visible', rep(FALSE, length(year)))
    args[[2]][i] = TRUE
    steps[[i]] <- list(label = year[[i]], method="restyle", args=args)
    sppRange <- cumulative.history$cum.spp[1:i]
    yearRange <- year[1:i]
    speciesPlot <- speciesPlot %>% add_lines(x=yearRange, y=sppRange, line=list(color='green'), type='scatter', mode='lines', visible = i == 1)
}

speciesPlot <- layout(speciesPlot, meta = list(mx_widgetId = "speciesPlot"), 
                      sliders = list(list(active=0, steps = steps)))

speciesPlot


# Most historical data collected in the 1920s and between the 1960s and 1980s;
# Large increase in observations and recorded species with the emergence of 
# iNaturalist in 2010s; no. species reported for Howe Sound has nearly doubled 
# over the last two decades

# Plot gridded choropleth illustrating historical timeline of plant surveys 1897-2022

# Note: this plot adds one layer, including cumulative series of features for multiple decades
# to illustrate historical timeline;
# alternatively, we might add multiple layers of gridded species richness, one for each decade

# Load map layers

# Layer 1: hillshade raster
hillshade <- raster("spatial_data/rasters/Hillshade_80m.tif")

# Layer 2: coastline
coastline <- mx_read("spatial_data/vectors/Islands_and_Mainland")

# Layer 3: watershed boundary
watershed.boundary <- mx_read("spatial_data/vectors/Howe_Sound")

# Layer 4: gridded history choropleth
gridded.history <- read.csv("tabular_data/gridded_history_1897-2024_cumulative.csv")
howegrid <- read_grid_frame("tabular_data/gridframe.json")

# Create color palette for species richness

richness <- gridded.history$richness
values <- richness %>% unique
values <- sort(values)
t <- length(values)
pal <- leaflet::colorFactor(viridis_pal(option = "D")(t), domain = values)

fetchLayer <- function (year) {
    rows <- filter(gridded.history, period == paste0("1890-", year))
    shape <- assign_cell_geometry_sf(rows, howegrid)
}

history.1900 <- fetchLayer(1900)
history.1910 <- fetchLayer(1910)
history.1920 <- fetchLayer(1920)
history.1930 <- fetchLayer(1930)
history.1940 <- fetchLayer(1940)
history.1950 <- fetchLayer(1950)
history.1960 <- fetchLayer(1960)
history.1970 <- fetchLayer(1970)
history.1980 <- fetchLayer(1980)
history.1990 <- fetchLayer(1990)
history.2000 <- fetchLayer(2000)
history.2010 <- fetchLayer(2010)
history.2020 <- fetchLayer(2020)
history.2024 <- fetchLayer(2024)

historyData <- list(mapTitle = "Map 2. Historical collection activities")
write_json(historyData, "viz_data/History-plotData.json")

# Plot map

heatMap <- leaflet(options=list(mx_mapId="History")) %>%
  setView(-123.2194, 49.66076, zoom = 8.5) %>%
  addTiles(options = providerTileOptions(opacity = 0.5)) %>%
  addRasterImage(hillshade, opacity = 0.8) %>%
  addPolygons(data = coastline, color = "black", weight = 1.5, fillOpacity = 0, fillColor = NA) %>%
  # These indices need to agree with positions in array "year"
  addPolygons(data = history.1900, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = list(mx_subLayerIndex = 0)) %>%
  addPolygons(data = history.1910, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = list(mx_subLayerIndex = 1)) %>%
  addPolygons(data = history.1920, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = list(mx_subLayerIndex = 2)) %>%
  addPolygons(data = history.1930, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = list(mx_subLayerIndex = 3)) %>%
  addPolygons(data = history.1940, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = list(mx_subLayerIndex = 4)) %>%
  addPolygons(data = history.1950, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = list(mx_subLayerIndex = 5)) %>%
  addPolygons(data = history.1960, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = list(mx_subLayerIndex = 6)) %>%
  addPolygons(data = history.1970, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = list(mx_subLayerIndex = 7)) %>%
  addPolygons(data = history.1980, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = list(mx_subLayerIndex = 8)) %>%
  addPolygons(data = history.1990, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = list(mx_subLayerIndex = 9)) %>%
  addPolygons(data = history.2000, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = list(mx_subLayerIndex = 10)) %>%
  addPolygons(data = history.2010, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = list(mx_subLayerIndex = 11)) %>%
  addPolygons(data = history.2020, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = list(mx_subLayerIndex = 12)) %>%
  addPolygons(data = history.2024, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = list(mx_subLayerIndex = 13)) %>%
  addLegend(position = 'topright',
            colors = viridis_pal(option = "D")(t),
            labels = values) %>%
  addPolygons(data = watershed.boundary, color = "black", weight = 4, fillOpacity = 0)

#Note that this statement is only effective in standalone R
print(heatMap)


