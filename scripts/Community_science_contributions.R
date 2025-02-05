# Map reporting status of vascular plants documented in the Átl’ka7tsem/Howe Sound Biosphere

# Load libraries

library(dplyr)
library(gapminder)
library(gganimate)
library(ggplot2)
library(ggthemes)
library(gifski)
library(hrbrthemes)
library(leaflet)
library(raster)
library(sf)
library(stringr)
library(tidyr)
library(viridis)
library(jsonlite)

library(plotly)

# Set relative paths (https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location-in-rstudio)
if (!isTRUE(getOption('knitr.in.progress'))) {
  setwd(paste0(dirname(rstudioapi::getActiveDocumentContext()$path), "/.."))
}

# Source dependencies

source("scripts/utils.R")
source("scripts/geomUtils.R")

# Read gridded records and summary

reporting.status.grid <- read.csv("tabular_data/gridded_reporting_status.csv")

howegrid <- read_grid_frame("tabular_data/gridframe.json")

fetchLayer <- function (requiredStatus) {
  rows <- filter(reporting.status.grid, status == requiredStatus)
  shape <- assign_cell_geometry_sf(rows, howegrid)
}

gridded.confirmed.records <- fetchLayer("confirmed")
gridded.historical.records <- fetchLayer("historical")
gridded.new.records <- fetchLayer("new")

# Create color palette for species richness

richness <- reporting.status.grid$richness
values <- richness %>% unique
values <- sort(values)
t <- length(values)
pal <- leaflet::colorFactor(viridis_pal(option = "D")(t), domain = values)

# Plot map

reportingStatusMap <- leaflet(options=list(mx_mapId="Status")) %>%
  setView(-123.2194, 49.66076, zoom = 8.5) %>%
  addTiles(options = providerTileOptions(opacity = 0.5)) %>%
  addRasterImage(hillshade, opacity = 0.8) %>%
  addPolygons(data = coastline, color = "black", weight = 1.5, fillOpacity = 0, fillColor = NA) %>%
  addPolygons(data = gridded.confirmed.records, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = labelToOption("confirmed")) %>%
  addPolygons(data = gridded.historical.records, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = labelToOption("historical")) %>%
  addPolygons(data = gridded.new.records, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = labelToOption("new")) %>%
  addLegend(position = 'topright',
            colors = viridis_pal(option = "D")(t),
            labels = values) %>%
  addPolygons(data = watershed.boundary, color = "black", weight = 4, fillOpacity = 0)

#Note that this statement is only effective in standalone R
print(reportingStatusMap)

# Stacked bar plot of historic vs confirmed vs new records

summary <- read.csv("tabular_data/vascular_plant_summary_resynthesized_2024-11-14.csv", fileEncoding = "UTF-8")

# Simplify "new" values of reportingStatus
summary <- summary %>%
  mutate(reportingStatus = ifelse(str_detect(reportingStatus, "new"), "new", reportingStatus))

summary <- summary %>%
  mutate(reportingStatus = ifelse(reportingStatus == "reported", "historical", reportingStatus))

new <- summary %>% filter(reportingStatus == "new")
confirmed <- summary %>% filter(reportingStatus == "confirmed")
historical <- summary %>% filter(reportingStatus == "historical")

y <- c('records')
confirmed.no <- c(nrow(confirmed))
historical.no <- c(nrow(historical))
new.no <- c(nrow(new))

reporting.status <- data.frame(y, confirmed.no, historical.no, new.no)

reportingStatusFig <- plot_ly(height = 140, reporting.status, x = ~confirmed.no, y = ~y, type = 'bar', orientation = 'h', name = 'confirmed',
                      
                      marker = list(color = '#5a96d2',
                             line = list(color = '#5a96d2',
                                         width = 1)))

# These names need to agree with those applied as labels to the map regions
reportingStatusFig <- reportingStatusFig %>% add_trace(x = ~historical.no, name = 'historical',
                         marker = list(color = '#decb90',
                                       line = list(color = '#decb90',
                                                   width = 1)))
reportingStatusFig <- reportingStatusFig %>% add_trace(x = ~new.no, name = 'new',
                         marker = list(color = '#7562b4',
                                       line = list(color = '#7562b4',
                                                   width = 1)))
reportingStatusFig <- reportingStatusFig %>% layout(barmode = 'stack', autosize=T, height=140, showlegend=FALSE,
                      xaxis = list(title = "Species Reporting Status"),
                      yaxis = list(title = "Records")) %>% 
  layout(meta = list(mx_widgetId = "reportingStatus")) %>%
  layout(yaxis = list(showticklabels = FALSE)) %>%
  layout(yaxis = list(title = "")) %>%
  config(displayModeBar = FALSE, responsive = TRUE)

reportingStatusFig

reportingPal <- list("confirmed" = "#5a96d2", "historical" = "#decb90", "new" = "#7562b4")

taxa.status <- summary %>% group_by(reportingStatus) %>% 
  summarize(taxa = paste(sort(unique(scientificName)),collapse=", "))

# Convert taxa to named list so that JSON can recognise it
statusTaxa <- split(x = taxa.status$taxa, f=taxa.status$reportingStatus)

# Write summarised plants to JSON file for viz 
# (selection states corresponding with bar plot selections: 'new', 'historical','confirmed')
statusData <- structure(list(palette = reportingPal, taxa = statusTaxa, mapTitle = "Map 3. Species Reporting Status"))

jsonStatus = jsonlite::toJSON(statusData, auto_unbox = TRUE, pretty = TRUE)

write_utf8(jsonStatus, "viz_data/Status-plotData.json")

# Export CSVs for confirmed, historical and new reports

plants <- timedFread("tabular_data/Howe_Sound_vascular_plant_records_consolidated_2024-11-14.csv")

confirmed.taxa.records <- plants %>% filter(scientificName %in% confirmed$scientificName)
new.taxa.records <- plants %>% filter(scientificName %in% new$scientificName)
historical.taxa.records <- plants %>% filter(scientificName %in% historical$scientificName)

timedWrite(confirmed.taxa.records, "outputs/AHSBR_vascular_plants_confirmed_taxa_records.csv")

timedWrite(new.taxa.records, "outputs/AHSBR_vascular_plants_new_taxa_records.csv")

timedWrite(historical.taxa.records, "outputs/AHSBR_vascular_plants_historical_taxa_records.csv")
