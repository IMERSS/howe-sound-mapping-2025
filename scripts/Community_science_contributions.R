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

# Read records and summary

plants <- read.csv("tabular_data/gridded_plants_WGS84.csv")

summary <- read.csv("tabular_data/vascular_plant_summary_resynthesized_2023-03-05.csv")

# Subset historic, confirmed and new records

new <- summary %>% filter(str_detect(reportingStatus, "new"))
confirmed <- summary %>% filter(reportingStatus == "confirmed")
reported <- summary %>% filter(reportingStatus == "reported")

# Create vectors of historic, confirmed and new taxa to query catalog of occurrence records

new.taxa <- unique(new$scientificName)
new.taxa <- new.taxa %>% paste(collapse = "|")

confirmed.taxa <- unique(confirmed$scientificName)
confirmed.taxa <- confirmed.taxa %>% paste(collapse = "|")

reported.taxa <- unique(reported$scientificName)
reported.taxa <- reported.taxa %>% paste(collapse = "|")

new.taxa.records <- plants %>% filter(str_detect(scientificName, new.taxa))
new.taxa.records$status <- 'new'

confirmed.taxa.records <- plants %>% filter(str_detect(scientificName,confirmed.taxa))
confirmed.taxa.records$status <- 'confirmed'

reported.taxa.records <- plants %>% filter(str_detect(scientificName, reported.taxa))
reported.taxa.records$status <- 'historic'

records <- rbind(new.taxa.records, confirmed.taxa.records, reported.taxa.records)

# Summarise plant species by reporting status

taxa.status <- records %>% group_by(status) %>% 
  summarize(taxa = paste(sort(unique(scientificName)),collapse=", "))

### PLOT GRIDDED HEATMAPS 

# Load map layers

# Hillshade raster
hillshade <- raster("spatial_data/rasters/Hillshade_80m.tif")

# Coastline
coastline <- mx_read("spatial_data/vectors/Islands_and_Mainland")

# Watershed boundary
watershed.boundary <- mx_read("spatial_data/vectors/Howe_Sound")

# Gridded records by reporting status
gridded.confirmed.records <- mx_read("spatial_data/vectors/gridded_confirmed_records")
gridded.historic.records <- mx_read("spatial_data/vectors/gridded_historic_records")
gridded.new.records <- mx_read("spatial_data/vectors/gridded_new_records")

# Combine records to create normalized palette

gridded.records <- rbind(gridded.confirmed.records, gridded.historic.records, gridded.new.records)

# Create color palette for species richness

richness <- gridded.records$richness
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
  addPolygons(data = gridded.historic.records, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = labelToOption("historic")) %>%
  addPolygons(data = gridded.new.records, fillColor = ~pal(richness), fillOpacity = 0.6, weight = 0, options = labelToOption("new")) %>%
  addLegend(position = 'topright',
            colors = viridis_pal(option = "D")(t),
            labels = values) %>%
  addPolygons(data = watershed.boundary, color = "black", weight = 4, fillOpacity = 0)

#Note that this statement is only effective in standalone R
print(reportingStatusMap)

# Stacked bar plot of historic vs confirmed vs new records

y <- c('records')
confirmed.no <- c(nrow(confirmed))
historic.no <- c(nrow(reported))
new.no <- c(nrow(new))

reporting.status <- data.frame(y, confirmed.no, historic.no, new.no)

reportingStatusFig <- plot_ly(height = 140, reporting.status, x = ~confirmed.no, y = ~y, type = 'bar', orientation = 'h', name = 'confirmed',
                      
                      marker = list(color = '#5a96d2',
                             line = list(color = '#5a96d2',
                                         width = 1)))

# These names need to agree with those applied as labels to the map regions
reportingStatusFig <- reportingStatusFig %>% add_trace(x = ~historic.no, name = 'historic',
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

reportingPal <- list("confirmed" = "#5a96d2", "historic" = "#decb90", "new" = "#7562b4")

# Convert taxa to named list so that JSON can recognise it
statusTaxa <- split(x = taxa.status$taxa, f=taxa.status$status)

# Write summarised plants to JSON file for viz 
# (selection states corresponding with bar plot selections: 'new', 'historic','confirmed')
statusData <- structure(list(palette = reportingPal, taxa = statusTaxa, mapTitle = "Map 3. Species Reporting Status"))

write(jsonlite::toJSON(statusData, auto_unbox = TRUE, pretty = TRUE), "viz_data/Status-plotData.json")

# Export CSVs for confirmed, historic and new reports

plants <- read.csv("tabular_data/Howe_Sound_vascular_plant_records_consolidated.csv")

confirmed.taxa.records <- plants %>% filter(str_detect(scientificName, confirmed.taxa))
new.taxa.records <- plants %>% filter(str_detect(scientificName, new.taxa))
reported.taxa.records <- plants %>% filter(str_detect(scientificName, reported.taxa))

write.csv(confirmed.taxa.records, "outputs/AHSBR_vascular_plants_confirmed_taxa_records.csv", row.names = FALSE)

write.csv(new.taxa.records, "outputs/AHSBR_vascular_plants_new_taxa_records.csv", row.names = FALSE)

write.csv(reported.taxa.records, "outputs/AHSBR_vascular_plants_historic_taxa_records.csv", row.names = FALSE)
