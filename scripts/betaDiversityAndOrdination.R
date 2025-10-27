# Exploratory analysis of Átl'ka7tsem/Howe Sound vascular plant diversity

 # Implementing gridded beta diversity analysis of Howe Sound vegetation data based on:
 # https://rfunctions.blogspot.com/2015/08/calculating-beta-diversity-on-grid.html

 # Set relative paths (https://stackoverflow.com/questions/13672720/r-command-for-setting-working-directory-to-source-file-location-in-rstudio)

 setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) 

 # Load packages

 library(betapart)
 library(CommEcol)
 library(dplyr)
 library(ecodist)
 library(gapminder)
 library(gganimate)
 library(ggplot2)
 library(ggrepel)
 library(ggthemes)
 library(gifski)
 library(hrbrthemes)
 library(picante)
 library(raster)
 library(RColorBrewer)
 library(rgdal)
 library(rgeos)
 library(stringr)
 library(terra)
 library(tidyr)
 library(vegan)


 #### FUNCTION: copy and paste into R ####

 betagrid<-function(gridshp, comp, xfeature, yfeature, radius, phylotree, phylobeta=F, index="sorensen"){
   data<-data.frame(gridshp[xfeature],gridshp[yfeature],comp)
   mean_turnover<-numeric(length(comp[,1]))
   mean_nestedness<-numeric(length(comp[,1]))
   mean_beta<-numeric(length(comp[,1]))
   for(i in 1:length(gridshp[[2]])){
     adj<-select.window(xf=data[i,1], yf=data[i,2], radius, xydata=data)[,-c(1,2)]
     if(phylobeta==F){
       ifelse(sum(nrow(adj))==0 || ncol(adj)==0, res<-0 , res<-beta.pair(adj, index.family=index))
     }else if(phylobeta==T){
       ifelse(sum(nrow(adj))==0 || ncol(adj)==0, res<-0 , res<-phylo.beta.pair(adj, phylotree, index.family=index))
     }
     ifelse(sum(nrow(adj))==0 || ncol(adj)==0, mean_turnover[i]<-0 , mean_turnover[i]<-mean(as.matrix(res[[1]])[2:length(as.matrix(res[[1]])[,1]),1],na.rm=TRUE) )
     ifelse(sum(nrow(adj))==0 || ncol(adj)==0, mean_nestedness[i]<-0 , mean_nestedness[i]<-mean(as.matrix(res[[2]])[2:length(as.matrix(res[[2]])[,1]),1],na.rm=TRUE) )
     ifelse(sum(nrow(adj))==0 || ncol(adj)==0, mean_beta[i]<-0 , mean_beta[i]<-mean(as.matrix(res[[3]])[2:length(as.matrix(res[[3]])[,1]),1],na.rm=TRUE) )
   }
   return(data.frame(cell=row.names(comp), mean_turnover, mean_nestedness, mean_beta))
 }


 #############
 ## MY DATA ##
 #############

 # Adapting the above code to implement gridded beta diversity analysis of my data:

 # Load the grid

 shape <- readOGR("spatial_data/vectors/1km_grid_WGS84_coordinates_x_vascular_plant_grid_NA_omit")

 # Read species occurrences

 data <- read.csv("tabular_data/1km_gridded_vascular_plant_records_2022-12-24_WGS84.csv")

 # Subset relevant fields

 data <- data %>% dplyr::select('scientific'|'id')

 # Add count field to generate matrix

 data$count <- 1

 # Generate matrix 

 matrix <- ecodist::crosstab(data$id, data$scientific, data$count)

 # Convert to presence / absence

 matrix[matrix > 0] <- 1 

 # Compare dimensions of matrix and grid cells

 nrow(matrix) # matrix with 771 rows matching 771 grid cells found in shape below
 nrow(shape)  # shape with 771 grid cells 

 # Which fields correspond with LONG (3) & LAT (2)? 

 names(shape)
 str(shape)

 # Call the function and get results! Let us calculate beta diversity for each focal cell. Note that the function will return results containing four columns: number of grid cell, the mean turnover partition of beta diversity, the mean nestedness partition of beta diversity, and the mean total beta diversity. Also, note that radius equals 0.25 degree, which is the same size as the resolution of our grid. This will make the function use only the 8 (or fewer) adjacent cells in relation to the focal cells. If you want more neighbor cells to be included in the analysis, you can use the double (0.5 in this example) or greater values.

 results <- betagrid(gridshp=shape, comp=matrix, xfeature=3, yfeature=2, radius=0.25, index="sorensen")

 # Standardize results fields to merge with grid in QGIS

 names(results) <- c('id','mean_turnover','mean_nestedness','mean_beta')

 # Add species richness to results

 matrix <- matrix %>% mutate(richness = rowSums(.))

 matrix$richness

 matrix$id <- row.names(matrix) 

 results$richness  <- matrix$richness[match(unlist(results$id), matrix$id)]

 # Output results

 write.csv(results,"outputs/betagrid_vascular_plants.csv", row.names = FALSE)


 ## Ordination of gridded vascular plant data

 # Read environmental metadata

 metadata <- read.csv("tabular_data/1km_grid_metadata.csv")

 # Construct label using environmental metadata to index grid cells in ordination analysis

 # Assign BEC unit based on spatial representation of BEC units within grid cells

 # First simplify BEC unit labels

 bec.labels <- unique(metadata$BGC_LABEL)
 bec.labels.new <- c("A","Mmm","Mmm","Hvm","Hvm","Hdm","Hxm","Cmm","Hms","Hds","Emw")

 grid.bec.labels <- data.frame(bec.labels,bec.labels.new)

 metadata$BGC_LABEL_NEW <- grid.bec.labels$bec.labels.new[match(unlist(metadata$BGC_LABEL), grid.bec.labels$bec.labels)]

 # Identify most well represented BEC units in each grid cell:

 grid.bec <- metadata %>% group_by(id, BGC_LABEL_NEW) %>% summarise(HA = max(HA))

 grid.bec <- grid.bec %>% group_by(id) %>% top_n(1, HA)

 # Reassign BEC labels

 metadata$BGC_LABEL_NEW <- grid.bec$BGC_LABEL_NEW[match(unlist(metadata$id), grid.bec$id)]

 # Concatenate grid cell labels

 metadata$GRID_ID <- paste(metadata$id, metadata$BGC_LABEL_NEW, metadata$REGION, sep = "", collapse = NULL)

 # Read betadiversity analysis results

 betadiversity <- read.csv("outputs/betagrid_vascular_plants.csv")

 # Merge betadiversity analysis results with metadata (note: betadiversity information only available for cells containing species)

 metadata <- left_join(metadata,betadiversity)

 # Reduce metadata frame to grid cells with biodiversity present

 metadata <- metadata %>% drop_na(richness)

 # Filter grid cells with species richness >30 (same is done with species dataframe below)

 metadata <- metadata %>% filter(richness>30)

 # How well are different biogeoclimatic units represented in this analysis?

 metadata %>% count(MAP_LABEL)

 # Table results: # grid cell representation of BEC units:

 # Coastal Western Hemlock dry maritime (48) - CWHdm

 # Coastal Western Hemlock (eastern) very dry maritime (29) - CWHxm1

 # Coastal Western Hemlock (montane) very wet maritime (12) - CWHvm2

 # Mountain henmlock (windward) moist maritime (10) - MHmm1

 # Alpine Tundra and Subalpine Parkland (7) - CMAunp

 # Mountain henmlock (leeward) moist maritime (5) - MHmm2

 # Coastal Western Hemlock (southern) dry sub-maritime (2) - CWHds1

 # Coastal Western Hemlock (southern)  moist maritime (2) - CWHms1


 # Read species occurrences

 data <- read.csv("tabular_data/1km_gridded_vascular_plant_records_2022-12-24_WGS84.csv")

 # Remove outlier grid cell (vascular plant records from early 1900s with problematic generalized coordinates)

 data <- data[!(data$id == '2716'),]

 # Subset relevant fields from biodiversity data

 data <- data %>% dplyr::select('scientific'|'id')

 # Add count field to generate matrix

 data$count <- 1

 # Generate matrix 

 matrix <- ecodist::crosstab(data$id, data$scientific, data$count)

 # Convert to presence / absence

 matrix[matrix > 0] <- 1 

 # Focus analysis on grid cells with species richness >30

 matrix <- matrix[rowSums(matrix[])>=30,]



 # Implement ordination analysis of grid cells with species richness >30

 # ( To consider ordinations against environmental variables, see:
 # http://kembellab.ca/r-workshop/biodivR/SK_Biodiversity_R.html )
 # OR : # https://www.rpubs.com/RGrieger/545184

 # Create MDS Jaccard similarity Index for Site Commmunity Data

 AHSBR.MDS <- metaMDS(matrix, distance = "jaccard", k = 2, trymax = 100, zerodist = "add")

 scores(matrix.MDS)

 stressplot(matrix.MDS)

 plot(matrix.MDS, type = 't')

 # Save NMDS results into dataframe

 site.scrs <- as.data.frame(scores(AHSBR.MDS, display = "sites")) 

 site.scrs <- cbind(site.scrs, id = rownames(site.scrs)) #add site names as variable if you want to display on plot

 site.scrs$id <- as.integer(site.scrs$id)

 # Bind NMDS scores with environmental metadata

 AHSBR.MDS.results <- left_join(metadata, site.scrs, by = 'id')

 # Prepare dataset for NMDS species scores

 AHSBR.spp.fit <- envfit(AHSBR.MDS, matrix, permutations = 999) # this fits species vectors

 spp.scrs <- as.data.frame(scores(AHSBR.spp.fit, display = "vectors")) #save species intrinsic values into dataframe
 spp.scrs <- cbind(spp.scrs, Species = rownames(spp.scrs)) #add species names to dataframe
 spp.scrs <- cbind(spp.scrs, pval = AHSBR.spp.fit$vectors$pvals) #add pvalues to dataframe so you can select species which are significant
 spp.scrs<- cbind(spp.scrs, abrev = abbreviate(spp.scrs$Species, minlength = 6)) #abbreviate species names
 sig.spp.scrs <- subset(spp.scrs, pval<=0.05) #subset data to show species significant at 0.05

 # Classic nMDS Plot

 ordiplot(matrix.MDS,type="p")
 ordiplot (matrix.MDS, display = 'sites', type = 'p')

 # Plot using GGPlot

 AHSBR.MDS.plot <- ggplot(AHSBR.MDS.results, aes(x=NMDS1, y=NMDS2, label=id))+ #sets up the plot
   geom_point(aes(NMDS1, NMDS2, colour = factor(AHSBR.MDS.results$MAP_LABEL), shape = factor(AHSBR.MDS.results$ZONE)), size = 2)+ #adds site points to plot
   geom_text(hjust=0, vjust=0)+
   coord_fixed()+
   theme_classic()+ 
   theme(panel.background = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"))+
   labs(colour = "BEC unit", shape = "ZONE")+ # add legend labels for Management and Landuse
   theme(legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12), axis.text = element_text(size = 10)) # add legend at right of plot

 AHSBR.MDS.plot + labs(title = "Átl'ka7tsem/Howe Sound Biosphere Vascular Plant Communities") #displays plot

 # With significant species:

 AHSBR.MDS.plot +
   geom_segment(data = sig.spp.scrs, aes(x = 0, xend=NMDS1, y=0, yend=NMDS2), arrow = arrow(length = unit(0.25, "cm")), colour = "grey10", lwd=0.3) + #add vector arrows of significant species
   ggrepel::geom_text_repel(data = sig.spp.scrs, aes(x=NMDS1, y=NMDS2, label = Species), cex = 3, direction = "both", segment.size = 0.25)+ #add labels for species, use ggrepel::geom_text_repel so that labels do not overlap
   labs(title = "Ordination with species vectors")