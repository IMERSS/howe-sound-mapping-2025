library(sf)
library(jsonlite)

# Generic utilities for working with geometric primitives including creating square gridded coordinate framed indexed by cell
# Antranig Basman, 2023-2025

# From https://en.wikipedia.org/wiki/Longitude#Length_of_a_degree_of_longitude
hortis.WGS84a = 6378137;
hortis.WGS84b = 6356752.3142;
hortis.WGS84e2 = (hortis.WGS84a * hortis.WGS84a - hortis.WGS84b * hortis.WGS84b) / (hortis.WGS84a * hortis.WGS84a);

# Length in metres for a degree of longitude at given latitude
hortis.longitudeLength <- function (latitude) {
    latrad <- pi * latitude / 180;
    sinrad <- sin(latrad);
    return (pi * hortis.WGS84a * cos(latrad) / (180 * sqrt(1 - hortis.WGS84e2 * sinrad * sinrad)))
}

# Length in metres for a degree of latitude at given latitude
hortis.latitudeLength = function (latitude) {
    latrad <- pi * latitude / 180;
    sinrad <- sin(latrad);
    return (pi * hortis.WGS84a * (1 - hortis.WGS84e2) / (180 * (1 - hortis.WGS84e2 * sinrad * sinrad) ^ 1.5));
}

# Convert a distance measure in longitude degrees to the same physical length in latitude degrees, given a particular latitude
hortis.longToLat <- function (lng, lat) {
    longLength <- hortis.longitudeLength(lat)
    latLength <- hortis.latitudeLength(lat)
    return (lng * longLength / latLength);
}

#' Generate a grid frame given spatial features and a cell size
#'
#' This function calculates a grid frame for a set of spatial points based on a specified cell size. 
#' It computes the bounding box of the input features and determines the longitudinal and latitudinal 
#' grid cell dimensions, as well as the number of cells required in each dimension.
#'
#' @param points An `sf` object or spatial dataset containing points or polygons from which the bounding box is derived.
#' @param cellsize A numeric value representing the desired size of each grid cell, in meters.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{bbox}{A list with the bounding box coordinates (`xmin`, `ymin`, `xmax`, `ymax`) as numeric values.}
#'   \item{bbox_b}{The bounding box as an `sf`-compatible object.}
#'   \item{longsize}{The longitudinal size of each grid cell, in decimal degrees.}
#'   \item{latsize}{The latitudinal size of each grid cell, in decimal degrees.}
#'   \item{longcount}{The number of grid cells along the longitudinal axis.}
#'   \item{latcount}{The number of grid cells along the latitudinal axis.}
#' }
#' @examples
#' library(sf)
#' points <- st_sfc(st_point(c(0, 0)), st_point(c(1, 1)), crs = 4326)
#' grid_frame <- make_grid_frame(points, 1000)
#' print(grid_frame)
#'
#' @export
make_grid_frame <- function (points, cellsize) {
  bbox_b <- st_bbox(points) # order xmin, ymin, xmax, ymax
  # Undo insufferable wrapping as "named numbers"
  bbox <- list(xmin = as.numeric(bbox_b$xmin), ymin = as.numeric(bbox_b$ymin), 
               xmax = as.numeric(bbox_b$xmax), ymax = as.numeric(bbox_b$ymax))
  midlat <- (bbox$ymax + bbox$ymin) / 2
  longdeg <- hortis.longitudeLength(midlat)
  longsize <- round(cellsize / longdeg, 6)
  latsize <- round(hortis.longToLat(longsize, midlat), 6)
  
  longcount <- as.integer(ceiling((bbox$xmax - bbox$xmin) / longsize))
  latcount <-  as.integer(ceiling((bbox$ymax - bbox$ymin) / latsize))
  
  list(bbox = bbox, bbox_b = bbox_b, longsize = longsize, latsize = latsize, longcount = longcount, latcount = latcount)
}

#' Write a grid frame to a JSON file
#'
#' This function saves a grid frame as a JSON file, excluding the \code{bbox_b} field.
#'
#' @param gridframe A list containing grid frame information, typically produced by the 
#'   \code{\link{make_grid_frame}} function.
#' @param filename A string specifying the path where the JSON file should be saved.
#'
#' @return No return value. The function writes the grid frame to a file.
#'
write_grid_frame <- function (gridframe, filename) {
  # https://stackoverflow.com/a/50143588/1381443
  filtered = within(gridframe, rm(bbox_b))
  write(jsonlite::toJSON(filtered, auto_unbox = TRUE, pretty = TRUE), filename)
}

#' Read a grid frame from a JSON file
#'
#' This function loads a grid frame from a JSON file.
#'
#' @param filename A string specifying the path to the JSON file.
#'
#' @return A list containing the grid frame data.
#'
#' @export
read_grid_frame <- function (filename) {
  gridframe <- jsonlite::read_json(filename)
}

#' Determine the grid cell index for a given geographical point in a supplied grid frame.
#'
#' @param gridframe A list containing grid frame information, typically produced by the \code{\link{make_grid_frame}} function. 
#'   Must include \code{bbox}, \code{longsize}, \code{latsize}, and \code{longcount}.
#' @param long A numeric value representing the longitude of the point.
#' @param lat A numeric value representing the latitude of the point.
#'
#' @return An integer representing the 0-based cell index if the point lies within the grid frame's bounding box.
#' If the point is outside the bounding box, the function returns \code{NA} and issues a warning.
#'
#' @details
#' The function checks whether the input point lies within the bounding box defined by \code{gridframe}. 
#' If so, it calculates the column and row indices (0-based) of the grid cell containing the point. 
#' The cell index is calculated using the formula:
#' \deqn{\text{cell_id} = \text{row_index} \times \text{longcount} + \text{col_index}}
#'
#' @note
#' The grid cell indices are 0-based, meaning the top-left cell is indexed as 0.
#' 
#' @examples
#' gridframe <- list(
#'   bbox = list(xmin = 0, ymin = 0, xmax = 10, ymax = 10),
#'   longsize = 1,
#'   latsize = 1,
#'   longcount = 10
#' )
#' cell_id <- point_to_cell(gridframe, long = 5, lat = 5)
#' print(cell_id) #55
#'
#' # Example with a point outside the bounding box
#' cell_id <- point_to_cell(gridframe, long = 15, lat = 15)
#' # Warning is issued, and NA is returned.
#'
#' @export
point_to_cell <- function (gridframe, long, lat) {
  # Check if the point is within the bbox
  if (long < gridframe$bbox$xmin || long > gridframe$bbox$xmax ||
      lat < gridframe$bbox$ymin || lat > gridframe$bbox$ymax) {
    warning("Point (", long, ", ", lat, ") lies outside the bounding box ", format(gridframe$bbox))
    return (NA)
  }
  
  # Calculate the column index (0-based)
  col_index <- as.integer(floor((long - gridframe$bbox$xmin) / gridframe$longsize))
  
  # Calculate the row index (0-based)
  row_index <- as.integer(floor((gridframe$bbox$ymax - lat) / gridframe$latsize))
  
  # Calculate the cell_id as a 0-based index
  cell_id <- row_index * gridframe$longcount + col_index
}

#' Given an sf frame with POINT geometries, derive a cell_id column containing the gridded cell id
#'
#' This function assigns a `cell_id` to each point in an `sf` dataframe based on a provided grid frame.
#'
#' @param points An `sf` dataframe with `POINT` geometries.
#' @param gridframe A list containing grid frame information, typically produced by the 
#'   \code{\link{make_grid_frame}} function.
#'
#' @return The input `sf` dataframe with an additional \code{cell_id} column, representing 
#' the grid cell associated with each point.
#'
#' @export
assign_cell_id <- function (points, gridframe) {
  assign.start <- Sys.time()
  coords <- sf::st_coordinates(points)
  longs <- coords[, "X"]
  lats <- coords[, "Y"]
  # Apply the calculate_cell_id function to all points in the sf object
  points$cell_id <- mapply(point_to_cell, longs, lats,
                           MoreArgs = list(gridframe = gridframe))
  assign.end <- Sys.time()
  cat ("Assigned ", nrow(points), " points in ", (assign.end - assign.start), "s\n")
  return (points)
}

#' Generate a `POLYGON` geometry representing the boundaries of a grid cell given its cell ID and a grid frame.
#'
#' @param gridframe A list containing grid frame information, typically produced by the \code{\link{make_grid_frame}} function. 
#'   Must include \code{bbox}, \code{longsize}, \code{latsize}, and \code{longcount}.
#' @param cell_id An integer representing the 0-based index of the grid cell.
#'
#' @return An `sf` \code{POLYGON} object representing the geographical boundaries of the specified grid cell.
#' @export
cell_id_to_polygon <- function (gridframe, cell_id) {
  # Calculate row and column indices
  row_index <- floor(cell_id / gridframe$longcount)
  col_index <- cell_id %% gridframe$longcount
  
  # Calculate the coordinates for the corners of the polygon
  xmin <- gridframe$bbox$xmin + col_index * gridframe$longsize
  xmax <- xmin + gridframe$longsize
  ymax <- gridframe$bbox$ymax - row_index * gridframe$latsize
  ymin <- ymax - gridframe$latsize
  
  # Create the POLYGON geometry
  polygon <- st_polygon(list(rbind(c(xmin, ymin), c(xmax, ymin), 
                                   c(xmax, ymax), c(xmin, ymax), 
                                   c(xmin, ymin))))
  return (polygon)
}

#' Calculate the centroid of a grid cell given its cell ID and a grid frame.
#'
#' @param gridframe A list containing grid frame information, typically produced by the \code{\link{make_grid_frame}} function. 
#'   Must include \code{bbox}, \code{longsize}, \code{latsize}, and \code{longcount}.
#' @param cell_id An integer representing the 0-based index of the grid cell.
#'
#' @return A numeric vector of length 2, representing the longitude (\code{x}) and latitude (\code{y}) 
#' coordinates of the centroid of the grid cell.
#'
cell_id_to_centroid <- function (gridframe, cell_id) {
  # Calculate row and column indices
  row_index <- floor(cell_id / gridframe$longcount)
  col_index <- cell_id %% gridframe$longcount
  
  # Calculate the coordinates for the corners of the polygon
  xmin <- gridframe$bbox$xmin + col_index * gridframe$longsize
  ymax <- gridframe$bbox$ymax - row_index * gridframe$latsize

  return (c(xmin + gridframe$longsize / 2, ymax - gridframe$latsize / 2))
}

#' Assign polygon geometry to a dataframe based on cell IDs
#'
#' This function converts a dataframe containing a `cell_id` column into an `sf` dataframe 
#' by assigning a polygon geometry to each row based on the corresponding grid cell.
#'
#' @param with_cell_id A dataframe that must include a column named \code{cell_id}, representing 
#'   the 0-based index of each grid cell.
#' @param gridframe A list containing grid frame information, typically produced by the 
#'   \code{\link{make_grid_frame}} function.
#'
#' @return An `sf` dataframe where each row retains the original dataframe attributes and includes 
#' a square `POLYGON` geometry representing the grid cell associated with each `cell_id`.
#'
assign_cell_geometry_sf <- function (with_cell_id, gridframe) {
  polygons <- mapply(cell_id_to_polygon, cell_id = with_cell_id$cell_id, 
                     MoreArgs = list(gridframe = gridframe), SIMPLIFY = FALSE)
  
  # Create an sf dataframe with POLYGON geometry
  sf_df <- st_sf(with_cell_id, geometry = st_sfc(polygons))
  st_crs(sf_df) <- "WGS84"
  
  return (sf_df)
}

#' Assign centroid coordinates to a dataframe
#'
#' This function adds `longitude` and `latitude` columns to a dataframe based on the centroids 
#' of grid cells corresponding to each `cell_id`.
#'
#' @param with_cell_id A dataframe that must include a column named \code{cell_id}, representing 
#'   the 0-based index of each grid cell.
#' @param gridframe A list containing grid frame information, typically produced by the 
#'   \code{\link{make_grid_frame}} function.
#'
#' @return A dataframe with the original columns, plus two new columns: \code{longitude} and 
#' \code{latitude}, representing the centroid coordinates of the grid cell for each row.
#'
assign_cell_centroids <- function (with_cell_id, gridframe) {
  # Apply calculate_centroid to each cell_id in the dataframe
  centroids <- mapply(cell_id_to_centroid, cell_id = with_cell_id$cell_id, 
                      MoreArgs = list(gridframe = gridframe))
  
  # Transpose the centroids matrix to get separate longitude and latitude vectors
  with_cell_id$longitude <- round(centroids[1, ], 6)
  with_cell_id$latitude <- round(centroids[2, ], 6)
  
  return (with_cell_id)
}

#' Assigns point sf geometry holding centroids to a dataframe
#'
#' This function adds centroid coordinates as `sf` point geometries to a dataframe based on 
#' the `cell_id` column.
#'
#' @param with_cell_id A dataframe that must include a column named \code{cell_id}, representing 
#'   the 0-based index of each grid cell.
#' @param gridframe A list containing grid frame information, typically produced by the 
#'   \code{\link{make_grid_frame}} function.
#'
#' @return An `sf` dataframe where each row includes a `POINT` geometry representing the centroid 
#' of the grid cell associated with each `cell_id`.
#'
#' @export
assign_cell_centroids_sf <- function (with_cell_id, gridframe) {
  with_coords <- assign_cell_centroids(with_cell_id, gridframe) %>% st_as_sf(coords=c("longitude", "latitude"))
  st_crs(with_coords) <- "WGS84"
  with_coords
}

# For a set of sf features which have cell_id defined, determine the cell_ids of them which intersect the supplied polygons
cells_for_polygons <- function (points, polygons) {
  intersections <- st_intersects(points, polygons, sparse = FALSE)
  # Get the indices of points that intersect with any polygon
  intersecting_points_indices <- apply(intersections, 1, any)
  intersecting_cell_ids <- points$cell_id[intersecting_points_indices]
}

# Expands a supplied region defined by a list of cell_ids by one square orthogonally in each direction
expandCells <- function (cell_ids, gridframe) {
    combined <- c(cell_ids, cell_ids + 1, cell_ids - 1, cell_ids + gridframe$longcount, cell_ids - gridframe$longcount)
    sort(unique(combined))
}