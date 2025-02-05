library(sf)

lat_lon <- function (data) {
  return (st_transform(data, "+proj=longlat +datum=WGS84"))
}

roundmulti <- function (multi, digits) {
  multi <- lapply(multi, function (matrix) {
    matrix <- lapply(matrix, function (coords) {
      round(coords, digits)
    })
  })
  return (st_multipolygon(multi))
}

roundpoly <- function (poly, digits) {
  poly <- lapply(poly, function (matrix) {
      round(matrix, digits)
  })
  return (st_polygon(poly))
}

round_sf <- function (fc, digits) {
  # https://gis.stackexchange.com/questions/329110/removing-empty-polygon-from-sf-object-in-r
  simple  <- fc %>% st_simplify(preserveTopology = TRUE, dTolerance = 5) %>% dplyr::filter(!st_is_empty(.))
  geom <- simple$geometry
  geom <- lapply(geom, function (one) {
    if (inherits(one, "MULTIPOLYGON")) {
      one <- roundmulti(one, digits)
    } else if (inherits(one, "POLYGON")) {
      one <- roundpoly(one, digits)
    } else if (inherits(one, "XY")) {
      one <- round(one)
    } else if (!st_is_empty(one)) {
      stop(paste("I don't know what it is ", class(one)))
    }
  })
  simple$geometry <- st_sfc(geom)
  simple
}

mx_read <- function (filename, digits = 4) {
  st_data <- st_read(filename, quiet=TRUE);
  dropped <- st_zm(st_data, drop = T, what = "ZM")
  trans <- lat_lon(dropped);
  rounded <- round_sf(trans, digits);
}

# Read a CSV file fast with the proper encoding - note that read.csv is generally faulty
timedFread <- function (toread) {
  start <- Sys.time()
  frame <- data.table::fread(toread, encoding = "UTF-8")
  end <- Sys.time()
  message("Read ", nrow(frame), " rows from ", toread, " in ", (end - start), "s")
  # Otherwise traditional R indexing notation fails
  as.data.frame(frame)
}

timedWrite <- function (x, towrite) {
  start <- Sys.time()
  # Approach for selective quoting taken from https://stackoverflow.com/a/25810538/1381443
  commas <- which(sapply(x, function(y) any(grepl(",",y))))
  write.csv(x, towrite, na = "", row.names = FALSE, quote = commas, fileEncoding = "UTF-8")
  end <- Sys.time()
  message("Written ", nrow(x), " rows to ", towrite, " in ", (end - start), "s")
}

# Attach the region's label as an "mx_regionId" option in the output data
labelToOption <- function (label) {
  return (list(mx_regionId = label))
}

# Cribbed from https://kevinushey.github.io/blog/2018/02/21/string-encoding-and-r/
write_utf8 <- function(text, f) {
  # step 1: ensure our text is utf8 encoded
  utf8 <- enc2utf8(text)
  
  # step 2: create a connection with 'native' encoding
  # this signals to R that translation before writing
  # to the connection should be skipped
  con <- file(f, open = "w+", encoding = "native.enc")
  
  # step 3: write to the connection with 'useBytes = TRUE',
  # telling R to skip translation to the native encoding
  writeLines(utf8, con = con, useBytes = TRUE)
  
  # close our connection
  close(con)
}