# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# From observations to gridded presences
# Input :
# obs_dt: observations , output of clean coordinates
# grid: grille chargée avec rast() au crs EPSG:3035
# resolution: resolution de la grille en mètres
# year: True if we want to aggregate by cell and year, false if we want to aggregate by cell only
# Output:
# data.table with columns: cell_id, species, lon, lat, (year if year=T)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

obs_to_grid <- function (obs_dt, grid, resolution, year=F){
  
  # Filter observations at precision uncertainty
  obs_dt<-obs_dt[coordinateUncertaintyInMeters <= resolution]
  
  # Vectorize
  obs_vec <- vect(
    obs_dt,
    geom = c("decimalLongitude", "decimalLatitude"),
    crs = "EPSG:4326"
  )

  # Reproject into the CRS of the grid
  obs_vec <- project(obs_vec, "EPSG:3035")
  
  # Associate observations to grid cells
  #get cell number in the grid for each observation based on its coordinates
  obs_dt[, cell_id := cellFromXY(grid, crds(obs_vec))] # create a new column cell_id in the observation data table
  
  # Columns to group by
  by_cols <- c("cell_id", "species")
  if (year) by_cols <- c(by_cols, "year")
  
  # Get unique observations by columns to group by
  obs_unique <- unique(obs_dt, by = by_cols)
  
  # Species name with _ and not a space
  obs_unique[, species := gsub("\\s+", "_", species)]

  # Add centroid coordinates
  coords <- xyFromCell(grid, obs_unique$cell_id)
  obs_unique[, `:=`(lon = coords[, 1], lat = coords[, 2])]
  
  # Columns to return
  cols <- c("cell_id", "species", "lon", "lat")
  if (year) cols <- c(cols, "year")
  return(obs_unique[, ..cols])
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Selection of background points 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

select_background_target <- function(presence_coords, nb.background, weighted_raster){
  
  # Find cell IDs at presences
  cell_IDs <- cellFromXY(object=weighted_raster, xy=presence_coords)
  
  # Set presence cells to NA to prevent selection
  weighted_raster[cell_IDs] <- NA
  cat("cell_id")
  
  # Random selection
  background_points <- spatSample(weighted_raster, size=nb.background, method = "weights", replace  = FALSE, na.rm = TRUE,
                                  xy = TRUE)
  cat("selected")
  return(background_points)
}
