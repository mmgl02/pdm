library(tidyverse)
library(terra) # pour gérer les raster
library(data.table) # pour gérer les gros df
library(readxl)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Preparation ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Paths
parameters_file<-file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "../../parameters.R") 
source(parameters_file)
source ("scripts/functions.R")

## for bioclim, landcover which are in epsg 4326
europe_4326<-vect(file.path(path_raw_mask_europe,"europe_eurostat.shp")) %>%
  makeValid()%>%
  aggregate() # aggregate the vector

# buffer, otherwise when projecting the env raster there are empty cells at the coast
europe_4326_buf <- europe_4326 %>%
  project("EPSG:3035") %>%          
  buffer(width = 20000) %>%     
  project("EPSG:4326") %>%          
  makeValid()  

# Grid
grid_1000<-rast (file.path(path_raw_gridded,"grid_1000_eu.tif")) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 1. Bioclim by year ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
years <- 1992:2020
vars<- c("bio01d","bio04d", "bio12d", "bio15d", "gdd5") #"scd" 

# Function to load, crop, project raster for one year
process_bioclim_year <- function(var_name, year, input_dir, mask_pre, grid) {
  
  if (var_name == "bio15d" & year == 2020) return(NULL) # because 2020 is deprecated
  
  file <- paste0(input_dir, "/CHELSA_EUR11_obs_", var_name, "_", year, "_V.2.1.nc")
  cat("Processing:", var_name, "year:", year, "\n")
  # Load
  rast(file) %>%
    crop(mask_pre) %>%
    mask(mask_pre) %>%
    project(grid, method = "bilinear") %>%
    setNames(var_name)
}

# Apply to all years # 1h
for (yr in years) {
  layers <- lapply(vars, function(v) process_bioclim_year(v, yr, path_raw_bioclim, europe_4326_buf, grid_1000))
  r_stack <- rast(Filter(Negate(is.null), layers)) # stack and suppress NULL layers (bio15d for 2020)
  writeRaster(r_stack,
              filename = file.path(path_prep_env, "M2", paste0("bioclim_", yr, ".tif")),
              overwrite = TRUE)
  cat("Saved year:", yr, "| layers:", paste(names(r_stack), collapse = ", "), "\n")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 2. LC by year ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
start<- Sys.time()

# Parameters
years <- 1992:2020
LC_files <- c(
  paste0("data_raw/landcover/byYear/ESACCI-LC-L4-LCCS-Map-300m-P1Y-", 1992:2015, "-v2.0.7.tif"),
  paste0("data_raw/landcover/byYear/C3S-LC-L4-LCCS-Map-300m-P1Y-",   2016:2020, "-v2.1.1.tif")
)

correspondances <- read_excel("data_raw/landcover/legend/broad_categories.xlsx")
lc_classnames <- unique(correspondances[1:35, "class"])$class

# function to binarize a pixel x of a raster (1 if the cell belongs to the class, 0 otherwise)
# input x =value of a pixel (1,2,3,4,5,6) = class of landcover
# output = vector of 6 values. ex: (1,0,0,0,0,0) if the pixel belongs to class 1
binarize_lc <- function(x){
  sapply(1:6, function(class) as.integer(x == class)) # test if pixel belongs to class 1, 2 ...
} 

# Function to load, crop, reclassify, binarize and project raster for one year 
process_lc_year <- function(file, mask_pre, grid) {
  
  rast(file) %>% # load
    crop(mask_pre) %>% # crop to save time
    mask(mask_pre) %>% # mask to save time
    classify(correspondances[, c("is", "becomes")]) %>% # transform to 6 classes (1,2,3,4,5,6)
    app(binarize_lc) %>% # binarize to have 6 layers with 1 if the cell belongs to the class, 0 otherwise
    project(grid, method = "mean") %>% # give % of each class in each cell
    setNames(lc_classnames)
}

# For each year
for (i in seq_along(years)) {
  r <- process_lc_year(LC_files[i], europe_4326_buf, grid_1000)
  writeRaster(r,
              filename = file.path(path_prep_env, "M2", paste0("LC_", years[i], ".tif")),
              overwrite = TRUE)
  cat("Saved year:", years[i], "| layers:", paste(names(r), collapse = ", "), "\n")
}

print(start - Sys.time())

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 3. Bioclim mean 1992-2020 if needed----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load all annual rasters
raster_files <- file.path(path_prep_env, "M2", paste0("bioclim_", years, ".tif"))

rasters_annual_list <- lapply(raster_files, rast)

# Compute mean across years for each variable and stack into a single raster bioclim_mean
bioclim_mean <- rast(lapply(vars, function(v) {
  cat ("Processing variable:", v, "\n")
  if (v == "bio15d") {rasters_annual_list <- rasters_annual_list[1:(length(rasters_annual_list)-1)]} # because 2020 is deprecated}
  
  # Extract layer of the variable v from each annual raster
  # --> Get one raster per variable with all years as layers
  rast_v_all_years <- rast(lapply(rasters_annual_list, function(r) {
    cat("year", names(r), "\n" )
    r[[v]]
  }))
  
  # Mean across years
  app( rast_v_all_years, mean, na.rm = TRUE) %>%
    setNames(v) # name of the variable v
}))

writeRaster(bioclim_mean,
            filename = file.path(path_prep_env, "M2", "bioclim_mean.tif"),
            overwrite = TRUE)

