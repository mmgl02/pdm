library(sf)
library(tidyverse)
library(ggplot2)
library(terra) # pour gérer les raster
library(data.table) # pour gérer les gros df
library(patchwork) #plot
library(readxl)
library(corrplot)

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
grid_10000<-rast (file.path(path_raw_gridded,"grid_10000_eu.tif")) 
# already mask and cropped in the eurostat europe study area. Otherwise need to crop and mask
# env variables again after having projected into the grid.

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Bioclim ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load covariables and do the mean over the years
years <- 1992:2020

vars <- c("gdd5","bio01d", "bio04d", "bio12d", "bio15d")
#vars<-c("bio01d")

process_bioclim <- function(var_name, years, input_dir, mask_pre, grid) {
  if (var_name=="bio15d") { years<-1992:2019} # because 2020 is deprecated
  
  # path
  files <- paste0(input_dir,"/CHELSA_EUR11_obs_",var_name,"_",years,"_V.2.1.nc")
  
  # Load: 1 multilayer raster with 1 layer= 1 year
  r <- rast(files)
  
  # Mean over the years 
  r_mean <- app(r, mean)
  
  # Crop and mask to Europe in 4326 before projecting (to avoid projecting empty cells)
  r_mean <- crop(r_mean, mask_pre)%>% 
    mask(mask_pre)
  
  # Project # crop again if the grid was not projected otherwise the extent of the raster is too big because of the projection # mask again to have the same NA cells as the other rasters
  r_mean <- project(r_mean, grid, method = "bilinear") #%>%
    #crop(mask_post) %>%
    #mask(mask_post)

  print(paste0(var_name, " projected"))
  return(r_mean)
}

results <- lapply(vars, function(v) { # 5 min
  process_bioclim(
    var_name = v,
    years = years,
    input_dir = path_raw_bioclim,
    mask_pre=europe_4326_buf,
    grid=grid_10000
  )
})

# Stack
r_stack <- rast(results)
names(r_stack) <- vars

# Plot to verify
plot(r_stack)

# Save
writeRaster(r_stack,filename = file.path(path_prep_env, "M1", "bioclim_mean_1992_2020.tif"),overwrite = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Remarques  ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# si je plot les rasters après avoir fait rast(files), les titres sont des jours qui correspondent à la moitié de l'année, mais c'est bien la moyenne qui est plotée
file<- paste0 (path_raw_bioclim,"/CHELSA_EUR11_obs_bio01d_1992_V.2.1.nc")
r <- rast(file)
r<-crop(r, europe_4326_buf) %>% mask(europe_4326_buf)
writeRaster(r, file.path(path_prep_env, "M1", "r_4326.tif"), overwrite = TRUE)

r_proj2<-project(r, grid_10000, method = "bilinear")
writeRaster(r_proj2, file.path(path_prep_env, "M1", "r_proj2.tif"), overwrite = TRUE)


g<-extend(grid,5)
writeRaster(g, file.path(path_prep_env, "M1", "g.tif"), overwrite = TRUE)

eu<-EUROSTAT_Countries_1k<- st_read(file.path(path_raw_mask_europe,"EUROSTAT_Countries_buffered1k_WGS84_nng (2024_07_11 13_01_52 UTC).shp"))
