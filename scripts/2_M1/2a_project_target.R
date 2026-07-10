# ---- Libraries ----
library(sf)
library(tidyverse)
library(ggplot2)
library(terra) # pour gérer les raster
library(data.table) # pour gérer les gros df
library(rnaturalearth) #fond de carte
library(tidyterra) # ggplot terra

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Preparation----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Paths
parameters_file<-file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "../../parameters.R") 
source(parameters_file)
source ("scripts/functions.R")

# coastline to plot
coastline <- ne_coastline(scale = "medium", returnclass = "sv")%>%
  project("EPSG:3035")%>%
  st_as_sf()

# Load obs if needed
target_dt <- fread("data_prep/target/targetgroup_gbif_clean.csv")

# Load grid
grid_10000<-rast (file.path(path_raw_gridded,"grid_10000_eu.tif")) 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Rasterize----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Select obs
target_dt_10000<-target_dt[coordinateUncertaintyInMeters<=10000 & year>=1992 & year<=2020]

# Vector
pts<-vect(target_dt_10000[, .(decimalLongitude, decimalLatitude)],
     geom = c("decimalLongitude","decimalLatitude"))
length(pts)
pts<-project(pts, "EPSG:3035")

# Rasterize
target_rast_10000<-rasterize(pts, grid_10000, field = 1, fun = "sum", background = NA) %>% # just project in the grid extent even where grid=NA
  mask(grid_10000) # to have the same NA cells as the grid raster
writeRaster(target_rast_10000, file = "data_prep/target/target_rast_10000.tif", overwrite = TRUE)

# Plot
ggplot()+
  geom_spatraster(data = target_rast_10000, maxcell= 5e7) +
  scale_fill_viridis_c(name='abundance',option='magma', trans="log10", na.value = "transparent")+
  geom_sf(data = coastline, fill = NA, color = "lightgrey", size = 0.1)+
  coord_sf(xlim = c(xmin(grid_1000), xmax(grid_1000)), ylim = c(ymin(grid_1000), ymax(grid_1000)), crs = st_crs(3035))+
  theme_minimal()+
  theme(panel.grid = element_blank())+
  ggtitle("Target group abundance 10 km - 1991 - 2020")


# Hist
hist(log(target_rast_10000), main = "Log (Target group abundance per cell) - Res = 10 km"  )


