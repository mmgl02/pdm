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
grid_1000<-rast (file.path(path_raw_gridded,"grid_1000_eu.tif"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 1. Rasterize - whole period ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# to plot statistics by comparing with 10k
target_dt_1000<-target_dt[coordinateUncertaintyInMeters<=1000& year>=1992 & year<=2020 ]

# Vectorize
pts<-vect(target_dt_1000[, .(decimalLongitude, decimalLatitude)],
          geom = c("decimalLongitude","decimalLatitude"))
length(pts)
pts<-project(pts, "EPSG:3035")

# Rasterize
target_rast_1000<-rasterize(pts, grid_1000, field = 1, fun = "sum", background = NA)%>% # just project in the grid extent even where grid=NA
  mask(grid_1000) # to have the same NA cells as the grid raster 
writeRaster(target_rast_1000, file = "data_prep/target/target_rast_1000.tif", overwrite = TRUE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 2. Rasterize - by year ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# to calibrate

# Create a list with one raster for each year
years <- sort(unique(target_dt_1000$year))  

for (y in years) {
  # Select year
  dt_year <- target_dt_1000[year == y]
  
  # Vectorize
  pts_year <- vect(dt_year[, .(decimalLongitude, decimalLatitude)],geom = c("decimalLongitude", "decimalLatitude"),crs = "EPSG:4326") %>%
    project("EPSG:3035")
  
  # Rasterize, mask and save
  rasterize(pts_year, grid_1000, field = 1, fun = "sum", background = NA) %>%
    mask(grid_1000) %>%
    writeRaster(paste0("data_prep/target/target_rast_1000_", y, ".tif"), overwrite = TRUE)
  
  cat("Saved year:", y, "\n")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 3. Plot ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
target_rast_1000_2020<-rast("data_prep/target/target_rast_1000_2020.tif")
plot(target_rast_1000_2020)

target_rast_1000_list<-lapply(years, function(y) rast(paste0("data_prep/target/target_rast_1000_", y, ".tif")))
names(target_rast_1000_list)<-years
plot(target_rast_1000_list[[10]])

# Plot
plot(log10(target_rast_1000), main= "Target group abundance 1 km - 1992 - 2020 (log10)")

# Number of observations per year
obs_per_year <- target_dt_1000[, .N, by = year][order(year)] # .N=number
ggplot(obs_per_year, aes(x = year, y = N)) +
  geom_point() +
  labs(title = "Number of target goup observations per year (Uncertainty < 1000 m)",
       x = "Year",
       y = "Number of observations")+
  theme_minimal()

# Number of grid observations per year
grid_obs_per_year <- data.table(
  year = names(target_rast_1000),
  n_cells = sapply(target_rast_1000_list, function(r) sum(values(r) > 0, na.rm = TRUE))
)
grid_obs_per_year$year<-as.numeric(grid_obs_per_year$year)

save(grid_obs_per_year, file = "results/2_targetgroup/grid_obs_per_year_1000.RData")

ggplot(grid_obs_per_year, aes(x = year, y = n_cells)) +
  geom_point()+
  labs(title = "Number of target group grid observations per year (Resolution= 1000 m)",
       x = "Year",
       y = "Number of observations")+
  theme_minimal()

