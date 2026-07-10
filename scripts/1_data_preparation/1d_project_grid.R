library(terra)
library(dplyr)

# Paths
parameters_file<-file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "../../parameters.R") 
source(parameters_file)

# Load grids
grid_10000<- rast (file.path(path_raw_gridded,"ReferenceGrid_Europe_bin_10000m.tif"))
grid_1000<- rast (file.path(path_raw_gridded,"ReferenceGrid_Europe_bin_1000m.tif"))
grid_100<- rast (file.path(path_raw_gridded,"ReferenceGrid_Europe_100m.tif"))
grid_50000<- rast("data_raw/Gridded/ReferenceGrid_Europe_bin_50000m.tif")

# Load mask europe
europe_3035<-vect(file.path(path_raw_mask_europe,"europe_eurostat.shp")) %>%
  aggregate() %>%
    project("EPSG:3035")

plot(europe_3035)
crs(europe_3035)

# Mask in Europe
grid_10000_eu<- crop(grid_10000, europe_3035) %>% mask(europe_3035)
plot(grid_10000_eu)
writeRaster(grid_10000_eu, file.path(path_raw_gridded,"grid_10000_eu.tif"), overwrite=TRUE)

grid_1000_eu<- crop(grid_1000, europe_3035) %>% mask(europe_3035)
plot(grid_1000_eu)
writeRaster(grid_1000_eu, file.path(path_raw_gridded,"grid_1000_eu.tif"), overwrite=TRUE)

grid_50000_eu<- crop(grid_50000, europe_3035) %>% mask(europe_3035)
plot(grid_50000_eu)
writeRaster(grid_50000_eu, file.path(path_raw_gridded,"grid_50000_eu.tif"), overwrite=TRUE)

# Load mask France
france_3035<-vect(file.path(path_raw_mask_france, "ARRONDISSEMENT.shp"))%>% # from  https://geoservices.ign.fr/geofla # EPSG 2154 Lambert 93
  aggregate () %>%
  project("EPSG: 3035")

# Mask in France
grid_10000_fr<- crop(grid_10000, france_3035) %>% mask(france_3035)
writeRaster(grid_10000_fr, file.path(path_raw_gridded,"grid_10000_fr.tif"), overwrite=TRUE)

grid_1000_fr<- crop(grid_1000, france_3035) %>% mask(france_3035)
plot(grid_1000_fr)
writeRaster(grid_1000_fr, file.path(path_raw_gridded,"grid_1000_fr.tif"), overwrite=TRUE)

grid_100_fr<- crop(grid_100, france_3035) %>% mask(france_3035)
writeRaster (grid_100_fr, file.path(path_raw_gridded,"grid_100_fr.tif"), overwrite=TRUE)

# Mask around french alps
alp_3035<-vect(file.path(path_raw_mask_france, "alps_fr.shp"))%>% project("EPSG: 3035")

grid_100_alp<- crop(grid_100, alp_3035)
plot(grid_100_alp)
writeRaster (grid_100_alp, file.path(path_raw_gridded,"grid_100_alp.tif"), overwrite=TRUE)

# Mask in jura
jura_3035<-vect(file.path(path_raw_mask_france, "jura.shp"))%>% project("EPSG: 3035")
grid_100_jura<- crop(grid_100, jura_3035) %>% mask(jura_3035)
plot(grid_100_jura)
writeRaster (grid_100_jura, file.path(path_raw_gridded,"grid_100_jura.tif"), overwrite=TRUE)

# Mask around South alps
alp_south_3035<-vect(file.path(path_raw_mask_france, "alp_sud.shp"))%>% project("EPSG: 3035")
grid_100_alp_south<- crop(grid_100, alp_south_3035) %>% mask(alp_south_3035)
plot(grid_100_alp_south)
writeRaster (grid_100_alp_south, file.path(path_raw_gridded,"grid_100_alp_s.tif"), overwrite=TRUE)
