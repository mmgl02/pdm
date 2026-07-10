# ---- Libraries ---- 
library(sf)
library(tidyverse)
library(ggplot2)
library(terra) # pour gérer les raster
library(data.table) # pour gérer les gros df
library(patchwork) #plot
library(rnaturalearth) #fond de carte
library(tidyterra) # ggplot terra
library(biomod2)
library(blockCV)

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

# Load grid
grid<-rast (file.path(path_raw_gridded,"grid_1000_eu.tif")) 
grid_50000<- rast("data_raw/Gridded/ReferenceGrid_Europe_bin_50000m.tif") # to plot

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Load observations----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# get species files
species_files <- list(
  "Bombus hortorum"   = file.path(path_prep_obs, "data_gbif_clean.csv") ,
  "Bombus pascuorum"  = file.path(path_prep_obs, "data_gbif_clean.csv"),
  "Bombus terrestris" = file.path(path_prep_obs, "data_gbif_clean.csv"),
  "Bombus alpinus"    = file.path(path_prep_obs, "balpins_data_gbif_clean.csv"),
  "Bombus pyrenaeus"  = file.path(path_prep_obs, "balpins_data_gbif_clean.csv"),
  "Bombus wurflenii"  = file.path(path_prep_obs, "balpins_data_gbif_clean.csv") 
)
# get each filename once
unique_files <- unique(unlist(species_files))

# Load each obs file and filter in 1992-2020
obs_list <- lapply(unique_files, function(f) {
  fread(f)[year >= 1992 & year <= 2020]
})
names(obs_list) <- unique_files # associate names

# Combine observations for all species into one datatable
obs_dt_all <- rbindlist(lapply(names(species_files), function(sp) {
  obs_list[[species_files[[sp]]]][species == sp] # get species file, and filter it by this species to keep only its observations
}))

# Number of observations per species 
obs_dt_all %>%
  group_by(species) %>%
  summarise(count = n()) %>%
  print()

# Get species names in the format "Bombus_hortorum" for later
species_names <- gsub("\\s+", "_", names(species_files))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Load target group ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# to select background points
years<- 1992: 2020

target_list<-lapply(years, function(y){
  rast(paste0("data_prep/target/target_rast_1000_", y,".tif"))
} )%>%
  setNames(years)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Load environment ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Bioclim by year
bioclim_list <- lapply(years, function(y){
  rast(paste0("data_prep/env/M2/bioclim_", y,".tif"))
} ) %>%
  setNames(years)

# LC by year
LC_list <- lapply(years, function(y){
  rast(paste0("data_prep/env/M2/LC_", y,".tif"))
} )%>%
  setNames(years)

# Variables to calibrate the model with, selected in c_choose_env
vars<- c("gdd5", "bio04d", "bio12d", "bio15d", "Agriculture", "Forest", "Grassland", "Wetland", "Settlement", "Others" )
vars_bioclim<- c("gdd5", "bio04d", "bio12d", "bio15d")
vars_LC<-c("Agriculture", "Forest", "Grassland", "Wetland", "Settlement", "Others" )

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 1. Project species data  ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# all species at the same time

# Associate observations to the grid at resolution 1 km in Europe, get a datatable
dt <- obs_to_grid (obs_dt_all, grid, 1000, year=T) # from functions.R

# Save
save (dt, file = "data_prep/obs/M2/dt_1000.RData")

# Suppress useless files
rm(obs_list, unique_files, obs_dt_all)
gc()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 2. Select background  ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# all species at the same time

# Annual target observation effort #10s
n_target_per_year<-rast (target_list) %>%
  global( sum, na.rm=TRUE)%>% #long
  setNames("years")

# Build weighted raster for each year
weighted_raster_list <- lapply(years, function(y) {
  r<-cover(target_list[[as.character(y)]], grid * 0) %>% # 0 in the continents, NA in the seas
    mask(bioclim_list[[1]][[1]])%>% # do not select where there is no bioclim value
    mask(LC_list[[1]][[1]]) # do not select where there is no LC value
  r<-r+1  
})%>%
  setNames(years)

# to run next time instead of what's above:
weighted_raster_list <- lapply(years, function(y) {
  r<-cover(target_list[[as.character(y)]], grid * 0) %>% # 0 in the continents, NA in the seas
    mask(bioclim_list[[as.character(y)]][[1]])%>% # do not select where there is no bioclim value
    mask(LC_list[[as.character(y)]][[1]]) # do not select where there is no LC value
  r<-r+1  
})%>%
  setNames(years)

# Function to select background points for one species for each year
selecting_background <- function(sp) {
  cat("Selecting background for species:", sp, "\n")
  
  # Select dt for the species
  dt_sp <- dt[species == sp, .(lon, lat, year)]
  nb.background <- ifelse(nrow(dt_sp) < 10000, 10000, nrow(dt_sp))
  
  # Sample years proportionally to annual observation effort
  y_sampled <- sample(rownames(n_target_per_year), nb.background, replace = TRUE,
                      prob = n_target_per_year$sum)
  y_table <- table(y_sampled)
  print(y_table)
  
  # Select background for each year
  background_points <- rbindlist(lapply(names(y_table), function(y) { 
    
    # number of background proportional to annual obs effort
    nb.background_y<- y_table[y]
    print(paste0(nb.background_y, " background points will be sampled for year ", y))
    
    # Presence coordinates
    presence_coords <- dt_sp[year == as.numeric(y), .(lon, lat)]
    
    pts <- select_background_target(presence_coords, nb.background_y, weighted_raster_list[[y]]) %>% #from functions.R
      as.data.table()
    pts$year <- as.numeric(y)
    pts
  }))
  
  # Proportion of background points outside the target group observation effort
  proportion <- background_points[sum == 1, .N] / nrow(background_points)
  cat("Species:", sp, "| n presences:", nrow(dt_sp), "| n background:", nb.background,
      "| prop outside target:", round(proportion, 3), "\n")
  
  # Bind background with presences
  rbind(
    dt_sp[, .(lon, lat, year, presence = 1L, species = sp)],
    background_points[, .(lon = x, lat = y, year, presence = 0L, species = sp)]
  )
}

# Apply to all species # 10 min for 3 sp
set.seed(123)
dt_background <- rbindlist(lapply(species_names, selecting_background))
save (dt_background, file = "data_prep/obs/M2/dt_1000_background_ori.RData")

rm(target_list, n_target_per_year, weighted_raster_list)
gc()
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 3. Associate to env  ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Bioclim
# extract values from the rasters
for (y in names(bioclim_list)) {

  for (var in vars_bioclim) {
    # Skip bio15d for 2020
    if (y == "2020" & var == "bio15d") next
    
    r <- bioclim_list[[y]][[var]]
    
    # extract values from the raster
    vals <- terra::extract(r, vect(dt_background[year==as.numeric(y), .(lon, lat)], geom = c("lon","lat"), crs = "EPSG:3035"))
    
    # add values to dt
    dt_background[year==as.numeric(y), (var) := vals[,2]]
    
    print(paste0("values added for ", var, y))
  }
}

# for 2020 and bio15 we will affect the values of 2019
vals <- terra::extract(bioclim_list$`2019`$bio15d, vect(dt_background[year==2020, .(lon, lat)], geom = c("lon","lat"), crs = "EPSG:3035"))
dt_background[year==2020, ('bio15d') := vals[,2]]

# LC # do the same
for (y in names(LC_list)) {
  for (var in vars_LC){
    r <- LC_list[[y]][[var]]
    vals <- terra::extract(r, vect(dt_background[year==as.numeric(y), .(lon, lat)], geom = c("lon","lat"), crs = "EPSG:3035"))
    dt_background[year==as.numeric(y), (var) := vals[,2]]
    print(paste0("values added for ", var, y))
  }
}

# Remove NA (FOR THE MOMENT, sinon on doit faire ça avant de selectionner les pseudo-absences, avec un mask de env) OK
dt_background <- dt_background[complete.cases(dt_background)]

dt_background%>%
  group_by(species, presence) %>%
  summarise(count = n()) %>%
  print()

# Save
save (dt_background, file = "data_prep/obs/M2/dt_1000_background.RData")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 4. Block CV for evaluation (function) ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#load("data_prep/obs/M2/dt_1000_background.RData")
sf_background<- st_as_sf(dt_background, coords = c("lon", "lat"), crs = "EPSG:3035") # for blockCV, we need an sf objectt

# Spatial blocking function
spatial_blocking <- function (sp){
  cat("spatial blocking\n")
  
  # Filter for one species
  sf_sp <- sf_background[sf_background$species == sp, ]
  
  # Load existing blocks
  load(paste0("data_prep/blocks/blocks_", sp, ".RData"))
  
  # Extract fold ID for each point
  folds <- st_join(sf_sp, scv$blocks["folds"])$folds
  
  # Convert to biomod2 user.table format: boolean matrix (n_points x n_folds)
  n_folds <- max(folds, na.rm = TRUE) #5
  user_table <- sapply(1:n_folds, function(k) folds != k)  # TRUE = in calibration, FALSE = in validation
  colnames(user_table) <- paste0("_allData_RUN", 1:n_folds)
  
  return(user_table)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 5. Format data (function) ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# build species modelling wrapper
formating <- function(sp){
  cat("formating")
  
  # get presence absence points
  dt_sp <- dt_background[dt_background$species == sp, ]
  
  # format data
  myBiomodData <- BIOMOD_FormatingData(resp.name = sp,
                                       resp.var = dt_sp[, presence],
                                       resp.xy = dt_sp [,.(lon,lat) ],
                                       expl.var = dt_sp [,..vars]
  )
  
  # Print summary
  print(myBiomodData)
  
  return (myBiomodData)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 6. Calibrate (function) ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

calibrating <-function(sp, myBiomodData, spatial_cv_folds){
  start <- Sys.time()
  cat("calibrating")
  
  # Define options
  myOptions<- bm_ModelingOptions(data.type = 'binary', models='RF', strategy = 'bigboss')
  print(myOptions)
  
  myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                      models="RF",
                                      modeling.id = "M2_1",
                                      CV.strategy = "user.defined",
                                      CV.user.table   = spatial_cv_folds,
                                      OPT.user = myOptions,
                                      metric.eval = c('TSS', 'AUCroc', "BOYCE"),
                                      var.import = 2,
                                      CV.do.full.models=FALSE,
                                      do.progress=TRUE)
  
  cat("Computation time :", Sys.time() - start, "\n")
  return (myBiomodModelOut)
  
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 7. Evaluate (function) ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
evaluating <- function(sp, myBiomodModelOut) {
  cat("evaluating:", sp, "\n")
  
  # Create plots directory
  plot_dir <- file.path(myBiomodModelOut@dir.name, myBiomodModelOut@sp.name, "models", myBiomodModelOut@modeling.id, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Evaluation scores
  cat ("Evaluation scores:\n")
  eval <- as.data.frame(get_evaluations(myBiomodModelOut))
  eval_summary <- eval %>% group_by(metric.eval) %>% summarise(mean = mean(validation))
  write.csv(eval_summary, file.path(plot_dir, "eval_summary.csv"), row.names = FALSE)
  
  # Eval boxplot
  cat ("Evaluation boxplot:\n")
  png(file.path(plot_dir, "eval_boxplot.png"))
  bm_PlotEvalBoxplot(bm.out = myBiomodModelOut, dataset = "validation", group.by = c("algo", "run"))
  dev.off()
  
  # Variable importance
  cat ("Variable importance:\n")
  png(file.path(plot_dir, "varimp_boxplot.png"))
  bm_PlotVarImpBoxplot(bm.out = myBiomodModelOut, group.by = c("expl.var", "algo", "algo"))
  dev.off()
  
  # Response curves
  cat ("Response curves:\n")
  mods <- get_built_models(myBiomodModelOut, run = "RUN1")
  png(file.path(plot_dir, "response_curves.png"))
  bm_PlotResponseCurves(bm.out = myBiomodModelOut, models.chosen = mods, fixed.var = "median")
  dev.off()
  
  # Observation and background points
  cat ("myBiomodData:\n")
  png(file.path(plot_dir, "myBiomodData.png"))
  plot(myBiomodData)
  dev.off()
  
  # Data of the species to plot
  dt_sp <- dt_background[dt_background$species == sp, ]
  
  # Plot raster data
  r<-rasterize(dt_sp, grid, values=dt_sp$presence)
  png(file.path(plot_dir, "myBiomodData_raster.png"))
  plot(r, col=c("grey","darkgreen"), main= paste("Presences and background for", sp))
  dev.off()
  
  # Plot raster data in a 50 k grid
  r<-rasterize(dt_sp, grid_50000, values=dt_sp$presence, fun='max')
  png(file.path(plot_dir, "obs_raster.png"))
  plot(r, col=c("grey","darkgreen"), main= paste("Presences and background rasterized in 50 km cells (max) \n for", sp))
  dev.off()
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----8. Ensemble modelling (function) ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# weighted by AUC

ensemble_modelling <- function(sp, myBiomodModelOut) {
  cat("Ensemble modelling:", sp, "\n")

  myBiomodEM <- BIOMOD_EnsembleModeling(
    bm.mod = myBiomodModelOut,
    models.chosen = "all",
    em.by = "all",
    em.algo = 'EMwmean',
    metric.select = "AUCroc",
    metric.eval = c("AUCroc", "TSS"),
    var.import = 0,
    EMwmean.decay = "proportional",
    do.progress = TRUE
  )
  
  return(myBiomodEM)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 9. Project (function) different from 2d in Europe ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# env raster # same for all species # 25 min for 1 projection in Europe
year<-"2019"
env <- c(bioclim_list[[year]][[vars_bioclim]], LC_list[[year]][[vars_LC]])

# function to project
projecting_eu <- function(sp, year, new.env) {
  cat("Projecting:", sp, "\n")
  
  myBiomodModel <- get(load( file.path(path_results_models, gsub("_", ".", sp), paste0(gsub("_", ".", sp), ".M2_1.models.out") )))
  myBiomodEM <- get(load( file.path(path_results_models, gsub("_", ".", sp), paste0(gsub("_", ".", sp), ".M2_1.ensemble.models.out") )))
  
  myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModel,
                                             proj.name = paste0(myBiomodModel@modeling.id, "_", year, "_Eu"),
                                             new.env = new.env, # only the variables used for calibration
                                             models.chosen = 'all',
                                             build.clamping.mask = FALSE,
                                            do.stack= FALSE,
                                            keep.in.memory = FALSE,
                                            overwrite=TRUE)
  myBiomodProj_EM<-BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, bm.proj=myBiomodProj)
  
  
  
  # Load raster and save it
  plot_dir <- file.path(myBiomodModel@dir.name, myBiomodModel@sp.name, "models", myBiomodModel@modeling.id, "plots")# same as in evaluating

  r <- rast(file.path(myBiomodModel@dir.name, myBiomodModel@sp.name, paste0("proj_", myBiomodProj@proj.name), "individual_projections", 
                      paste0(myBiomodModel@sp.name, "_EMwmeanByAUCroc_mergedData_mergedRun_mergedAlgo.tif")))
  pdf(file.path(plot_dir, "projected_distribution_eu.pdf"))
  plot(r, main = paste("HS2", year, sp))
  dev.off()
  
  return(myBiomodProj)
}

species_names<-  c("Bombus_alpinus",
                   "Bombus_pyrenaeus",
                   "Bombus_wurflenii",
                   "Bombus_pascuorum",
                   "Bombus_terrestris",
                   "Bombus_hortorum")

for (sp in species_names) {
  projecting_eu(sp, year, env)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Apply 4, 5 and 6, 7, 8, 9  to each species ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# setwd(root)
gc()
species_names<-  c("Bombus_alpinus",
                   "Bombus_pyrenaeus",
                   "Bombus_wurflenii")

for (sp in species_names) {
  
  set.seed(123)
  spatial_cv_folds<-spatial_blocking(sp)
  
  setwd(path_results_models)
  myBiomodData<-formating(sp)
  myBiomodModelOut<-calibrating(sp, myBiomodData, spatial_cv_folds)
  evaluating(sp, myBiomodModelOut)
  myBiomodEM<-ensemble_modelling(sp, myBiomodModelOut)
  projecting(sp, myBiomodEM, "2019", env)
  setwd(root)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 9. Project (function) in France different from 2d  ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# env raster # same for all species # 25 min for 1 projection in Europe
year<-"2019"
env <- c(bioclim_list[[year]][[vars_bioclim]], LC_list[[year]][[vars_LC]])

# better to project only in France, so crop env in France # projecting: 2 min
## France
france_3035<-vect(file.path(path_raw_mask_france, "ARRONDISSEMENT.shp"))%>% # from  https://geoservices.ign.fr/geofla # EPSG 2154 Lambert 93
  aggregate () %>%
  project("EPSG: 3035")

## Crop and mask env
env<-crop (env, france_3035) %>%
  mask(france_3035)

# function to project
projecting_fr <- function(sp, myBiomodEM, year, new.env) {
  cat("Projecting:", sp, "\n")
  
  myBiomodProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                             proj.name = paste0(myBiomodEM@modeling.id, "_", year, "_France"),
                                             new.env = new.env, # only the variables used for calibration
                                             models.chosen = 'all',
                                             build.clamping.mask = TRUE)
  
  # Load raster and save it
  plot_dir <- file.path(myBiomodEM@dir.name, myBiomodEM@sp.name, "models", myBiomodEM@modeling.id, "plots")# same as in evaluating
  
  r <- rast(file.path(myBiomodEM@dir.name, myBiomodEM@sp.name, paste0("proj_", myBiomodProj@proj.name), paste0("proj_", myBiomodProj@proj.name,"_", myBiomodEM@sp.name, "_ensemble.tif")))
  # %>% mask(grid) #unecessary to mask for france
  
  pdf(file.path(plot_dir, "projected_distribution.pdf"))
  plot(r, main = paste("HS2", year, sp))
  dev.off()
  
  return(myBiomodProj)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----Load results ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load presence background
dt_background<-get(load("data_prep/obs/M2/dt_1000_background.RData"))

# load my biomodmodel out
myBiomodModelOut<-get(load( file.path(path_results_models, "Bombus.alpinus/Bombus.alpinus.M2_1.models.out") ))
myBiomodModelOut<-get(load( file.path(path_results_models, "Bombus.hortorum/Bombus.hortorum.M2_1.models.out") ))

#load mybiomodEM
myBiomodEM<-get(load( file.path(path_results_models, "Bombus.alpinus/Bombus.alpinus.M2_1.ensemble.models.out") ))
myBiomodEM<-get(load( file.path(path_results_models, "Bombus.hortorum/Bombus.hortorum.M2_1.ensemble.models.out") ))

# sp
sp<-"Bombus_alpinus"
sp<-"Bombus_hortorum"

#blocks
load(file.path("data_prep/blocks/blocks_Bombus_alpinus.RData"))
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----TESTS  ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

r<-rast(file.path(path_results_models, gsub("_", ".", sp), "proj_Current", paste0("proj_Current_", gsub("_", ".", sp), "_ensemble.tif")))%>%
  mask(grid)


plot (r, main = paste("HS1 ", sp))

na_points <- dt_background[is.na(bio15d)]
plot(grid)
points(na_points[, .(lon, lat)], 
       col = ifelse(na_points$presence == 1, "red", "blue"), 
       pch = 16)



print(blocks)
plot(blocks)

# plot dt sp
r<-rasterize(dt_sp, grid, values=dt_sp$presence)
plot(r)
plot(r, col=c("grey","darkgreen"))

# aggregate
ragg<-aggregate(r, fact=10, fun='max')
plot(ragg, col=c("white","darkgreen"))
freq(ragg)
plot(ragg)
plot(r)
print(Sys.time)
