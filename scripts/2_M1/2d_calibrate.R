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
grid<-rast (file.path(path_raw_gridded,"grid_10000_eu.tif")) 
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
target<-rast(file.path(path_prep_target, "target_rast_10000.tif"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Load environment ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
bioclim<-rast (file.path(path_prep_env, "M1", "bioclim_mean_1992_2020.tif"))

# Variables to calibrate the model with, selected in c_choose_env
vars<- c("gdd5", "bio04d", "bio12d", "bio15d")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 1. Project species data  ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# all species at the same time

# Associate observations to the grid at resolution 10km in Europe, get a datatable
dt <- obs_to_grid (obs_dt_all, grid, 10000) # from functions.R

# Save
save (dt, file = "data_prep/obs/M1/dt_10000.RData")

# Suppress useless files
rm(obs_list, unique_files, obs_dt_all)
gc()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 2. Select background  ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# all species at the same time

# Get the weighted raster to select background points from the target group observation effort raster
# Do not select background in the seas :  NA in the seas, 0 on the continents
weighted_raster<-cover(target, grid*0)

# Do not select background where there is no environmental value
weighted_raster<-mask(weighted_raster,  bioclim[[1]])

# Add 1 to sample everywhere in the continents even where there is no target group observation
weighted_raster<-weighted_raster+1

# Function to select background points for one species
selecting_background <- function(sp) {
  
  # Select presence coordinates for the species
  dt_sp <- dt[species == sp, .(lon, lat)]
  
  # Number of background points
  nb.background <- ifelse(nrow(dt_sp) < 10000, 10000, nrow(dt_sp))
  
  background_points <- select_background_target(
    presence_coords = dt_sp,
    nb.background   = nb.background,
    weighted_raster = weighted_raster
  ) %>% as.data.table()
  
  # Proportion of background poinst not on a target group observation cell
  proportion <- background_points[sum == 1, .N] / nrow(background_points)
  cat("Species:", sp, "| n presences:", nrow(dt_sp), "| n background:", nb.background,
      "| prop outside target:", round(proportion, 3), "\n")
  
  # Bind background with presences
  rbind(
    dt_sp[, .(lon, lat, presence = 1L, species = sp)],
    background_points[, .(lon = x, lat = y, presence = 0L, species = sp)]
  )
}

# Apply to all species
set.seed(123)
dt_background <- rbindlist(lapply(species_names, selecting_background))
save (dt_background, file = "data_prep/obs/M1/dt_10000_background_ori.RData")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 3. Associate to env  ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
for (var in vars) {
  # Select layer
  r <- bioclim[[var]]
  
  # extract values from the raster
  vals <- terra::extract(r, vect(dt_background[, .(lon, lat)], geom = c("lon","lat"), crs = "EPSG:3035"))
  
  # add to dt
  dt_background[, (var) := vals[,2]]
}

## Test association df # verify with QGIS
# df_sf=st_as_sf(dt_background, coords=c("lon", "lat"), crs= st_crs(3035))
# raster_hor<-rasterize(df_hor_sf, grid, field='gdd5', fun='mean')
# par(mfrow=c(1,2))
# plot(raster_hor)
# plot(bioclim$gdd5)

# Remove NA (FOR THE MOMENT, sinon on doit faire ça avant de selectionner les pseudo-absences, avec un mask de env OK)
dt_background <- dt_background[complete.cases(dt_background)]

dt_background%>%
  group_by(species, presence) %>%
  summarise(count = n()) %>%
  print()

# Save
save (dt_background, file = "data_prep/obs/M1/dt_10000_background.RData")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 4. Block CV for evaluation (function) ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#load("data_prep/obs/M1/dt_10000_background.RData")

# BlockCV: https://cran.r-project.org/web/packages/blockCV/vignettes/tutorial_2.html
# convert to sf
sf_10000<-st_as_sf(dt_background, coords=c("lon", "lat"), crs= st_crs(3035))

# raster mask
crs(grid) <- "EPSG:3035"

# # Stats to decide the size of spatial blocks
# ## estimation of the range with sac of the response variable
# sac1 <- cv_spatial_autocor(x=sf_10000,
#                            column='presence',
#                            plot = FALSE)
# 
# sac1$range #4989293 so min this range 
# ### we could do the same for the predictor but it is too long so first aggregate # 30s
# results_agg<- lapply(results, function(r){
#   aggregate(r, fact=c(10,10), fun='mean')%>% project("EPSG:3035")
# })
# res(results_agg[[1]]) #resolution
# 
# sac2 <- cv_spatial_autocor(r= results_agg,num_sample=5000,
#                            plot = FALSE)
# sac2$range_table
# sac2$plots
# # very few spatial autocorrelation??

# Spatial blocking function
spatial_blocking <- function (sp){
  cat("spatial blocking")
  
  # Filter for one species
  sf_sp <- sf_10000 [sf_10000$species==sp,]
  
  # Create spatial blocks and folds
  scv <- cv_spatial(
    x = sf_sp,
    column = "presence", # the response column (binary)
    r = grid,
    k = 5, # number of folds
    size = 100000, # size of the blocks in metres= 100 km
    selection = "random",
    iteration = 10, # to create folds with balanced records e.g. similar number of presence and absence records in each fold
    progress = TRUE,
    biomod2=TRUE,
    raster_colors = terrain.colors(10, rev = TRUE) # options from cv_plot for a better colour contrast
  )
  
  # Save blocks for evaluation of model 2
  save (scv, file = paste0("data_prep/blocks/blocks_", sp, ".RData"))
  
  # Plot
  p <- cv_plot(scv, r = grid, raster_colors = 'lightgrey')
  svg(file = paste0("data_prep/blocks/blocks_", sp, ".svg"), width  = 16, height = 16) # plus on augmente resol, plus chiffres petits
  print(p) 
  dev.off()
  
  # Format for biomod
  spatial_cv_folds<-scv$biomod_table
  colnames(spatial_cv_folds) <- paste0("_allData_RUN", 1:ncol(spatial_cv_folds))
  
  return(spatial_cv_folds)
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
                                      modeling.id = "M1_1",
                                      CV.strategy = "user.defined",
                                      CV.user.table   = spatial_cv_folds,
                                      OPT.user = myOptions,
                                      metric.eval = c('TSS', 'AUCroc', "BOYCE"),
                                      var.import = 1,
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
  plot_dir<-file.path(myBiomodModelOut@dir.name, myBiomodModelOut@sp.name, "models", myBiomodModelOut@modeling.id, "plots")
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
# ---- 9. Project (function) ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
projecting <- function(sp) {
  cat("Projecting:", sp, "\n")
  
  myBiomodModel <- get(load( file.path(path_results_models, gsub("_", ".", sp), 
                                       paste0(gsub("_", ".", sp), ".M1_1.models.out") )))
  myBiomodEM <- get(load( file.path(path_results_models, gsub("_", ".", sp), 
                                    paste0(gsub("_", ".", sp), ".M1_1.ensemble.models.out") )))
  
  myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModel,
                                    proj.name = "M1",
                                    new.env = bioclim[[vars]],
                                    models.chosen = 'all',
                                    build.clamping.mask = FALSE,
                                    do.stack= FALSE,
                                    keep.in.memory = FALSE,
                                    overwrite=TRUE)
  myBiomodProj_EM<-BIOMOD_EnsembleForecasting(bm.em = myBiomodEM, bm.proj=myBiomodProj)
  
  # Load raster and save it
  plot_dir <- file.path(myBiomodModel@dir.name, myBiomodModel@sp.name, "models", myBiomodModel@modeling.id, "plots")# same as in evaluating
  
  r <- rast(file.path(myBiomodModel@dir.name, myBiomodModel@sp.name, paste0("proj_", myBiomodProj@proj.name), "individual_projections", 
                      paste0(myBiomodModel@sp.name, "_EMwmeanByAUCroc_mergedData_mergedRun_mergedAlgo.tif")))%>%
    mask(grid) # because the env raster are 20k buffered around the grid
  pdf(file.path(plot_dir, "projected_distribution_eu.pdf"))
  plot(r, main = paste("HS1", sp))
  dev.off()
  
  return(myBiomodProj)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Apply 4, 5 and 6, 7, 8, 9  to each species ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

for (sp in species_names) {
  set.seed(123)
  # spatial_cv_folds<-spatial_blocking(sp)
  # 
  # setwd(path_results_models)
  # myBiomodData<-formating(sp)
  # myBiomodModelOut<-calibrating(sp, myBiomodData, spatial_cv_folds)
  # evaluating(sp, myBiomodModelOut)
  # myBiomodEM<-ensemble_modelling(sp, myBiomodModelOut)
  projecting(sp)
  setwd(root)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----Load results ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dt_10000_background<-get(load("data_prep/obs/M1/dt_10000_background.RData"))
grid<-rast (file.path(path_raw_gridded,"grid_10000_eu.tif"))

#load mybiomodEM
myBiomodEM<-get(load( file.path(path_results_models, "Bombus.alpinus/Bombus.alpinus.M1_1.ensemble.models.out") ))

myBiomodModelOut<-get(load( file.path(path_results_models, "Bombus.hortorum/Bombus.hortorum.M1_1.models.out") ))
sp<-"Bombus.hortorum"
evaluating(sp, myBiomodModelOut)

load(file.path("data_prep/blocks/blocks_Bombus_alpinus.RData"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----TESTS  ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
sp<-"Bombus_alpinus"
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

# presence and absences
sf_10000 %>%
  filter(presence == 0) %>%
  group_by(species) %>%
  summarise(count = n())



