library(tidyverse)
library(ggplot2)
library(biomod2)
library(blockCV)
library(terra)
library(sf)
library(data.table)
library(modEvA) # Boyce # before
library(pROC) #AUC
library(viridis)
library(ecospat) #TSS sensitivity specificity

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Preparation ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
gc()
# Paths
parameters_file<-file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "../../parameters.R") 
source(parameters_file)
source ("scripts/functions.R")

# Load grids
grid_10000<-rast (file.path(path_raw_gridded,"grid_10000_eu.tif")) 
grid_1000<-rast (file.path(path_raw_gridded,"grid_1000_eu.tif")) 
#grid_50000<- rast("data_raw/Gridded/ReferenceGrid_Europe_bin_50000m.tif") # to plot

# France
france_3035<-vect(file.path(path_raw_mask_france, "ARRONDISSEMENT.shp"))%>% # from  https://geoservices.ign.fr/geofla # EPSG 2154 Lambert 93
  aggregate () %>%
  project("EPSG: 3035")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Load observations and background---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dt_10000_background<-get(load("data_prep/obs/M1/dt_10000_background.RData"))
dt_1000_background<-get(load("data_prep/obs/M2/dt_1000_background.RData"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Load variables and associate to the presence-background dataset ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Variables
vars_M1<-setdiff(names(dt_10000_background), c("species", "lon", "lat", "presence"))
vars_M2<-setdiff(names(dt_1000_background), c("species","lon", "lat", "year", "presence"))
vars_M1_M1<-paste0(vars_M1, "_M1")

# Load
bioclim_M1<-rast (file.path(path_prep_env, "M1", "bioclim_mean_1992_2020.tif")) [[vars_M1]]

years<-1992:2020

bioclim_M2_list <- lapply(years, function(y){
  if (y==2020){vars_M1<-setdiff(vars_M1, "bio15d")} # bio15d missing in 2020
  rast(paste0("data_prep/env/M2/bioclim_", y,".tif")) [[vars_M1]]} )  %>% # do not select bio01d
  setNames(years)

LC_M2_list <- lapply(years, function(y){
  rast(paste0("data_prep/env/M2/LC_", y,".tif"))} )%>%
  setNames(years)

## Associate
dt_M1_M2 <- copy(dt_1000_background)
for (var in vars_M1) {
  
  r <- bioclim_M1[[var]]
  
  # extract values from the raster
  vals <- terra::extract(r, vect(dt_M1_M2[, .(lon, lat)], geom = c("lon","lat"), crs = "EPSG:3035"))
  
  # add to dt
  dt_M1_M2[, (paste0(var, "_M1")) := vals[,2]]
}  

save (dt_M1_M2, file="data_prep/obs/dt_M1_M2.RData")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 1. Predict the 1k dataset with the 5-CV models ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load ("data_prep/obs/dt_M1_M2.RData")

predicting_for_evaluation<-function(sp){
  cat("Predicting for species:", sp, "/n")
  
  # Select species in dataset
  dt_sp<-dt_M1_M2[species==sp]
  
  # Load models
  M1<-get(load( file.path(path_results_models, gsub ("_", ".", sp), paste0(gsub ("_", ".", sp), ".M1_1.models.out"))))
  M2<-get(load( file.path(path_results_models, gsub ("_", ".", sp), paste0(gsub ("_", ".", sp), ".M2_1.models.out"))))
  M1_models<-get_built_models(M1) # returns the names of the models built for each fold, e.g. _allData_RUN1_RF
  
  # Get calibration table i.e. the calibration-test folds of M2 (obs 1k on which the combined model is evaluated)
  cal_M2<-get_calib_lines(M2)
  
  ## Predict HS2
  # Predict HS2 for each obs
  pred_M2<-get_predictions(M2, model.as.col=TRUE) # returns a dt with line=obs and col= CV model (*5)
  
  # Set NA where obs used in calibration set
  pred_M2[cal_M2] <-NA
  
  # Get the prediction where obs in the test set (only one per row)
  dt_sp$HS2 <- rowMeans(pred_M2, na.rm = TRUE) /1000
  
  ## Predict HS1
  dt_sp$HS1<-NA
  # Data frame of calibration lines
  cal_M2_df <- as.data.frame(cal_M2)  # [n_obs_1k × 5] 
  
  # For each run (each 5-CV fold)
  for (run in 1:5) {

    # Rename vars M1
    env_M1 <- dt_sp[, ..vars_M1_M1] %>%
      setnames(vars_M1_M1, vars_M1)
    
    # Predict with the M1 model of the corresponding run 
    proj_M1 <- BIOMOD_Projection(
      bm.mod              = M1,
      proj.name           = paste0(sp, "_HS1_RUN", run),
      new.env             = env_M1,
      new.env.xy          = dt_sp[, .(lon, lat)],
      models.chosen       = M1_models[run],
      build.clamping.mask = FALSE,
      overwrite = TRUE
    )
    
    # Update HS1 for the test set observations of the run
    run_col  <- paste0("_allData_RUN", run)
    idx_test <- which(!cal_M2_df[, run_col])  # obs 1k in test set for this run: cal=FALSE
    
    HS1 <- get_predictions(proj_M1)$pred / 1000 # in fact, for the Boyce index as it is a correlation it does not change anything to divide by 1000 or not
    HS1<-HS1[idx_test]
    dt_sp$HS1[idx_test] <- HS1
  }
  
  # Combine Permanence of ratios
  dt_sp$HS_PR<-1/ (  1 + ( (1-dt_sp$HS1)/dt_sp$HS1 ) * ((1-dt_sp$HS2)/dt_sp$HS2) )
  
  # Combine multiply
  dt_sp$HS_multiply<-sqrt(dt_sp$HS1 * dt_sp$HS2)
  
  dir.create(file.path(path_results_models, gsub ("_", ".", sp), "models/M1_M2"), recursive = TRUE, showWarnings = FALSE)
  save (dt_sp, file=file.path(path_results_models, gsub ("_", ".", sp),"models/M1_M2/dt_sp_M1_M2.RData"))
  return (dt_sp)}

for (sp in species_names) {
  cat("Predicting for species:", sp, "/n")
  predicting_for_evaluation(sp)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 2.1. Evaluate on the 1k dataset ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load(file.path(path_results_models, gsub ("_", ".", "Bombus_hortorum"), "models/M1_M2/dt_sp_M1_M2.RData"))
sp<-"Bombus_hortorum"

species_names<-c("Bombus_hortorum", "Bombus_pascuorum", "Bombus_terrestris", "Bombus_alpinus", "Bombus_pyrenaeus", "Bombus_wurflenii")
models<- c("HS1", "HS2", "HS_PR", "HS_multiply")

evaluating<-function(dt_sp, sp){
  cat("Evaluating for species:", sp, "/n")
  
  # Get calibration lines
  M2<-get(load( file.path(path_results_models, gsub ("_", ".", sp), paste0(gsub ("_", ".", sp), ".M2_1.models.out"))))
  cal_M2<-get_calib_lines(M2)
  cal_df <- as.data.frame(cal_M2)
  
  # AUC and Boyce for each run
  eval_runs <- rbindlist(lapply(1:5, function(run) {
    
    run_col  <- paste0("_allData_RUN", run)
    idx_test <- which(!cal_df[, run_col])
    
    obs <- dt_sp$presence[idx_test]
    
    rbindlist(lapply(models, function(pred) {
      
      pred_vals <- dt_sp[idx_test, get(pred)]
      
      # Binary evaluation
      tss<-ecospat.max.tss(Pred=pred_vals, Sp.occ=obs)
      threshold<-tss$max.threshold
      eval<-ecospat.meva.table(Pred=pred_vals, Sp.occ=obs, th=threshold)$EVALUATION_METRICS
  
      data.table(
        model = pred,
        run   = run,
        AUC   = as.numeric(pROC::auc(response = obs, predictor = pred_vals)),
        Boyce = Boyce(obs = obs, pred = pred_vals, plot=FALSE, rm.duplicate=FALSE)$Boyce,
        TSS_max= tss$max.TSS,
        Sensitivity= eval[eval$Metric=="Sensitivity", "Value"],
        Specificity= eval[eval$Metric=="Specificity", "Value"],
        threshold=threshold
      )
    }))
  }))
  print(eval_runs)
  
    # save
    write.csv(eval_runs, file.path(path_results_models, gsub ("_", ".", sp), "models/M1_M2/evaluation_table_M1_M2_all_runs.csv"), row.names = FALSE)
}

# Evaluating
for (sp in species_names) {
  cat("Evaluating species:", sp, "/n")
  dt_sp <- get(load(file.path(path_results_models, gsub ("_", ".", sp), "models/M1_M2/dt_sp_M1_M2.RData")))
  evaluating(dt_sp, sp)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 3. Variable importance ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# method used in biomod: for each model, for each variable, times = 1, maybe répéter 3 fois
#1) shuffle, 
#2) predict with the shuffled variable, 
#3) calculate the correlation between the new predictions and the original predictions. 
#4) importance is 1 - correlation.

predicting_for_varimp <- function(sp) {
  
  dt_sp <- dt_M1_M2[species == sp]
  
  # Load models
  M1 <- get(load(file.path(path_results_models, gsub("_", ".", sp), paste0(gsub("_", ".", sp), ".M1_1.models.out"))))
  M2 <- get(load(file.path(path_results_models, gsub("_", ".", sp), paste0(gsub("_", ".", sp), ".M2_1.models.out"))))
  
  # Calibration lines
  cal_M2_df <- as.data.frame(get_calib_lines(M2))
  
  # Predict HS2 for all runs at once
  pred_M2 <- get_predictions(M2, model.as.col = TRUE)
  
  # Build one dt per run
  dt_sp_training <- lapply(1:5, function(run) {
    
    run_col   <- paste0("_allData_RUN", run)
    idx_train <- which(cal_M2_df[, run_col])  # training obs
    
    # HS2 : prediction of the run model on its training obs
    HS2 <- pred_M2[idx_train, paste0(gsub("_", ".", sp),run_col, "_RF")] / 1000
    
    # HS1 : load projection already calculated for this run
    proj_M1 <- get(load(file.path(path_results_models, gsub("_", ".", sp),
                                  paste0("proj_", sp, "_HS1_RUN", run, "/",
                                         gsub("_", ".", sp), ".", sp, "_HS1_RUN", run, ".projection.out"))))
    HS1 <- get_predictions(proj_M1)$pred / 1000
    HS1 <- HS1[idx_train]
    
    # Combine
    HS_PR       <- 1 / (1 + ((1 - HS1) / HS1) * ((1 - HS2) / HS2))
    HS_multiply <- sqrt(HS1 * HS2)
    
    # Assemble dt with obs, variables and predictions
    dt<- copy(dt_sp[idx_train])
    dt[, `:=`(
      run         = run,
      HS1         = HS1,
      HS2         = HS2,
      HS_PR       = HS_PR,
      HS_multiply = HS_multiply
    )]
  })
  
  names(dt_sp_training) <- paste0("run", 1:5)
  save(dt_sp_training, file=file.path(path_results_models, gsub ("_", ".", sp), "models/M1_M2/dt_sp_training.RData"))
  return(dt_sp_training)
}

for (sp in species_names) {
  cat("Predicting for variable importance for species:", sp, "/n")
  predicting_for_varimp(sp)
}

# Function
get_varimp <-function(var, sp, run, M1, M2, M1_models, M2_models, dt_sp_training, seed) {
  
  # define model
  if (var %in% vars_M1_M1){
    model <- M1
    models<-M1_models
    env_vars<-vars_M1_M1
    mod<-"M1"
    col_unchanged<-"HS2"
  } else { 
      model <- M2
      models<-M2_models
      env_vars<-vars_M2
      mod<-"M2"
      col_unchanged<-"HS1"
      }
  
  # Shuffle
  dt_sp<-dt_sp_training[[run]]
  dt_shuffled <- copy(dt_sp)
  dt_shuffled[[var]] <- sample(dt_shuffled[[var]])

  # Select env variables, and rename if it is M1
  env <- dt_shuffled[, ..env_vars]
  if (mod=="M1"){
    setnames(env, old = vars_M1_M1, new = vars_M1)
  }

  #Predict with the shuffled variable
  proj<-BIOMOD_Projection(
    bm.mod              = model,
    proj.name           = paste0('varimp_', var, '_', run, '_', seed),
    new.env             = env,
    new.env.xy          = dt_shuffled[, c("lon", "lat")],
    models.chosen       = models[run],
    build.clamping.mask = FALSE,
    overwrite=TRUE,
  )
  
  HS_s<-get_predictions(proj)$pred / 1000
  
  # Unchanged col
  HS_unchanged<-dt_shuffled[[col_unchanged]]
  
  # Predict the combination
  HS_PR_s       <- 1 / (1 + ((1 - HS_s) / HS_s) * ((1 - HS_unchanged) / HS_unchanged))
  HS_multiply_s <- sqrt(HS_s * HS_unchanged)
  
  
  # Importance = 1 - cor(original, shuffled)
  dt<-data.table(
    var         = var,
    seed        = seed,
    run=run, 
    HS1         = if (mod=="M1") 1 - cor(dt_sp$HS1,         HS_s, use="complete.obs") else NA, # 2 na for pascu, 1 for alpinus, verify why
    HS2         = if (mod=="M2")  1 - cor(dt_sp$HS2,         HS_s, use="complete.obs") else NA,
    HS_PR       =                     1 - cor(dt_sp$HS_PR,       HS_PR_s, use="complete.obs"),
    HS_multiply =                     1 - cor(dt_sp$HS_multiply, HS_multiply_s, use="complete.obs")
  )
  
  print(paste0("Variable: ", var, " - Seed: ", seed, " done"))
  return (dt)
}

# Function to apply to each species
## Define variables
vars_combined<-c(vars_M1_M1, vars_M2)

getting_variable_importance<-function (sp) {
  cat("Getting variable importance for species:", sp, "/n")
  
  # Load prediction for training
  dt_sp_training<-get(load(file.path(path_results_models, gsub ("_", ".", sp), "models/M1_M2/dt_sp_training.RData")))
  
  # Load models for predicting after shuffling
  M1<-get(load( file.path(path_results_models, gsub ("_", ".", sp), paste0(gsub ("_", ".", sp), ".M1_1.models.out"))))
  M2<-get(load( file.path(path_results_models, gsub ("_", ".", sp), paste0(gsub ("_", ".", sp), ".M2_1.models.out"))))
  M1_models<-get_built_models(M1)
  M2_models<-get_built_models(M2)
  
  results<-list()
  k<-1
  
  # Get variable importance for each variable and each seed)
  for (var in vars_combined){
    for (seed in 1:3){
      set.seed(seed)
      for (run in 1:5){
        results[[k]] <- get_varimp(var, sp, run, M1, M2, M1_models, M2_models, dt_sp_training, seed)
        k <- k + 1
      }
    }
  }
  var_imp_table <- rbindlist(results)
  ## Save
   write.csv(var_imp_table, file.path(path_results_models, gsub ("_", ".", sp), "models/M1_M2/variable_importance_M1_M2.csv"), row.names = FALSE)
  return(var_imp_table)
}


for (sp in species_names) {
  getting_variable_importance(sp)
}

# TRY WITH SHAP to have the direction of change ? NO, because the response curve might be non linear

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 4. Project in Europe ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
projecting_2019<-function(sp, run, ensemble=FALSE){
  cat("Comparing projections for species:", sp, "/n")
  sp_label<-gsub("_", ".", sp)
  
   # Load projections
  if (ensemble==FALSE){
    rast_M1<-rast(paste0("results/models/", sp_label, "/proj_M1/individual_projections/proj_M1_", 
                          sp_label, "_allData_RUN", run, "_RF.tif"))%>%
    '/' (1000)
  
  rast_M2<-rast(paste0( "results/models/", sp_label, "/proj_M2_1_2019_Eu/individual_projections/proj_M2_1_2019_Eu_", 
                        sp_label, "_allData_RUN", run, "_RF.tif"))%>%
    '/' (1000)
  }
 
  if (ensemble== TRUE){
    rast_M1<-rast(paste0("results/models/", sp_label, "/proj_M1/individual_projections/",
                         sp_label, "_EMwmeanByAUCroc_mergedData_mergedRun_mergedAlgo.tif"))  %>%
      '/' (1000)
    rast_M2<- rast(paste0( "results/models/", sp_label, "/proj_M2_1_2019_Eu/individual_projections/",
                         sp_label, "_EMwmeanByAUCroc_mergedData_mergedRun_mergedAlgo.tif"))%>%
      '/' (1000)
  }

  # Disaggregate rast 10km
  rast_M1_disagg<-disagg(rast_M1, fact=c(10,10), method="near")%>% # or resample but takes more time
    crop(rast_M2)
  
  # Good extent by cropping rast M2
  rast_M2<-crop(rast_M2, rast_M1_disagg)
  
  # Combine with PR
  rast_PR<-1/ (  1 + ( (1-rast_M1_disagg)/rast_M1_disagg ) * ((1-rast_M2)/rast_M2) )
  
  # Combine with multiply
  rast_multiply<-sqrt(rast_M1_disagg * rast_M2)
  
  # Save
  dir.create(file.path(path_results_models, sp_label, "models/M1_M2/projections_2019"), recursive = TRUE, showWarnings = FALSE)
  writeRaster(rast_PR, file.path(path_results_models, sp_label, "models/M1_M2/projections_2019", paste0("proj_PR_RUN", run,".tif")), overwrite=TRUE)
  writeRaster(rast_multiply, file.path(path_results_models, sp_label, "models/M1_M2/projections_2019", paste0("proj_multiply_RUN", run,".tif")), overwrite=TRUE)
}

for (sp in species_names) {
  for (run in 1:5) {
    projecting_2019(sp, run)
  }
  projecting_2019(sp, run=0, ensemble=TRUE)
}

# Plot and save
pdf(file.path(path_results_models, gsub ("_", ".", sp), "models/M1_M2/comparison_proj_2019.pdf"))
par(mfrow=c(2,2), oma=c(0,0,2,0)) 
plot(rast_M1, main='HS 1 (1992-2019) ', range=c(0,1))
plot(rast_M2, main='HS 2', range=c(0,1))
plot(rast_PR, main='PR (HS1, HS2)', range=c(0,1))
plot(rast_multiply, main= expression(bold(sqrt('HS1 * HS2'))), range=c(0,1))
title(paste0(sp, " habitat suitability 2019") , outer=TRUE)
dev.off()

# Plot densities
vals <- lapply(list(rast_M1, rast_M2, rast_PR, rast_multiply), function(r) values(r, na.rm=TRUE))
dens <- lapply(vals, density)
cols <- c("red","blue","green","orange")
ylim_max <- max(sapply(dens, function(d) max(d$y))) # to have the highest y-axis limit

pdf(file.path(path_results_models, gsub ("_", ".", sp), "models/M1_M2/density_proj_2019.pdf"))
plot(dens[[1]], col=cols[1], lwd=2, xlim=c(0,1), ylim=c(0,ylim_max), main="Habitat suitability scores (HS) distributions", xlab="HS")
for(i in 2:4) lines(dens[[i]], col=cols[i], lwd=2)
legend("topright", c("HS1","HS2","PR","√(HS1×HS2)"), col=cols, lwd=2, cex=0.8)
dev.off()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 5. Evaluate on the whole study area with Boyce ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
models<-c("HS1", "HS2", "HS_PR", "HS_multiply")
year<- "2019"

evaluating_in_europe <- function(sp, years = year) {
  cat("Evaluating in Europe for species:", sp, "/n")
  sp_label <- gsub("_", ".", sp)
  
  # Select species observations
  dt_sp <- dt_M1_M2[species == sp]
  
  # Calibration lines
  M2 <- get(load(file.path(path_results_models, sp_label, paste0(sp_label, ".M2_1.models.out"))))
  cal_df <- as.data.frame(get_calib_lines(M2))
  
  # Load spatial blocks
  load(paste0("data_prep/blocks/blocks_", sp, ".RData"))
  
  # Rasterize blocks
  rast_blocks_1km  <- rasterize(vect(scv$blocks), bioclim_M2_list[[1]],  field="folds")
  
  # Calculate boyce in the test set for each run and model
  eval_table <- rbindlist(lapply(1:5, function(run) {
    cat("Run", run)
    
    # Y= Presences
    ## Select test set
    run_col  <- paste0("_allData_RUN", run)
    idx_test <- which(!cal_df[, run_col])
    dt_test <- dt_sp[idx_test]
    
    ## Presences in the year in the test set
    obs<-dt_test[presence==1 & year==year, .(lon, lat)]
    
    # ^Y= Predicted raster for each model
    ## test mask
    test_mask <- rast_blocks_1km == run
    
    ## Load raster
    rbindlist(lapply(models, function(model) {
      cat("Model:", model)
      if (model == "HS1") {
        proj_raster <- rast(paste0("results/models/", sp_label, "/proj_M1/individual_projections/proj_M1_", sp_label, "_allData_RUN", run, "_RF.tif"))%>%
        disagg(fact=c(10,10), method="near")%>%
          crop(test_mask)
      } else if (model == "HS2") {
        proj_raster <- rast(paste0("results/models/", sp_label, "/proj_M2_1_2019_Eu/individual_projections/proj_M2_1_2019_Eu_", sp_label, "_allData_RUN", run, "_RF.tif"))
      } else if (model == "HS_PR") {
        proj_raster <- rast(paste0("results/models/", sp_label, "/models/M1_M2/projections_2019/proj_PR_RUN", run, ".tif"))
      } else if (model == "HS_multiply") {
        proj_raster <- rast(paste0("results/models/", sp_label, "/models/M1_M2/projections_2019/proj_multiply_RUN", run, ".tif"))
      }
      
    ## Mask raster to keep only the test set blocks
    test_mask<-crop(test_mask, proj_raster)
    proj_raster_test <-mask(proj_raster, test_mask, maskvalues= FALSE)
    
    ## Boyce
    boyce<-Boyce(obs = obs,  pred=proj_raster_test, plot=FALSE)$Boyce
    
    data.table(
      species = sp,
      model   = model,
      run     = run,
      Boyce   = boyce
    )
    }))
  }))
  print(eval_table)
  
  # save
  write.csv(eval_table, file.path(path_results_models, sp_label, "models/M1_M2/evaluation_table_M1_M2_europe.csv"), row.names = FALSE)
  return(eval_table)
}  
start<-Sys.time()
for (sp in species_names) {
  
  evaluating_in_europe(sp)
}
print(Sys.time()-start)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Apply 1,2,3,4 to each species ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
species_names<-c("Bombus_pascuorum", "Bombus_hortorum","Bombus_terrestris", "Bombus_alpinus", "Bombus_pyrenaeus", "Bombus_wurflenii")

for (sp in species_names) {
  cat("Processing species:", sp, "/n")
  dt_sp <- predicting_for_evaluation(sp)
  evaluating(dt_sp, sp)
  #getting_variable_importance(sp)
  #comparing_projections_2019(sp)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ----Load results ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
load (file.path(path_results_models, gsub ("_", ".", "Bombus_hortorum"), "dt_sp_M1_M2.RData"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- TESTS ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Nb of obs per year  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
obs_per_year<-dt_sp[, .N, by=year][order(year)]
plot(obs_per_year$year, obs_per_year$N, xlab="Year", ylab="Number of observations", main="Observations per year")

# Calibration ok 
cal_M2_df$lon<-dt_sp$lon
cal_M2_df$lat<-dt_sp$lat

plot(cal_M2_df$lon, cal_M2_df$lat, col=ifelse(cal_M2_df$`_allData_RUN2`, "red", "blue"), pch=16, cex=0.5,
     xlab="Longitude", ylab="Latitude", main=paste0("Calibration lines for ", sp))

# plot nombre de données selon l'année
count_year<-dt_M1_M2 %>%
  group_by(species, presence) %>%
  summarise(count = n()) %>%
  print()
ggplot(count_year, aes(x = year, y = count, color = species, group = species)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(x = "Year", y = "Number of records", color = "Species",
       title = "Number of records per year")

# na in dt_M1_M2 per species
na_count<-dt_M1_M2 %>%
  group_by(species) %>%
  summarise(across(everything(), ~ sum(is.na(.)))) %>%
  print()

# nb pres
dt_sp[presence==1, .N]
