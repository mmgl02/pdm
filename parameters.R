# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Paths ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Root
root <- "D:/stage/3_code"
setwd(root)

# for the folders and subfolders I sometimes assume that we are already in the root and that we have the same structure

# Folders
path_data_raw     <- file.path(root, "data_raw")
path_data_prep    <- file.path(root, "data_prep")
path_results      <- file.path(root, "results")
path_scripts      <- file.path(root, "scripts")

# Sub folders data_raw
path_raw_obs      <- file.path(path_data_raw, "obs")
path_raw_bioclim  <- file.path(path_data_raw, "bioclim")
path_raw_landcover <- file.path(path_data_raw, "landcover")
path_raw_soil <- file.path(path_data_raw, "soil")
path_raw_mask_europe <- file.path(path_data_raw, "mask_europe")
path_raw_mask_france <- file.path(path_data_raw, "mask_france")
path_raw_gridded <- file.path(path_data_raw, "gridded")
path_raw_target<- file.path(path_data_raw, "target")

# Sub folders data_prep
path_prep_obs     <- file.path(path_data_prep, "obs")
path_prep_target     <- file.path(path_data_prep, "target")
path_prep_species <- file.path(path_data_prep, "species")
path_prep_env      <- file.path(path_data_prep, "env")

# Sub folders results
path_results_obs     <- file.path(path_results, "obs")
path_results_target       <- file.path(path_results, "target")
path_results_env      <- file.path(path_results, "env")
path_results_models <- file.path(path_results, "models")
