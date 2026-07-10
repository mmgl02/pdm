# search for the parameter file in the root. In the root: folder scripts.
path_parameters<-file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "../parameters.R") 
source(path_parameters)

dirs<- c(path_raw_obs, path_raw_bioclim,
         path_prep_obs, path_prep_target, path_prep_species, path_prep_env,
         path_results_obs, path_results_target, path_results_env, path_results_species )

lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE)

