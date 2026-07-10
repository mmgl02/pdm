# ---- Libraries ---- 
library(sf)
library(tidyverse)
library(ggplot2)
library(rgbif)
library(data.table)
library(CoordinateCleaner)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Preparation ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Paths
parameters_file<-file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "../../parameters.R") 
source(parameters_file)

# To change
species_to_download<-c("Bombus alpinus", "Bombus pyrenaeus", "Bombus wurflenii" )
filename_raw<- "balpins_data_gbif.csv"
filename_prep <-"balpins_data_gbif_clean.csv"

# GBIF options
options(gbif_user = "moea")
options(gbif_pwd = "") # add the password
options(gbif_email = "moea.geffardlemaitre@gmail.com")

# Eurostat countries for CoordinateCleaner
EUROSTAT_Countries_1k<- st_read(file.path(path_raw_mask_europe,"EUROSTAT_Countries_buffered1k_WGS84_nng (2024_07_11 13_01_52 UTC).shp"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 1. Download studied species ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

name <- species_to_download

key <- map_chr(name, ~ name_backbone(name = .x)[[1]])

gbif_dl <- occ_download(
  pred_and(
    pred_in("taxonKey", key),
    pred("geometry", "POLYGON((-45 21,65 21,65 73,-45 73,-45 21))"),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred("occurrenceStatus","PRESENT"),
    pred("basisOfRecord", "HUMAN_OBSERVATION"),
    pred_gte("year", 1992),
    pred_lte("year", 2025), # we keep only until 2020, but 2025 to validate
    pred_lte("coordinateUncertaintyInMeters", 50000) # finally we only kept the 10k uncertainty obs; filtred in 1b_obs_processing
    ),
  format = "SIMPLE_CSV")

occ_download_wait(gbif_dl, status_ping = 5) 

dl <- occ_download_get(gbif_dl) %>% #retrieve a download from GBIF to your computer
  occ_download_import() #load a download from your computer to R

nrow(dl)

# Save
write.csv(dl, file.path(path_raw_obs, filename_raw), row.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 2. Clean ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Import observations if needed
#df <- read.csv(file.path(path_raw_obs, "balpins_data_gbif.csv"))
df<-dl

### Suppression colonnes inutiles
cols_to_remove <- c(
  "issue", 
  "mediaType", 
  "lastInterpreted", 
  "establishmentMeans", 
  "license",
  "rightsHolder", 
  "identifiedBy", 
  "dateIdentified", 
  "catalogNumber",
  "collectionCode", 
  "institutionCode", 
  "basisOfRecord",
  "infraspecificEpithet",
  "stateProvince",
  "typeStatus",
  "publishingOrgKey",
  "verbatimScientificNameAuthorship")


n0 <- nrow(df)
df <- df %>% 
  dplyr::select(-any_of(cols_to_remove)) %>%
  drop_na(year)

n1 <- nrow(df)

flags_cc <- clean_coordinates(x = df,
                              lon = "decimalLongitude",
                              lat = "decimalLatitude",
                              value = "spatialvalid",
                              verbose = FALSE, species = "species", 
                              seas_ref = EUROSTAT_Countries_1k,capitals_rad=10,
                              tests = c("equal",
                                        "gbif", "capitals",
                                        "institutions", 
                                        "zeros", "
                                       seas"))

print(summary(flags_cc))
df <- df[flags_cc$.summary, ]

write.csv(df, 
          file.path(path_prep_obs, filename_prep), 
          row.names = FALSE)