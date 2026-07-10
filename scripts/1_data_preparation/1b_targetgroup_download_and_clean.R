#### Librairies ####

library(sf)
library(tidyverse)
library(ggplot2)
library(rgbif)
library(CoordinateCleaner)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Préparation de l'espace de travail ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Paths
parameters_file<-file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "../../parameters.R") 
source(parameters_file)

# To change
begin_year<-1991 # begin downloading at year superior or equal to
end_year <- 2020 # end downloading at year
max_uncertainty<- 10000 # in meters
filename_raw<- "targetgroup_gbif.csv"
filename_prep<- "targetgroup_gbif_clean.csv"

# Options de gbif
options(gbif_user = "moea")
options(gbif_pwd = "") # add the password
options(gbif_email = "moea.geffardlemaitre@gmail.com")

# Eurostat countries for CoordinateCleaner
EUROSTAT_Countries_1k<- st_read(file.path(path_raw_mask_europe,"EUROSTAT_Countries_buffered1k_WGS84_nng (2024_07_11 13_01_52 UTC).shp"))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 1. Download all families of the target group ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

name <- c("Andrenidae",
          "Apidae",
          "Colletidae",
          "Halictidae",
          "Megachilidae",
          "Melittidae",
          "Stenotritidae"
)

key <- map_chr(
  name,
  ~ as.character(
    name_backbone(name = .x, rank = "family")$usageKey
  )
) %>%
  as.numeric()

# Download
gbif_dl <- occ_download(
  pred_and(
    pred_in("taxonKey", key),
    pred("geometry", "POLYGON((-45 21,65 21,65 73,-45 73,-45 21))"),
    pred("hasCoordinate", TRUE),
    pred("hasGeospatialIssue", FALSE),
    pred("occurrenceStatus","PRESENT"),
    pred("basisOfRecord", "HUMAN_OBSERVATION"),
    pred_gte("year", begin_year),
    pred_lte("year", end_year),
    pred_lte("coordinateUncertaintyInMeters", max_uncertainty) 
    ),
  format = "SIMPLE_CSV")

occ_download_wait(gbif_dl, status_ping = 5) 

dl <- occ_download_get(gbif_dl) %>% #retrieve a download from GBIF to your computer
  occ_download_import() #load a download from your computer to R

# Save
write.csv(dl, file.path(path_raw_obs, filename_raw), row.names = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### 2. Cleaning ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load data
#df <- read.csv(file.path(path_raw_obs, "targetgroup_gbif.csv"), sep="\t")
df<-dl

### Remove useless columns
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
print(n0)
df <- df %>% 
  dplyr::select(-any_of(cols_to_remove)) %>%
  drop_na(year)

n1 <- nrow(df)
print(n1)

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
          file = file.path(path_prep_target, filename_prep), 
          row.names = FALSE)
