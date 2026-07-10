library(tidyverse)
library(ggplot2)
library(terra)
library(sf)
library(data.table)
library(viridis)
library(patchwork)
library(tidyterra)
library(ggspatial) # anootation scale

library(xtable)# save the figures and table in documents, export the tables in latex

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Preparation ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Paths
parameters_file<-file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "../../parameters.R") 
source(parameters_file)
source ("scripts/functions.R")

# Load grids
grid_10000<-rast (file.path(path_raw_gridded,"grid_10000_eu.tif")) 
grid_1000<-rast (file.path(path_raw_gridded,"grid_1000_eu.tif")) 
grid_50000<- rast(file.path(path_raw_gridded,"grid_50000_eu.tif"))
grid_50000_poly <- as.polygons(grid_50000) %>% st_as_sf()
grid_10000_poly <- as.polygons(grid_10000) %>% st_as_sf()


# Species
species_names<-c( "Bombus_hortorum","Bombus_pascuorum","Bombus_terrestris", "Bombus_alpinus", "Bombus_pyrenaeus", "Bombus_wurflenii")

species_labels <- c(
  Bombus_hortorum   = "B. hortorum",
  Bombus_pascuorum  = "B. pascuorum",
  Bombus_terrestris = "B. terrestris",
  Bombus_alpinus    = "B. alpinus",
  Bombus_pyrenaeus  = "B. pyrenaeus",
  Bombus_wurflenii  = "B. wurflenii"
)

species_colors <- c(
  Bombus_hortorum   = "#E0313D",   #"#31E083",
  Bombus_pascuorum  = "#E0D431",
  Bombus_terrestris = "#E0A031",
  Bombus_alpinus    = "#31E0D4",
  Bombus_pyrenaeus  = "#3194E0",
  Bombus_wurflenii  = "#313DE0"
)

# Models
models_labels <- c("M 10km", "M 1km", "Nested PR", "Nested GM")
models<-c("HS1", "HS2", "HS_PR", "HS_multiply")
models_corresp <- setNames(models_labels, models)

# extent for the Alps
ext_alp<-ext(c(3.9e+6, 4.5e+6, 2.3e+6,2.8e+6))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 1a. Observation effort per species ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load observations
## get species files
species_files <- list(
  "Bombus hortorum"   = file.path(path_prep_obs, "data_gbif_clean.csv") ,
  "Bombus pascuorum"  = file.path(path_prep_obs, "data_gbif_clean.csv"),
  "Bombus terrestris" = file.path(path_prep_obs, "data_gbif_clean.csv"),
  "Bombus alpinus"    = file.path(path_prep_obs, "balpins_data_gbif_clean.csv"),
  "Bombus pyrenaeus"  = file.path(path_prep_obs, "balpins_data_gbif_clean.csv"),
  "Bombus wurflenii"  = file.path(path_prep_obs, "balpins_data_gbif_clean.csv") 
)

## get each filename once
unique_files <- unique(unlist(species_files))

## Load each obs file and filter in 1992-2020
obs_list <- lapply(unique_files, function(f) {
  fread(f)[year >= 1992 & year <= 2020]
})
names(obs_list) <- unique_files # associate names

## Combine observations for all species into one datatable
obs_dt_all <- rbindlist(lapply(names(species_files), function(sp) {
  obs_list[[species_files[[sp]]]][species == sp] # get species file, and filter it by this species to keep only its observations
}))

## Number of observations per species 
obs_dt_all %>%
  group_by(species) %>%
  summarise(count = n()) %>%
  print()

## Rename species
obs_dt_all[, species := gsub(" ", "_", species)]

# 10km uncertainty in 50k cell
obs_dt_10000<-obs_dt_all[coordinateUncertaintyInMeters<=10000]

obs_rast_10000_list <- lapply(species_names, function(sp) {
  dt_sp <- obs_dt_10000[species == sp]
  pts <- vect(dt_sp[, .(decimalLongitude, decimalLatitude)],
                   geom = c("decimalLongitude", "decimalLatitude"),
                   crs = "EPSG:4326")
  pts <- project(pts, "EPSG:3035")
  rasterize(pts, grid_50000, field = 1, fun = "sum", background = NA)
})
names(obs_rast_10000_list) <- species_names


# 1 km uncertainty in 50k cell
obs_dt_1000<-obs_dt_all[coordinateUncertaintyInMeters<=1000]

obs_rast_1000_list <- lapply(species_names, function(sp) {
  dt_sp <- obs_dt_1000[species == sp]
  pts <- vect(dt_sp[, .(decimalLongitude, decimalLatitude)],
              geom = c("decimalLongitude", "decimalLatitude"),
              crs = "EPSG:4326")
  pts <- project(pts, "EPSG:3035")
  rasterize(pts, grid_50000, field = 1, fun = "sum", background = NA)
})
names(obs_rast_1000_list) <- species_names


# Ratio 1k abundance/ 10 k abundance
ratio<-function(r1k, r10k) {
  ifel(is.na(r1k) & !is.na(r10k), 0, r1k / r10k) # if there is no observation at 1k but there is at 10k, put 0, otherwise it would have been NA
}
obs_rast_ratio_1k_10k_list <- mapply(ratio, obs_rast_1000_list,obs_rast_10000_list,SIMPLIFY = FALSE) 

# Stack
ratio_stack <- rast(obs_rast_ratio_1k_10k_list)

# Plot ratio
grid_50000_poly <- as.polygons(grid_50000) %>% st_as_sf()

p <- ggplot() +
  geom_sf(data = grid_50000_poly, fill = "#E9E9F2", color = NA) +
  geom_spatraster(data = ratio_stack, maxcell = 5e7) +
  scale_fill_viridis_c(
    option = 'magma', end = 0.8, na.value = "transparent",
    name = expression(paste("Ratio ", N[1*km] / N[10*km]))
  ) +
  facet_wrap(~lyr, nrow = 3, ncol = 2, labeller = as_labeller(species_labels)) +
  coord_sf(
    xlim = c(xmin(grid_50000), xmax(grid_50000)),
    ylim = c(ymin(grid_50000), ymax(grid_50000)),
    crs = st_crs(3035)
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "italic", size = 12),  
    axis.text  = element_text(size = 8),
    legend.title = element_text(size = 12)
  )
p

# save
ggsave(filename = file.path(path_results, "obs/obs_ratio_1k_10k.pdf"), 
       plot = p,  dpi = 300,width=8, device="pdf")

#test
obs_dt_1000%>%
  group_by(species)%>%
  summarise(count = n())%>%
  print()

obs_dt_10000%>%
  group_by(species)%>%
  summarise(count = n())%>%
  print()

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 1b. Observation effort target group  ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#raster 10k
target_rast_10000<-rast( "data_prep/target/target_rast_10000.tif")

#raster 1km uncertainty aggregated in 10km cells
# Load obs if needed
target_dt <- fread("data_prep/target/targetgroup_gbif_clean.csv")

# Select obs
target_dt_1000<-target_dt[coordinateUncertaintyInMeters<=1000 & year>=1992 & year<=2020]

# Vector
pts<-vect(target_dt_1000[, .(decimalLongitude, decimalLatitude)],
          geom = c("decimalLongitude","decimalLatitude"))
length(pts)
pts<-project(pts, "EPSG:3035")

# Rasterize
target_rast_1000<-rasterize(pts, grid_10000, field = 1, fun = "sum", background = NA) %>% # just project in the grid extent even where grid=NA
  mask(grid_10000) # to have the same NA cells as the grid raster

# Ratio
target_ratio<-ratio(target_rast_1000, target_rast_10000)

# absolute 10k
p1<-ggplot() +
  geom_sf(data = grid_10000_poly, fill = "#E9E9F2", color = NA) +
  geom_spatraster(data = target_rast_10000, maxcell = 5e7) +
  scale_fill_viridis_c(
    option = 'viridis', na.value = "transparent",
    trans = "log10",
    labels = scales::label_number(),
    name = expression(paste("Observation effort ", N[10*km]))
  ) +
  coord_sf(
    xlim = c(xmin(grid_50000), xmax(grid_50000)),
    ylim = c(ymin(grid_50000), ymax(grid_50000)),
    crs = st_crs(3035)
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.2),
    strip.text = element_text(face = "italic", size = 12),  
    axis.text  = element_text(size = 8),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 12)
  )+
  labs(title = "(a)")

# absolute 1k
p2<-ggplot() +
  geom_sf(data = grid_10000_poly, fill = "#E9E9F2", color = NA) +
  geom_spatraster(data = target_rast_1000, maxcell = 5e7) +
  scale_fill_viridis_c(
    option = 'viridis', na.value = "transparent",
    trans = "log10",
    labels = scales::label_number(),
    name = expression(paste("Observation effort ", N[1*km]))
  ) +
  coord_sf(
    xlim = c(xmin(grid_50000), xmax(grid_50000)),
    ylim = c(ymin(grid_50000), ymax(grid_50000)),
    crs = st_crs(3035)
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.2),
    strip.text = element_text(face = "italic", size = 12),  
    axis.text  = element_text(size = 8),
    axis.text.x = element_blank(),
    legend.title = element_text(size = 12)
  )+
  labs(title = "(b)")

# ratio
p3 <- ggplot() +
  geom_sf(data = grid_10000_poly, fill = "#E9E9F2", color = NA) +
  geom_spatraster(data = target_ratio, maxcell = 5e7) +
  scale_fill_viridis_c(
    option = 'magma', end = 0.8, na.value = "transparent",
    name = expression(paste("Ratio ", N[1*km] / N[10*km]))
  ) +
  coord_sf(
    xlim = c(xmin(grid_50000), xmax(grid_50000)),
    ylim = c(ymin(grid_50000), ymax(grid_50000)),
    crs = st_crs(3035)
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.2),
    strip.text = element_text(face = "italic", size = 12),  
    axis.text  = element_text(size = 8),
    legend.title = element_text(size = 12)
  )+
  labs(title = "(c)")

p<-p1/p2/p3

# Save
ggsave(filename = file.path(path_results, "obs/target_obs_effort.pdf"), 
       plot = p,  dpi = 300,width=7, height=10, device="pdf") #the more we reduce widht and height, the bigger the character size

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 1c. Observations for framework figure ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Run beginning of 1a until: names(obs_rast_1000_list) <- species_names

# Binarize obs_rast_1000
r<-as.factor(ifel(is.na(obs_rast_10000_list[[1]]), NA, 1))

p1<-ggplot() +
  geom_sf(data = grid_10000_poly, fill = "#E9E9F2", color = NA) +
  geom_spatraster(data = r, maxcell = 5e7) +
  scale_fill_manual(
    values = c("darkgreen"),
    na.value = "transparent",
    guide = "none"  # retirer si vous voulez garder la légende
  ) +
  coord_sf(
    xlim = c(xmin(grid_50000), xmax(grid_50000)),
    ylim = c(ymin(grid_50000), ymax(grid_50000)),
    crs = st_crs(3035)
  ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "grey", fill = NA, linewidth = 0.2),
    strip.text = element_text(face = "italic", size = 12),  
    axis.text  = element_text(size = 8),
    legend.title = element_text(size = 12)
  )
p1

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 2a. Species occurences map Europe ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1 plot per specie: left panel a): obs 10k gridded 10k/ right panel: obs 1k gridded 10k, with pseudo absences
dt_10000_background <- get(load("data_prep/obs/M1/dt_10000_background.RData"))
dt_1000_background  <- get(load("data_prep/obs/M2/dt_1000_background.RData"))

for (sp in species_names) {
  
  dt_10000_sp <- dt_10000_background[dt_10000_background$species == sp, ]
  dt_1000_sp  <- dt_1000_background[dt_1000_background$species == sp, ]
  
  # rasterize in 50m grid and convert into factor
  r_10k <- as.factor(rasterize(dt_10000_sp, grid_50000, values = dt_10000_sp$presence, fun = 'max'))
  r_1k  <- as.factor(rasterize(dt_1000_sp,  grid_50000, values = dt_1000_sp$presence,  fun = 'max'))
  
  # Panel a) 10 km obs gridded to 50 km
  p_a <- ggplot() +
    geom_spatraster(data = r_10k) +
    scale_fill_manual(
      values = c("0" = "grey", "1" = "darkgreen"),
      labels = c("0" = "Pseudo-absence", "1" = "Presence"),
      na.value = "transparent",
      na.translate = FALSE,
      name = NULL   
    ) +
    labs(title = "i) 10 km resolution") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          plot.title=element_text(hjust=0.5, size=12, face="bold"))
  
  # Panel b) 1 km obs gridded to 50 km
  p_b <- ggplot() +
    geom_spatraster(data = r_1k) +
    scale_fill_manual(
      values = c("0" = "grey", "1" = "darkgreen"),
      labels = c("0" = "Pseudo-absence", "1" = "Presence"),
      na.value = "transparent",
      na.translate = FALSE,
      name = NULL
    ) +
    labs(title = "ii) 1 km resolution") +
    theme_minimal() +
    theme(panel.grid = element_blank(),axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
          plot.title=element_text(hjust=0.5, size=12, face="bold"))
  
  # Combine
  p_combined <- (p_a + p_b) +
    plot_layout(guides = "collect")
  
  ggsave(
    filename = file.path(path_results_obs, paste0("species/obs_raster_", sp, ".png")),
    plot = p_combined,
    width = 12, height = 6, dpi = 300
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 2b. Species occurences map Alps ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1 plot per specie: left panel a): obs 10k gridded 10k/ right panel: obs 1k gridded 10k, with pseudo absences
dt_10000_background <- get(load("data_prep/obs/M1/dt_10000_background.RData"))
dt_1000_background  <- get(load("data_prep/obs/M2/dt_1000_background.RData"))

for (sp in species_names) {
  
  dt_10000_sp <- dt_10000_background[dt_10000_background$species == sp, ]
  dt_1000_sp  <- dt_1000_background[dt_1000_background$species == sp, ]
  
  # rasterize in 10 km grid and convert into factor
  r_10k <- as.factor(rasterize(dt_10000_sp, grid_10000, values = dt_10000_sp$presence, fun = 'max'))
  r_1k  <- as.factor(rasterize(dt_1000_sp,  grid_10000, values = dt_1000_sp$presence,  fun = 'max'))
  
  # Crop alp extent
  r_10k <- crop(r_10k, ext_alp)
  r_1k  <- crop(r_1k,  ext_alp)
  
  # Panel a) 10 km obs gridded 
  p_a <- ggplot() +
    geom_spatraster(data = r_10k) +
    scale_fill_manual(
      values = c("0" = "grey", "1" = "darkgreen"),
      labels = c("0" = "Pseudo-absence", "1" = "Presence"),
      na.value = "transparent",
      na.translate = FALSE,
      name = NULL   
    ) +
    coord_sf(expand = FALSE) + # pour avoir le cadre direct autour du raster sans espace
    labs(title = "i) 10 km resolution") +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          plot.title=element_text(hjust=0.5, size=12, face="bold"),
          panel.border = element_rect(color = "black", size = 0.1))
  
  # Panel b) 1 km obs gridded
  p_b <- ggplot() +
    geom_spatraster(data = r_1k) +
    scale_fill_manual(
      values = c("0" = "grey", "1" = "darkgreen"),
      labels = c("0" = "Pseudo-absence", "1" = "Presence"),
      na.value = "transparent",
      na.translate = FALSE,
      name = NULL
    ) +
    coord_sf(expand = FALSE) +
    labs(title = "ii) 1 km resolution") +
    theme_minimal() +
    theme(panel.grid = element_blank(),axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
          plot.title=element_text(hjust=0.5, size=12, face="bold"),
          panel.border = element_rect(color = "black", size = 0.1))
  
  # Combine
  p_combined <- (p_a + p_b) +
    plot_layout(guides = "collect")
  
  ggsave(
    filename = file.path(path_results_obs, paste0("species/obs_raster_alp_", sp, ".png")),
    plot = p_combined,
    width = 12, height = 6, dpi = 300
  )
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 3. Obs count between resolutions table ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
obs_count <- obs_dt_all[,.(
  n_50000 = sum(coordinateUncertaintyInMeters <= 50000),
  n_10000 = sum(coordinateUncertaintyInMeters <= 10000),
  n_1000  = sum(coordinateUncertaintyInMeters <= 1000),
  n_100   = sum(coordinateUncertaintyInMeters <= 100)
),
by = species
]

obs_count<-obs_count [, ":="(
  ratio_10000_50000 = n_10000 / n_50000,
  ratio_1000_50000 = n_1000 / n_50000,
  ratio_100_50000 = n_100 / n_50000) ]
#save (obs_count , file = "results/1_obs/obs_count.RData")
obs_count 

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 4a. Projections disagregated 10 km ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to plot without axes
plotting <- function(rast, title, show_scale = FALSE) {
  p<-ggplot() +
    geom_spatraster(data = rast, maxcell=5e7) +
    scale_fill_viridis_c(name = "Habitat suitability", na.value = "transparent", limits=c(0,1)) +  # Titre de la légende
    theme_minimal() +
    theme(plot.title=element_text(hjust=0.5, size=12, face="bold"),
          panel.background = element_blank(),  # Fond transparent
          panel.grid = element_blank(),
          axis.text = element_blank(),          # Supprime les textes des axes
          axis.ticks = element_blank(),        # Supprime les ticks
          axis.title = element_blank(),       # Supprime les titres des axes
          legend.title = element_text(size = 12),   # Style du titre de la légende
          panel.border = element_rect(color = "black", size = 0.1))+
    
    labs(title=title)
  
  p
}
# function to plot with axes
plotting_alp <- function(rast, title, show_scale = FALSE, axisx = FALSE, axisy = FALSE) {
  p <- ggplot() +
    geom_spatraster(data = rast, maxcell = 5e7) +
    scale_fill_viridis_c(name = "Habitat suitability", na.value = "transparent", limits = c(0, 1)) +
    theme_minimal() +
    theme(plot.title    = element_text(hjust = 0.5, size = 12, face = "bold"),
          panel.background = element_blank(),
          panel.grid    = element_blank(),
          axis.text     = element_blank(),    
          axis.ticks    = element_blank(),
          axis.title    = element_blank(),
          legend.title  = element_text(size = 12))+ #,
          #panel.border = element_rect(color = "black", linewidth = 0.1)) +
    labs(title = title)
  
  if (axisx) {
    p <- p + theme(axis.text.x = element_text(size = 7), axis.ticks.x = element_line())
  }
  if (axisy) {
    p <- p + theme(axis.text.y = element_text(size = 7), axis.ticks.y = element_line())
  }
  
  p
}

comparing_projections_eu<-function(sp){
  cat("Comparing projections for species:", sp, "/n")
  
  # Load projections
  rast_M1<-rast(file.path(path_results_models, gsub ("_", ".", sp), "proj_M1/individual_projections", 
                          paste0(gsub ("_", ".", sp), "_EMwmeanByAUCroc_mergedData_mergedRun_mergedAlgo.tif")))%>%
    '/' (1000)
  
  rast_M2<-rast(file.path(path_results_models, gsub ("_", ".", sp), "proj_M2_1_2019_Eu/individual_projections", 
                          paste0(gsub ("_", ".", sp), "_EMwmeanByAUCroc_mergedData_mergedRun_mergedAlgo.tif")))%>%
    '/' (1000)
  # crop if needed
  
  # Disaggregate rast 10km
  rast_M1_disagg<-disagg(rast_M1, fact=c(10,10), method="near")%>% # or resample but takes more time
    crop(rast_M2)
  
  # Good extent
  rast_M2<-crop(rast_M2, rast_M1_disagg)
  
  # Combine with PR
  rast_PR<-1/ (  1 + ( (1-rast_M1_disagg)/rast_M1_disagg ) * ((1-rast_M2)/rast_M2) )
  
  # Combine with multiply
  rast_multiply<-sqrt(rast_M1_disagg * rast_M2)
  
  rast_M2_agg<-aggregate(rast_M2, fact=c(10,10), fun="mean", na.rm=TRUE)
  rast_PR_agg<-aggregate(rast_PR, fact=c(10,10), fun="mean", na.rm=TRUE)
  rast_multiply_agg<-aggregate(rast_multiply, fact=c(10,10), fun="mean", na.rm=TRUE)
  
  
  p1 <- plotting_alp(rast_M1, "M 10 km - 1992-2020", axisy = TRUE)
  p2 <- plotting_alp(rast_M2_agg, "M 1 km - 2019")
  p3 <- plotting_alp(rast_PR_agg, "Nested PR", axisx = TRUE, axisy = TRUE)
  p4 <- plotting_alp(rast_multiply_agg, "Nested GM", show_scale=TRUE, axisx = TRUE)
  
  final_plot <- wrap_plots(p1, p2, p3, p4, ncol = 2, guides = "collect") +
    theme(legend.position = "right") +
    annotation_scale(data = rast_M1, location = "br")
  
  ggsave(filename = file.path(path_results, paste0("models/all/projections/without_border/",sp,".pdf")), 
         plot = final_plot,  dpi = 300,width=10, height=8, device="pdf")

}

for (sp in species_names) {
  comparing_projections_eu(sp)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 4b. Projections alps----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

comparing_projections_alp<-function(sp){
  cat("Comparing projections for species:", sp, "/n")
  
  # Load projections
  rast_M1<-rast(file.path(path_results_models, gsub ("_", ".", sp), "proj_M1/individual_projections", 
                          paste0(gsub ("_", ".", sp), "_EMwmeanByAUCroc_mergedData_mergedRun_mergedAlgo.tif")))%>%
    '/' (1000)
  
  rast_M2<-rast(file.path(path_results_models, gsub ("_", ".", sp), "proj_M2_1_2019_Eu/individual_projections", 
                          paste0(gsub ("_", ".", sp), "_EMwmeanByAUCroc_mergedData_mergedRun_mergedAlgo.tif")))%>%
    '/' (1000)
  # crop if needed

  
  # Disaggregate rast 10km
  rast_M1_disagg<-disagg(rast_M1, fact=c(10,10), method="near")%>% # or resample but takes more time
    crop(ext_alp)
  
  # Good extent
  rast_M2<-crop(rast_M2, ext_alp)
  rast_M1<-crop(rast_M1, ext_alp)
  
  # Combine with PR
  rast_PR<-1/ (  1 + ( (1-rast_M1_disagg)/rast_M1_disagg ) * ((1-rast_M2)/rast_M2) )
  
  # Combine with multiply
  rast_multiply<-sqrt(rast_M1_disagg * rast_M2)
  
  p1 <- plotting_alp(rast_M1, "M 10 km - 1992-2020", axisy = TRUE)
  p2 <- plotting_alp(rast_M2, "M 1 km - 2019")
  p3 <- plotting_alp(rast_PR, "Nested PR", axisx = TRUE, axisy = TRUE)
  p4 <- plotting_alp(rast_multiply, "Nested GM", axisx = TRUE)
  
  final_plot <- wrap_plots(p1, p2, p3, p4, ncol = 2, guides = "collect") +
    theme(legend.position = "right")
  
  ggsave(filename = file.path(path_results, paste0("models/all/projections/",sp,"alp.pdf")), 
         plot = final_plot,  dpi = 300,width=10,height=8, device="pdf")
}

for (sp in species_names) {
  comparing_projections_alp(sp)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 4c. Projections alps with presence points ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# function to plot with axes
plotting_alp <- function(rast, title, show_scale = FALSE, axisx = FALSE, axisy = FALSE) {
  p <- ggplot() +
    geom_spatraster(data = rast, maxcell = 5e7) +
    scale_fill_viridis_c(name = "Habitat suitability", na.value = "transparent", limits = c(0, 1)) +
    theme_minimal() +
    theme(plot.title    = element_text(hjust = 0.5, size = 12, face = "bold"),
          panel.background = element_blank(),
          panel.grid    = element_blank(),
          axis.text     = element_blank(),    
          axis.ticks    = element_blank(),
          axis.title    = element_blank(),
          legend.title  = element_text(size = 12))+ #,
    #panel.border = element_rect(color = "black", linewidth = 0.1)) +
    labs(title = title)
  
  if (axisx) {
    p <- p + theme(axis.text.x = element_text(size = 7), axis.ticks.x = element_line())
  }
  if (axisy) {
    p <- p + theme(axis.text.y = element_text(size = 7), axis.ticks.y = element_line())
  }
  
  p
}

comparing_projections_alp<-function(sp){
  cat("Comparing projections for species:", sp, "/n")
  
  # Load projections
  rast_M1<-rast(file.path(path_results_models, gsub ("_", ".", sp), "proj_M1/individual_projections", 
                          paste0(gsub ("_", ".", sp), "_EMwmeanByAUCroc_mergedData_mergedRun_mergedAlgo.tif")))%>%
    '/' (1000)
  
  rast_M2<-rast(file.path(path_results_models, gsub ("_", ".", sp), "proj_M2_1_2019_Eu/individual_projections", 
                          paste0(gsub ("_", ".", sp), "_EMwmeanByAUCroc_mergedData_mergedRun_mergedAlgo.tif")))%>%
    '/' (1000)
  # crop if needed
  
  
  # Disaggregate rast 10km
  rast_M1_disagg<-disagg(rast_M1, fact=c(10,10), method="near")%>% # or resample but takes more time
    crop(ext_alp)
  
  # Good extent
  rast_M2<-crop(rast_M2, ext_alp)
  rast_M1<-crop(rast_M1, ext_alp)
  
  # Combine with PR
  rast_PR<-1/ (  1 + ( (1-rast_M1_disagg)/rast_M1_disagg ) * ((1-rast_M2)/rast_M2) )
  
  # Combine with multiply
  rast_multiply<-sqrt(rast_M1_disagg * rast_M2)
  
  p1 <- plotting_alp(rast_M1, "M 10 km - 1992-2020", axisy = TRUE)
  p2 <- plotting_alp(rast_M2, "M 1 km - 2019")
  p3 <- plotting_alp(rast_PR, "Nested PR", axisx = TRUE, axisy = TRUE)
  p4 <- plotting_alp(rast_multiply, "Nested GM", axisx = TRUE)
  
  final_plot <- wrap_plots(p1, p2, p3, p4, ncol = 2, guides = "collect") +
    theme(legend.position = "right")
  
  ggsave(filename = file.path(path_results, paste0("models/all/projections/with_presences/",sp,"alp.pdf")), 
         plot = final_plot,  dpi = 300,width=10,height=8, device="pdf")
}

for (sp in species_names) {
  comparing_projections_alp(sp)
}
