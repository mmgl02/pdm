library(sf)
library(tidyverse)
library(ggplot2)
library(terra) # pour gérer les raster
library(data.table) # pour gérer les gros df
library(patchwork) #plot
library(readxl)
library(corrplot)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Preparation ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list=ls())
# Paths
parameters_file<-file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "../../parameters.R") 
source(parameters_file)
source ("scripts/functions.R")

# Load bioclim variables # fast
vars<- c("bio01d","bio04d", "bio12d", "bio15d", "gdd5")
bioclim<- rast (file.path(path_prep_env, "M1", "bioclim_mean_1992_2020.tif"))

# Load target rast 10 km
target<-rast("data_prep/target/target_rast_10000.tif")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- Associate variables with 10 000 target group points ---- 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
set.seed(123)

# sample 10 000 points from the target raster
sample<-spatSample(target, size = 10000, method = "random", na.rm = TRUE, xy=TRUE) 

# extract values from the bioclim raster
for (varname in names(bioclim)) {
  
  # get the raster
  r <- bioclim[[varname]]
  
  # extract values from the raster
  vals <- terra::extract(r, vect(sample[, c("x", "y")], geom = c("x","y"), crs = "EPSG:3035"))
  
  # add to the sample dataframe
  sample [[varname]]<- vals[,2]
}

# Correlation
vars <- c("gdd5", "bio01d", "bio04d", "bio12d", "bio15d")
cor<-cor(sample[vars], use = "complete.obs", method = "spearman")
corrplot(cor, method='number',  col = rev(COL2("RdBu", 200)))

# histogrammes
par(mfrow = c(2, 2), mar = c(4, 4, 2, 1))
hist(dt_hor$gdd5, main = "gdd5", xlab = "",  breaks = 30)
hist(dt_hor$bio04d, main = "bio04d", xlab = "", breaks = 30)
hist(dt_hor$bio12d, main = "bio12d", xlab = "",  breaks = 30)
hist(dt_hor$bio15d, main = "bio15d", xlab = "",  breaks = 30)
