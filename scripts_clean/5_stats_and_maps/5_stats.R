library(tidyverse)
library(ggplot2)
library(biomod2)
library(blockCV)
library(terra)
library(sf)
library(data.table)
library(modEvA) # Boyce
library(pROC) #AUC
library(viridis)
library(patchwork)

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
grid_50000<- rast(file.path(path_raw_gridded,"grid_1000_eu.tif"))

# Data
dt_10000_background<-get(load("data_prep/obs/M1/dt_10000_background.RData"))
dt_1000_background<-get(load("data_prep/obs/M2/dt_1000_background.RData"))

# Variables
vars_M1<-setdiff(names(dt_10000_background), c("species", "lon", "lat", "presence"))
vars_M2<-setdiff(names(dt_1000_background), c("species","lon", "lat", "year", "presence"))
vars_M1_M1<-paste0(vars_M1, "_M1")

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

# Variables
vars_colors<-c(
  "Climate (global)" = "#D6BEF9",
  "Climate (annual)"   = "#C113D1",
  "Land cover (annual)"     = "#23D113"
)

vars_shape<-c(
  "Climate (global)" = 15,
  "Climate (annual)"   = 16,
  "Land cover (annual)"     = 17
)

var_labels <- c(gdd5_M1 = "gdd5_g", bio04d_M1="bio04d_g", bio12d_M1="bio12d_g", bio15d_M1="bio15d_g")

# Run
run_colors<<-c(
  "1" = "#F0E911",
  "2"= "#11F079",
  "3"= "#F0A011",
  "4"= "#F01188",
  "5"= "#1118F0"
)


run_shapes<-c(
  "1" = 15,
  "2" = 17,
  "3" = 18,
  "4" = 19,
  "5" = 20)

# Models
models_labels <- c("M 10km", "M 1km", "Nested PR", "Nested GM")
models<-c("HS1", "HS2", "HS_PR", "HS_multiply")
models_corresp <- setNames(models_labels, models)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 1a. Evaluate on the 1k dataset : boxplots ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# all evaluation metrics * model * species * run
all_eval <- rbindlist(lapply(species_names, function(sp) {
  dt <- fread(file.path(path_results_models, gsub("_", ".", sp), "models/M1_M2/evaluation_table_M1_M2_all_runs.csv"))
  dt[, species := sp]
  dt
}))

all_eval$run<-as.factor(all_eval$run)

#rename and set order
all_eval[, model   := factor(model, levels = models, labels = models_labels)]
all_eval[, species := factor(species, levels = species_names, labels = species_labels[species_names])]

# function to plot 
scatterplot <- function(metric, title, show_strip = TRUE) {
  ggplot(all_eval, aes(x = model, y = get(metric), shape= run, col=run)) +
    geom_point(size = 2, alpha = 1) +
    facet_wrap(~ species, scales = "free_y", ncol = 1, strip.position = "left")+
    theme_minimal() +
    scale_shape_manual(values = run_shapes, name="Run") +
    scale_color_manual(values = run_colors, name="Run") +
    theme(
      plot.title         = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.text.x        = element_text(angle = 40, hjust = 1, size=11),
      strip.text.y.left  = if (show_strip) element_text(face = "italic", size = 10, angle = 0) else element_blank(),
      strip.background   = element_blank(),
      strip.placement    = "outside",
      panel.border       = element_rect(color = "black") # rectangle around each panel
    ) +
    scale_y_continuous(breaks = scales::pretty_breaks(n = 3))+ # for 3 graduation in the y scale, does not really work
    labs(x = NULL, y = NULL, title = title)
}

# Plot AUC, Boyce, TSS : noms d'espèces uniquement sur le premier (colonne de gauche)
plots <- list(
  scatterplot("AUC",     "AUC",   show_strip = TRUE),
  scatterplot("Boyce",   "Boyce", show_strip = FALSE),
  scatterplot("TSS_max", "TSS",   show_strip = FALSE)
)

final_plot <- wrap_plots(plots, ncol = 3) +
  plot_layout(guides = "collect")

final_plot

ggsave(filename = file.path(path_results_models, "all/evaluation_scatterplot.pdf"),
       plot     = final_plot, width    = 9, height   = 12,  
       dpi      = 300, device   = "pdf")

# Plot TSS treshold, sensitivity, specificity
plots2 <- list(
  scatterplot("threshold", "Max TSS threshold", show_strip = TRUE),
  scatterplot("Sensitivity",   "Sensitivity",   show_strip = FALSE),
  scatterplot("Specificity",   "Specificity",   show_strip = FALSE)
)
final_plot2<- wrap_plots(plots2, ncol = 3) +
  plot_layout(guides = "collect")

final_plot2

ggsave(filename = file.path(path_results_models, "all/evaluation_scatterplot_TSS_sens_spec.pdf"),
       plot     = final_plot2, width    = 9, height   = 12,  
       dpi      = 300, device   = "pdf")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 1b. Evaluate on the 1k dataset: heatmap ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Moyennes par espèce x modèle pour les 3 métriques (with all_eval_summary from 1a)

# mean and sd per species * model
all_eval_summary <- all_eval[, .(AUC_mean    = mean(AUC, na.rm = TRUE),
                                 AUC_sd      = sd(AUC, na.rm = TRUE),
                                 Boyce_mean  = mean(Boyce, na.rm = TRUE),
                                 Boyce_sd    = sd(Boyce, na.rm = TRUE),
                                 TSS_mean    = mean(TSS_max, na.rm = TRUE),
                                 TSS_sd      = sd(TSS_max, na.rm = TRUE)),
                             by = .(species, model)]

# rename and set order
all_eval_summary[, model   := factor(model, levels = models, labels = models_labels)]
all_eval_summary[, species := factor(species, levels = rev(species_names),
                                     labels = species_labels[rev(species_names)])]

## Function to plot heatmap with a metric
make_heatmap <- function(metric_col, title, show_legend = TRUE, show_y_labels = TRUE) {
  
  d <- all_eval_summary[, .(species, model, value = get(metric_col))]
  
  ggplot(d, aes(x = model, y = species, fill = value)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = round(value, 3)), size = 5) +
    scale_fill_viridis(discrete = FALSE, name = title) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title        = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.text.y       = if (show_y_labels) element_text(face="italic", size= 12) else element_blank(),
      axis.text.x       = element_text(angle = 40, size= 12,hjust = 1),
      panel.grid        = element_blank(),
      legend.position   = if (show_legend) "right" else "none",
      legend.title      = element_text(size = 9)
    ) +
    labs(x = NULL, y = NULL, title = title)
}

# 3 heatmaps : noms d'espèces uniquement sur le premier
heatmap_AUC   <- make_heatmap("AUC_mean",   "AUC",   show_y_labels = TRUE)
heatmap_Boyce <- make_heatmap("Boyce_mean", "Boyce", show_y_labels = FALSE)
heatmap_TSS   <- make_heatmap("TSS_mean",   "TSS",   show_y_labels = FALSE)

# Assemblage en 3 colonnes with patchwork
final_heatmaps <- wrap_plots(heatmap_AUC, heatmap_Boyce, heatmap_TSS, ncol = 3) +
  plot_layout(guides = "collect")

final_heatmaps

ggsave(filename = file.path(path_results_models, "all/evaluation_heatmaps.svg"),
       plot     = final_heatmaps, width    = 12, height   = 6,  
       dpi      = 300, device   = "svg")
ggsave (filename = file.path(path_results_models, "all/evaluation_heatmaps.pdf"),
       plot     = final_heatmaps, width    = 12, height   = 6,  
       dpi      = 300, device   = "pdf")

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 2a. Variable importance per species ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load and plot variable importance for each species

# Panel for one model and one species
make_species_boxplot_panel <- function(var_imp_long, model_name, letter, show_legend = FALSE) {
  
  d <- var_imp_long[model == model_name]
  
  # Order variables by mean importance (décroissant), propre à ce modèle
  var_order <- d[, .(mean_imp = mean(imp, na.rm = TRUE)), by = var][order(mean_imp), var]
  d[, var := factor(var, levels = var_order)]
  
  ggplot(d, aes(x = var, y = imp, color = category, fill = category)) +
    geom_boxplot(alpha = 0.3, outlier.size = 1) +
    coord_flip() +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_discrete(labels = function(x) ifelse(x %in% names(var_labels), var_labels[x], x)) +
    scale_color_manual(values = vars_colors, breaks = names(vars_colors), name = "Category") +
    scale_fill_manual(values = vars_colors, breaks = names(vars_colors), name = "Category") +
    theme_minimal() +
    theme(
      plot.title  = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.text.y = element_text(size = 9),
      plot.tag    = element_text(size = 10),
      plot.tag.position = c(0, 1),
      legend.position = if (show_legend) "right" else "none"
    ) +
    labs(title = models_corresp[[model_name]], x = NULL, y = NULL, 
         color = NULL, fill = NULL, tag = paste0("(", letter, ")"))
}

# for one species
make_species_boxplot <- function(sp) {
  
  # Load var_imp_table
  var_imp_table <- fread(file.path(path_results_models, gsub("_", ".", sp), "models/M1_M2/variable_importance_M1_M2.csv"))
  
  # Long format
  var_imp_long <- melt(var_imp_table, id.vars = c("var", "seed"),
                       measure.vars = c("HS1", "HS2", "HS_PR", "HS_multiply"),
                       variable.name = "model", value.name = "imp")
  
  # assign categories
  var_imp_long[, category := fcase(
    var %in% vars_M1_M1, "Climate (global)",
    var %in% vars_M1,    "Climate (annual)",
    default = "Land cover (annual)"
  )]
  var_imp_long <- var_imp_long[!is.na(imp)]
  
  # one panel per model, legend only on last one
  show_legend_vec <- c(rep(FALSE, length(models) - 1), TRUE)
  
  plots_list <- Map(make_species_boxplot_panel,
                    MoreArgs = list(var_imp_long = var_imp_long),
                    model_name = models,
                    letter = letters[1:length(models)],
                    show_legend = show_legend_vec)
  
  # height
  n_vars <- sapply(models, function(m) uniqueN(var_imp_long[model == m, var]))
  
  final_plot <- wrap_plots(plots_list, ncol = 1) +
    plot_layout(heights = n_vars) +
    plot_annotation(
      caption = paste0("Variable importance for ", species_labels[[sp]])
    ) &
    theme(plot.caption = element_text(hjust = 0.5, size = 11, face = "bold"))
  
  # Save
  total_vars <- sum(n_vars)
  ggsave(filename = file.path(path_results_models, gsub("_", ".", sp), "models/M1_M2/variable_importance_M1_M2_plot.svg"),
         plot     = final_plot, width    = 10, height   = total_vars * 0.3,  
         dpi      = 300, device   = "svg")
  
  final_plot
}

# Plot and save svg for all species
plots_species <- lapply(species_names, make_species_boxplot)
names(plots_species) <- species_names

# show
plots_species[["Bombus_hortorum"]]

# Save pdf all species
for (sp in species_names) {
  ggsave(filename = paste0(path_results_models,"/all/var_imp_sp/variable_importance_", sp,".pdf"),
         plot     = plots_species[[sp]], width    = 10, height   = 8,  
         dpi      = 300, device   = "pdf")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 2b. Variable importance cumulative barplots ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Load var imp for each species and join
all_var_imp <- rbindlist(lapply(species_names, function(sp) {
  
  var_imp_table <- fread(file.path(path_results_models, gsub("_", ".", sp),
                                   "models/M1_M2/variable_importance_M1_M2.csv"))
  
  var_imp_long <- melt(var_imp_table, id.vars = c("var", "seed"),
                       measure.vars = c("HS1", "HS2", "HS_PR", "HS_multiply"),
                       variable.name = "model", value.name = "imp")
  
  var_imp_long[, category := fcase(
    var %in% vars_M1_M1, "Climate (global)",
    var %in% vars_M1,    "Climate (annual)",
    default = "Land cover (annual)"
  )]
  var_imp_long <- var_imp_long[!is.na(imp)]
  var_imp_long[, species := sp]
  
  var_imp_long
}))

# Mean per species
mean_imp <- all_var_imp[, .(mean_imp = mean(imp, na.rm = TRUE)),
                        by = .(model, var, category, species)]

#change order
mean_imp[, species := factor(species, levels = rev(species_names))]


#function to make a cumulative barplot

species_generalists <- c("Bombus_hortorum", "Bombus_pascuorum", "Bombus_terrestris")
species_localized   <- c("Bombus_alpinus", "Bombus_pyrenaeus", "Bombus_wurflenii")

make_cumulative_barplot <- function(model_name, species_subset = species_names, letter) {
  
  d <- mean_imp[model == model_name & species %in% species_subset]
  
  #order vars (décroissant, propre au sous-ensemble d'espèces)
  var_order <- d[, .(total_imp = sum(mean_imp, na.rm = TRUE)), by = var][order(total_imp), var]
  d[, var := factor(var, levels = var_order)]
  
  # Add a point for the category on the left of the barplot
  var_cat <- unique(d[, .(var, category)])
  var_cat[, x_pos := -max(d$mean_imp, na.rm = TRUE) * 0.04] #for the point
  
  # plot
  ggplot(d, aes(x = var, y = mean_imp, fill = species)) +
    geom_col(position = "stack", width = 0.7) +
    geom_point(data = var_cat, aes(x = var, y = x_pos, color = category, shape = category),
               inherit.aes = FALSE, size = 2.5) +
    coord_flip() +
    scale_x_discrete(labels = function(x) ifelse(x %in% names(var_labels), var_labels[x], x)) +
    scale_fill_manual(values = species_colors, breaks = species_names, labels = species_labels,
                      drop = TRUE) +
    scale_color_manual(values = vars_colors, breaks = names(vars_colors), name = "Variable type") +
    scale_shape_manual(values = vars_shape, breaks = names(vars_colors), name = "Variable type") +
    theme_minimal() +
    theme(
      plot.title  = element_text(face = "bold", size = 11, hjust = 0.5),
      axis.text.y = element_text(size = 9),
      plot.tag = element_text( size = 10),
      plot.tag.position = c(0, 1),
      panel.border=element_rect(color = "grey", linewidth=0.2), #with rectangles, TO REMOVE IF NOT NEEDED
    ) +
    labs(title = models_corresp[[model_name]], x = NULL, y = NULL, fill = "Species",
         tag = paste0("(", letter, ")"))
}

# Function to do the figure according to species subset
make_full_figure <- function(species_subset) {
  
  plots_list <- Map(make_cumulative_barplot,
                    model_name = models,
                    species_subset = rep(list(species_subset), length(models)),
                    letter = letters[1:length(models)])
  
  wrap_plots(plots_list, ncol = 1) +
    plot_layout(guides = "collect", heights = n_vars) +
    plot_annotation(
      caption = paste0("Cumulative mean variable importance")
    ) &
    theme(plot.caption = element_text(hjust = 0.5, size = 11, face = "bold"))
}

# Les trois figures
final_plot_all          <- make_full_figure(species_names)
# final_plot_generalists  <- make_full_figure(species_generalists)
# final_plot_localized    <- make_full_figure(species_localized)

# final_plot_all
# final_plot_generalists
# final_plot_localized

# save
n_vars <- sapply(models, function(m) uniqueN(mean_imp[model == m & species %in% species_names, var]))

total_vars <- sum(n_vars)

ggsave(filename = file.path(path_results_models, "all/variable_importance_cumulative_barplot_rect.svg"),
       plot     = final_plot_all, width    = 10, height   = total_vars * 0.3,  
       dpi      = 300,device   = "svg")
# ggsave(filename = file.path(path_results_models, "all/variable_importance_cumulative_barplot_point_generalists.svg"),
#        plot     = final_plot_generalists, width    = 10, height   = total_vars * 0.3,  
#        dpi      = 300,device   = "svg")
# ggsave(filename = file.path(path_results_models, "all/variable_importance_cumulative_barplot_point_localized.svg"),
#        plot     = final_plot_localized, width    = 10, height   = total_vars * 0.3,  
#        dpi      = 300, device   = "svg")


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 3. Descriptive table observation ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# test

# M1

results_list <- list()

for (sp in species_names){
  
  # Load model
  M1 <- get(load(file.path(path_results_models,gsub("_", ".", sp),paste0(gsub("_", ".", sp), ".M1_1.models.out"))))
  
  # Calib lines: matrix [obs,run] with TRUE/FALSE for each obs and run
  cal_M1 <- get_calib_lines(M1)
  
  # resp_var for which we want to calculate the number used in the evaluation set : Presence or absence
  resp_var <- get_formal_data(M1, "resp.var")
  
  for (run_name in colnames(cal_M1)){
    
    lines_run <- cal_M1[, run_name]  
    resp_var_run   <- resp_var[!lines_run]
    
    nb_pres <- sum(resp_var_run == 1, na.rm = TRUE)
    nb_abs  <- sum(resp_var_run == 0, na.rm = TRUE)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      species = sp,
      run = as.integer(gsub("_allData_RUN", "", run_name)),
      nb_pres = nb_pres,
      nb_abs  = nb_abs
    )
  }
}

table_obs <- rbindlist(results_list)
table_obs

# verification of the total number of presences per species
## with all pres
load("data_prep/obs/M1/dt_10000_background.RData")
sum(is.na(dt_background))

t1<- dt_background[, .(nb_pres = sum(presence == 1, na.rm = TRUE),
                          nb_abs  = sum(presence == 0, na.rm = TRUE)),
                      by = species]

## with the calib lines done above
t2<-table_obs[, .(nb_pres = sum(nb_pres, na.rm = TRUE),
                  nb_abs  = sum(nb_abs, na.rm = TRUE)),
              by = species]
# --> OK same
write.csv(table_obs, file.path(path_results, "obs/table_obs_per_species_per_run_M1.csv"), row.names = FALSE)
print(xtable(table_obs,
             table.placement = "H"),
      include.rownames = FALSE)

# M2
load("data_prep/obs/M2/dt_1000_background.RData")
sum(is.na(dt_background))

results_list <- list()

for (sp in species_names){
  
  # Load model
  M2 <- get(load(file.path(path_results_models,gsub("_", ".", sp),paste0(gsub("_", ".", sp), ".M2_1.models.out"))))
  
  # Calib lines: matrix [obs,run] with TRUE/FALSE for each obs and run
  cal_M2 <- get_calib_lines(M2)
  
  # resp_var for which we want to calculate the number used in the evaluation set : Presence or absence
  resp_var <- get_formal_data(M2, "resp.var")
  
  for (run_name in colnames(cal_M2)){
    
    lines_run <- cal_M1[, run_name]  
    resp_var_run   <- resp_var[!lines_run]
    
    nb_pres <- sum(resp_var_run == 1, na.rm = TRUE)
    nb_abs  <- sum(resp_var_run == 0, na.rm = TRUE)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      species = sp,
      run = as.integer(gsub("_allData_RUN", "", run_name)),
      nb_pres = nb_pres,
      nb_abs  = nb_abs
    )
  }
}

table_obs <- rbindlist(results_list)
table_obs

# add total per species
table_tot<- table_obs[, .(nb_pres = sum(nb_pres, na.rm = TRUE), nb_abs=sum(nb_abs)), by=species]

write.csv(table_obs, file.path(path_results, "obs/table_obs_per_species_per_run_M2.csv"), row.names = FALSE)
print(xtable(table_obs,
             table.placement = "H"),
      include.rownames = FALSE)



# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ---- 4. Plot blocks ----
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_blocks<- function (sp){
  load(paste0("data_prep/blocks/blocks_", sp, ".RData"))
  # Plot
  p <- cv_plot(scv, r=grid,  raster_colors = '#F5EDF1', max_pixels = 1e+04, label_size=2)
  
  ggsave(filename = paste0("data_prep/blocks/blocks_", sp, ".pdf"),
         plot = p, dpi = 100, width = 7, height = 7)
}
for (sp in species_names) {
  plot_blocks(sp)
}
