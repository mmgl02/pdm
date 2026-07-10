#setwd("D:/stage/3_code")

# Libraries
library(ggplot2)
library(viridis)
library(patchwork)
library(metR)


# Load parameters
parameters_file<-file.path(dirname(rstudioapi::getActiveDocumentContext()$path), "../../parameters.R") 
source(parameters_file)


# Grid
grid <- expand.grid(
  HS1 = seq(0, 1, length.out = 200),
  HS2 = seq(0, 1, length.out = 200)
)

# Calculate
grid$HS_PR         <- 1 / (1 + ((1 - grid$HS1) / grid$HS1) * ((1 - grid$HS2) / grid$HS2))
grid$HS_multiply <- sqrt(grid$HS1 * grid$HS2)
grid$HS_sum      <- (grid$HS1 + grid$HS2) / 2

# Plots
p1 <- ggplot(grid, aes(x = HS1, y = HS2, fill = HS_PR, z = HS_PR)) +
  geom_raster() +
  geom_contour(aes(z = HS_PR), breaks = seq(0, 0.9, by = 0.1), col = 'black') +
  geom_label_contour(breaks = 0.5, label.placer = label_placer_fraction(frac = 0.5), colour = "black", fill = "transparent", label.size = 0) +
  scale_fill_viridis_c(name = "f(HS1,HS2)") +
  labs(title = "PR(HS1, HS2)") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text()
  ) +
  coord_fixed()

p2<-ggplot(grid, aes(x = HS1, y = HS2, fill = HS_multiply, z=HS_multiply)) +
  geom_raster() +
  geom_contour(aes(z = HS_multiply), breaks = seq(0, 0.9, by = 0.1), col='black') +
  geom_label_contour(breaks=0.5,label.placer = label_placer_fraction(frac = 0.5) ,colour = "black", fill = "transparent", label.size = 0)+
  scale_fill_viridis_c(name = "f(HS1,HS2)") +
  labs(title = expression(bold(sqrt('HS1 * HS2')))) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text()
  ) +
  coord_fixed()

p3<-ggplot(grid, aes(x = HS1, y = HS2, fill = HS_multiply, z=HS_multiply)) +
  geom_raster() +
  geom_contour(aes(z = HS_sum), breaks = seq(0, 0.9, by = 0.1), col='black') +
  geom_label_contour(breaks=0.5,label.placer = label_placer_fraction(frac = 0.5) ,colour = "black", fill = "transparent", label.size = 0)+
  scale_fill_viridis_c(name = "f(HS1,HS2)") +
  labs(title = "(HS1+HS2) /2") +
  theme_minimal() +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank(),
    axis.ticks = element_line(color = "black"),
    axis.text = element_text()
  ) +
  coord_fixed()

# Plot the 3 graphs (appendix)
plots<-wrap_plots(p1, p2, p3, ncol=3)+ plot_layout(guides="collect")
plots
ggsave(filename = file.path(path_results_models, "all/comp_nested_methods_plus.pdf"),
       plot     = plots, dpi      = 300, height= 4, device   = "pdf")

# Plot only 2 graphs (methods section)
plots<-wrap_plots(p1, p2, ncol=2)+ plot_layout(guides="collect")
plots
ggsave(filename = file.path(path_results_models, "all/comp_nested_methods.pdf"),
       plot     = plots, dpi      = 300, height= 4, device   = "pdf")
