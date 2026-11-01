####
#### Edit LOQ figures
####

# Load packages
library(tidyverse); packageVersion("tidyverse") # 2.0.0
library(patchwork); packageVersion("patchwork") # 1.3.2
library(macam); packageVersion("macam") # 0.2.6

# Load figure objects
g1 <- readRDS("m_all_LoQ-calculator_DecayOut/DecayPlot.obj")
g2 <- readRDS("m_all_LoQ-calculator_LinearOut/DecayPlot.obj")
g3 <- readRDS("m_all_LoQ-calculator_P2Out/DecayPlot.obj")
g4 <- readRDS("m_all_LoQ-calculator_P3Out/DecayPlot.obj")
g5 <- readRDS("m_all_LoQ-calculator_P4Out/DecayPlot.obj")
g6 <- readRDS("m_all_LoQ-calculator_P5Out/DecayPlot.obj")
g7 <- readRDS("m_all_LoQ-calculator_P6Out/DecayPlot.obj")

(outdir <- outdir_create())


# ----------------------------------------------------------------------------- #
# Extract coordinates from the grey polygon layer and remove it
# ----------------------------------------------------------------------------- #
edit_fig <- function(ggobj, title, loq_x = 8.5) {
  poly_idx <- which(sapply(ggobj$layers, function(x) inherits(x$geom, "GeomPolygon")))
  
  # Extract the coordinate data (PDAT) stored inside the polygon layer before deleting it
  pdat <- ggobj$layers[[poly_idx]]$data
  x_min <- min(pdat$x, na.rm = TRUE) # Left edge (lowest standard concentration)
  x_loq <- max(pdat$x, na.rm = TRUE) # Right edge (LOQ concentration)
  y_min <- min(pdat$y, na.rm = TRUE) # Bottom edge
  y_loq <- max(pdat$y, na.rm = TRUE) # Top edge (CV threshold at LOQ)
  
  # Remove the grey polygon layer from the plot
  ggobj$layers[[poly_idx]] <- NULL
  
  # Add the dashed segments to the plot
  seg_df <- data.frame(x    = c(x_min - 0.5, x_loq),
                       xend = c(x_loq, x_loq),
                       y    = c(y_loq, -0.1),
                       yend = c(y_loq, y_loq))
  #l3$layers[[4]]$aes_params$linewidth <- 0.5
  gg_modified <- ggobj +
    annotate("text", y = 0.1, x= 2, label="LOD", color = "red", angle=90) +
    annotate("text", y = 0.1, x= loq_x, label="LOQ", color = "black", angle=90) +
    geom_segment(data = seg_df,
                 aes(x = x, xend = xend, y = y, yend = yend),
                 linetype = "dotted",
                 color = "black",
                 linewidth = 0.5,
                 inherit.aes = FALSE) +
    labs(title = title) +
    coord_cartesian(xlim = c(0.8,2*10^6),
                    ylim = c(0,0.6))
  
  # Display the updated plot
  return(gg_modified)
}


# Compile figures
g_all <- edit_fig(g1, "Exponential decay model (LOQ = 14 copies/reaction)", 9) +
  edit_fig(g2, "Linear model (LOQ = 50 copies/reaction)", 30) +
  edit_fig(g3, "2nd polynomial (LOQ = 24 copies/reaction)", 14) + 
  edit_fig(g4, "3rd polynomial (LOQ = 11 copies/reaction)", 6.5) +
  edit_fig(g5, "4th polynomial (LOQ = 8 copies/reaction)", 5) +
  edit_fig(g6, "5th polynomial (LOQ = 15 copies/reaction)", 9) +
  edit_fig(g7, "6th polynomial (LOQ = 13 copies/reaction)", 8.5) +
  plot_spacer() + plot_layout(ncol = 2)

# Save figures
ggsave(file.path(outdir, "LOQ_all_models.pdf"), g_all, width = 10, height = 10)
saveRDS(list(edit_fig(g1, "Exponential decay model (LOQ = 14 copies/reaction)", 9),
             edit_fig(g2, "Linear model (LOQ = 50 copies/reaction)", 30),
             edit_fig(g3, "2nd polynomial (LOQ = 24 copies/reaction)", 14), 
             edit_fig(g4, "3rd polynomial (LOQ = 11 copies/reaction)", 6.5),
             edit_fig(g5, "4th polynomial (LOQ = 8 copies/reaction)", 5),
             edit_fig(g6, "5th polynomial (LOQ = 15 copies/reaction)", 9),
             edit_fig(g7, "6th polynomial (LOQ = 13 copies/reaction)", 8)), "../../04_FigCode/data_robj/Fig_LOQ_all_models.obj")
