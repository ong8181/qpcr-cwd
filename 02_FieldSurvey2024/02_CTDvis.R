####
#### Lantau CWD qPCR - CTD data -
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0
library(patchwork); packageVersion("patchwork") # 1.2.0
library(ggrepel); packageVersion("ggrepel") # 0.9.5
library(ggsci); packageVersion("ggsci") # 3.0.0
library(cols4all); packageVersion("cols4all") # 0.7.1
library(cowplot); packageVersion("cowplot") # 1.1.3
library(oce); packageVersion("oce") # 1.8.3
library(Cairo); packageVersion("Cairo") # 1.6.2
theme_set(theme_bw())

# Create output directory
dir.create("02_CTDvisOut")

# Load data
d <- read.csv("data/data_ctd.csv")

# Calculate depth
d <- d %>% mutate(depth = swDepth(pressure = d$pressure_dbar, latitude = 22.2))

# Remove erroneous data
## Surface data
d <- d %>% filter(depth > 0.2)
## Get labels for the locations
d_label <- d %>% group_by(location) %>% slice(which.max(depth))


# ------------------------------------------ #
# Visualize depth profiles of environmental parameters
# ------------------------------------------ #
g1 <- d %>% ggplot(aes(x = temperature, y = depth, color = location)) +
  geom_path() +
  geom_text_repel(data = d_label, aes(y = depth, x = temperature, label = location), color = "black") +
  scale_color_manual(values = c4a("tol.rainbow")) +
  scale_y_reverse() + 
  xlab(expression(paste("Temperature (", degree, "C)"))) +
  ylab("Depth (m)") +
  theme(legend.position = "none") +
  NULL

g2 <- d %>% ggplot(aes(x = salinity_permil, y = depth, color = location)) +
  geom_path() +
  geom_text_repel(data = d_label, aes(y = depth, x = salinity_permil, label = location), color = "black") +
  scale_color_manual(values = c4a("tol.rainbow")) +
  scale_y_reverse() + 
  xlab(expression(paste("Salinity (\u2030)"))) +
  ylab("Depth (m)") +
  theme(legend.position = "none") +
  NULL

g3 <- d %>% ggplot(aes(x = DO_mg_L, y = depth, color = location)) +
  geom_path() +
  geom_text_repel(data = d_label, aes(y = depth, x = DO_mg_L, label = location), color = "black") +
  scale_color_manual(values = c4a("tol.rainbow")) +
  scale_y_reverse() + 
  xlab("DO (mg/L)") +
  ylab("Depth (m)") +
  NULL


# ------------------------------------------ #
# Save figures
# ------------------------------------------ #
## cairo_pdf() for unicode
CairoPDF("02_CTDvisOut/EnvProfile.pdf", width = 12, height = 6)
plot_grid(g1, g2, g3, nrow = 1, rel_widths = c(1,1,1.15), align = "nv", axis = "btlr")
dev.off()

CairoPDF("02_CTDvisOut/EnvProfile_w_points.pdf", width = 12, height = 6)
plot_grid(g1 + geom_point(size = 1),
          g2 + geom_point(size = 1),
          g3 + geom_point(size = 1),
          nrow = 1, rel_widths = c(1,1,1.15), align = "nv", axis = "btlr")
dev.off()

saveRDS(list(g1, g2, g3), "../04_FigCode/data_robj/Fig_CTD_summer.obj")

# Save session information
macam::save_session_info()