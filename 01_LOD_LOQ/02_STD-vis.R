####
#### qPCR standard curve
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0
library(patchwork); packageVersion("patchwork") # 1.2.0
library(cols4all); packageVersion("cols4all") # 0.8, 2025.08.05
library(ggsci); packageVersion("ggsci") # 3.0.0

# Create output
outdir <- macam::outdir_create()

# Read data
d_amp <- read.csv("data_lod/data_amp.csv")
d_amp$quantity <- as.character(d_amp$quantity)


# ------------------------------------------------------- #
# Visualize
# ------------------------------------------------------- #
(g1 <- d_amp %>% #filter(plate == "P2") %>% 
    ggplot(aes(x = cycle, y = delta_rn, group = id, color = factor(quantity))) +
    geom_line() + ggtitle("Delta Rn") + scale_color_igv(name = "Copies/µL") +
    facet_wrap(~plate, labeller = as_labeller(c("P1" = "Conc. 0.5-100 copies/µL x 24 replicates",
                                                "P2" = "Conc. 1,000-1,000,000 copies/µL x 24 replicates"))) +
    theme(legend.position = "bottom"))

g2 <- d_amp %>% 
  ggplot(aes(x = Cycle, y = Rn, group = Well, color = factor(conc))) +
  geom_line() + scale_color_igv(name = "STD conc.") + ggtitle("Rn")

g3 <- d_amp %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = factor(fg))) +
  geom_line() + scale_color_igv(name = "DNA (fg)") + ggtitle("Delta Rn")

# Combine and save
saveRDS(g1, "../04_FigCode/data_robj/Fig_AmpCurve.obj")
ggsave(paste0(outdir, "/AmplificationPlot.pdf"), plot = g1, width = 8, height = 6)

macam::save_session_info()

