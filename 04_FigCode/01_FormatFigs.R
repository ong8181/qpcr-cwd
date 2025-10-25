####
#### Visualize CWD qPCR results
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2025.08.05
library(patchwork); packageVersion("patchwork") # 1.3.5, 2025.08.05
library(cols4all); packageVersion("cols4all") # 0.8, 2025.08.05
library(ggtext); packageVersion("ggtext") # 0.1.2, 2025.02.21
library(Cairo); packageVersion("Cairo") # 1.6.2, 2025.08.05
library(ggbeeswarm); packageVersion("ggbeeswarm") # 0.7.2, 2025.09.27

# Create output directory
set.seed(1234)

# ------------------------------------------------ #
# Load the figure objects
# ------------------------------------------------ #
# LOD, LOQ and Amplification Curves
gglod <- readRDS("data_robj/Fig_LOD_LOQ.obj")
l1 <- gglod[[1]]
l2 <- gglod[[4]]
l3 <- gglod[[2]]
l4 <- readRDS("data_robj/Fig_AmpCurve.obj")
lod_df <- read.csv("../01_LOD_LOQ/01_LoD-calculatorOut/Assay summary.csv")

# Other related species and samples
o1 <- readRDS("data_robj/Fig_OtherSpp.obj")

# Lantau qPCR result
r1 <- readRDS("data_robj/Fig_Survey_summer.obj")
r2 <- readRDS("data_robj/Fig_Survey_winter.obj")
# CTD results
c1 <- readRDS("data_robj/Fig_CTD_summer.obj")
c2 <- readRDS("data_robj/Fig_CTD_winter.obj")
# qPCR v.s. Environment
e1 <- readRDS("data_robj/Fig_qPCRenv_summer.obj")
e2 <- readRDS("data_robj/Fig_qPCRenv_summer_logis.obj")
e3 <- readRDS("data_robj/Fig_qPCRenv_winter.obj")
e4 <- readRDS("data_robj/Fig_qPCRenv_winter_logis.obj")


# ------------------------------------------------ #
# LOD and LOQ
# ------------------------------------------------ #
# LOD
l1$layers[[5]]$data$y <- 20
l1$layers[[6]]$data$y <- 20
l1$layers[[7]]$data$x <- 7000
l1$layers[[7]]$data$y <- 38
l1 <- l1 + labs(title = NULL, x = "Standard DNA conc. (copies/reaction)",
                y = "Ct")

# LOD plot
l2 <- l2 +
  annotate("text", y = 0.15, x= 1.8,label="LOD", color = "red", angle=90) +
  theme(legend.position = c(0.95, 0.05), legend.justification = c(1, 0)) +
  NULL

# LOQ
## Re-draw regression
l3$layers[[3]] <- NULL # regression
l3$layers[[4]] <- NULL # polygon
l3 <- l3 + geom_segment(data = lod_df, aes(x = 0.5, y = 0.35, xend = LOQ, yend = 0.35), linetype = 2, linewidth = 0.2) +
  geom_segment(data = lod_df, aes(x = LOQ, y = 0.35, xend = LOQ, yend = -0.1), linetype = 2, linewidth = 0.2) +
  stat_smooth(method = "lm", formula = y ~ poly(x,6),se=FALSE, size = 0.5) +
  geom_point(data = lod_df, aes(x = LOQ, y = 0.35), color = "red3", shape = 8, size = 2) +
  annotate("text", y = 0.1, x= 2, label="LOD", color = "red", angle=90) +
  annotate("text", y = 0.1, x= 8.5, label="LOQ", color = "black", angle=90) +
  labs(title = NULL,
       y = "Coefficent of variation",
       x = "Standard DNA conc. (copies/reaction)") +
  coord_cartesian(xlim = c(1,2000000), ylim = c(0,0.7)) +
  theme(text = element_text(family = "Arial")) + theme_bw()

# Combine LOD and LOQ plots
fig_lodloq <- (l1) /
  ((l2) + (l3)) +
  plot_layout(height = c(1.4, 1)) + 
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))

# Amplification curve
l4 <- l4 + scale_color_manual(values = c4a("carto.bold"),
                              labels = c("0.5","1.0","10","100","1,000",
                                         "10,000","100,000","1,000,000"),
                              name = "Copies/µL") +
  xlim(10,55) + labs(title = NULL, y = "\u0394Rn", x = "Cycle") +
  theme(axis.text.y = element_markdown())


# ------------------------------------------------ #
# Amplification of other related species
# ------------------------------------------------ #
std_pattern <- o1[[1]] + labs(title = NULL, x = "Cycle", y = "\u0394Rn") +
  guides(color = guide_legend(ncol = 2)) +
  theme(axis.text.y = element_markdown())
stenella <- o1[[2]] + labs(title = NULL, x = "Cycle", y = "\u0394Rn") +
  theme(axis.text.y = element_markdown(), legend.position = "top")
delphinus <- o1[[3]] + labs(title = NULL, x = "Cycle", y = "\u0394Rn") +
  theme(axis.text.y = element_markdown(), legend.position = "top")
tursiops <- o1[[4]] + labs(title = NULL, x = "Cycle", y = "\u0394Rn") +
  theme(axis.text.y = element_markdown(), legend.position = "top")
tissue <- o1[[5]] + labs(title = NULL, x = "Cycle", y = "\u0394Rn") +
  scale_color_manual(values = c("gray60", "pink2"),
                     labels = c("tissue" = "CWD tissue")) +
  theme(axis.text.y = element_markdown(), legend.position = "top")

# Combine all the figures
fig_otherspp <- (std_pattern +
                   theme(legend.position = "top", plot.tag.position = c(0,0.9)) +
                   guides(color = guide_legend(nrow = 2))) /
  ((tissue) + (stenella)) /
  ((delphinus) + (tursiops)) +
  plot_layout(height = c(1.2, 1, 1)) + 
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))


# ------------------------------------------------ #
# Lantau results
# ------------------------------------------------ #
# Main qPCR results
## Surface and bottom combined
fig_qpcr_qual <- 
    ((r1[[3]] + ggtitle("2024-07, Summer survey") + theme(legend.position = "none"))) +
    ((r2[[3]] + ggtitle("2025-01, Winter survey") + theme(legend.position = "right"))) +
    plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))
fig_qpcr_conc <-
  ((r1[[1]] + ggtitle("2024-07, Surface layer") + theme(legend.position = "none")) +
     r1[[2]] + ggtitle("2024-07, Bottom layer") ) /
  ((r2[[1]] + ggtitle("2025-01, Surface layer") + theme(legend.position = "none")) +
     r2[[2]] + ggtitle("2025-01, Bottom layer") ) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))

# CTD sensor
fig_ctd <- (c1[[1]] + ggtitle("2024-07, Summer survey") + c1[[2]] + c1[[3]]) /
  (c2[[1]] + ggtitle("2025-01, Winter survey") + c2[[2]] + c2[[3]])  +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))

# Environment v.s. qPCR
## Customize plots
e2_1 <- e2[[1]]; e2_1$layers[[2]] <- NULL
e2_1 <- e2_1 + geom_beeswarm(aes(x = sur_bot, y = as.numeric(detection1)), cex = 2.2)
e4_1 <- e4[[1]]; e4_1$layers[[2]] <- NULL
e4_1 <- e4_1 + geom_beeswarm(aes(x = sur_bot, y = as.numeric(detection1)), cex = 2.2)
# Adjust text positions
e2_2 <- e2[[2]]; e2_2$layers[[3]]$data$x <- 27.2
e2_3 <- e2[[3]]; e2_3$layers[[3]]$data$x <- 22
e2_4 <- e2[[4]]; e2_4$layers[[3]]$data$x <- 4.5
e4_2 <- e4[[2]]; e4_2$layers[[3]]$data$x <- 17.25; e4_2$layers[[3]]$data$y <- 0.375
e4_3 <- e4[[3]]; e4_3$layers[[3]]$data$x <- 31.625
e4_4 <- e4[[4]]; e4_4$layers[[3]]$data$y <- 0.3
# Linetype
e2_2$layers[[2]]$aes_params$linetype <- "dashed"
e2_3$layers[[2]]$aes_params$linetype <- "dashed"
e2_4$layers[[2]]$aes_params$linetype <- "dashed"
e4_2$layers[[2]]$aes_params$linetype <- "dashed"
e4_3$layers[[2]]$aes_params$linetype <- "dashed"
e4_4$layers[[2]]$aes_params$linetype <- "dashed"

fig_env_summer <-
  (e2_1 + e2[[2]] + scale_color_manual(values = c("royalblue","red3"), name = NULL)) /
  (e2[[3]] + scale_color_manual(values = c("royalblue","red3")) + theme(legend.position = "none") +
     e2[[4]] + scale_color_manual(values = c("royalblue","red3"), name = NULL)) +
  plot_annotation(tag_levels = "a") + theme(plot.tag = element_text(face = "bold"))
fig_env_winter <-
  (e4_1 + e4[[2]] + scale_color_manual(values = c("royalblue","red3"), name = NULL)) /
  (e4[[3]] + scale_color_manual(values = c("royalblue","red3")) + theme(legend.position = "none") +
     e4[[4]] + scale_color_manual(values = c("royalblue","red3"), name = NULL)) +
  plot_annotation(tag_levels = "a") + theme(plot.tag = element_text(face = "bold"))


# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
# Save figures
## Main figures
ggsave("formatted_figs/Fig02_StdCurve.pdf", fig_lodloq, height = 8, width = 9)
ggsave("formatted_figs/Fig03_SurveyqPCR.pdf", fig_qpcr_qual, height = 6, width = 12)

## Supporting figures
CairoPDF("formatted_figs/FigS01_RelatedSppAmp.pdf", height = 14, width = 12); fig_otherspp; dev.off()
CairoPDF("formatted_figs/FigS02_SurveyCTD.pdf", height = 13, width = 11); fig_ctd; dev.off()
ggsave("formatted_figs/FigS03_SurveyqPCR_conc.pdf", fig_qpcr_conc, height = 10, width = 12)
CairoPDF("formatted_figs/FigS04_SurveyEnvLogis1.pdf", height = 8, width = 10); fig_env_summer; dev.off()
CairoPDF("formatted_figs/FigS05_SurveyEnvLogis2.pdf", height = 8, width = 10); fig_env_winter; dev.off()

# Save sessioninfo
macam::save_session_info()
