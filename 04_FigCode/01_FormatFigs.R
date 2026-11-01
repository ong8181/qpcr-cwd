####
#### Visualize CWD qPCR results
####

# Load libraries
library(tidyverse); packageVersion("tidyverse") # 2.0.0, 2025.08.05
library(patchwork); packageVersion("patchwork") # 1.3.2, 2026.10.25
library(cols4all); packageVersion("cols4all") # 0.10, 2026.10.22
library(ggtext); packageVersion("ggtext") # 0.2.0, 2026.10.22
library(Cairo); packageVersion("Cairo") # 1.7.0, 2026.10.22
library(ggbeeswarm); packageVersion("ggbeeswarm") # 0.7.3, 2026.10.22

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
lod_df <- read.csv("../01_LOD_LOQ/01_LoDLoQ-calculatorOut/Assay summary.csv")

# Other LOQ models
all_loq <- readRDS("data_robj/Fig_LOQ_all_models.obj")
q1 <- all_loq[[1]]
q2 <- all_loq[[2]]
q3 <- all_loq[[3]]
q4 <- all_loq[[4]]
q5 <- all_loq[[5]]
q6 <- all_loq[[6]]
q7 <- all_loq[[7]]

# Other related species and samples
o1 <- readRDS("data_robj/Fig_OtherSpp.obj")

# Lantau qPCR result
r1 <- readRDS("data_robj/Fig_Survey_summer.obj")
#r2 <- readRDS("data_robj/Fig_Survey_winter.obj")
# CTD results
c1 <- readRDS("data_robj/Fig_CTD_summer.obj")
#c2 <- readRDS("data_robj/Fig_CTD_winter.obj")
# qPCR v.s. Environment
e1 <- readRDS("data_robj/Fig_qPCRenv_summer.obj")
e2 <- readRDS("data_robj/Fig_qPCRenv_summer_logis.obj")
#e3 <- readRDS("data_robj/Fig_qPCRenv_winter.obj")
#e4 <- readRDS("data_robj/Fig_qPCRenv_winter_logis.obj")


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
  labs(title = "Two-parameter Weibull II\n(LOD = 3.226 copies/reaction for 1 rep)") +
  NULL

# LOQ
## Re-draw regression
l3$layers[[3]] <- NULL # delete polygon
l3$layers[[4]] <- NULL # to re-draw regression line
l3 <- l3 + geom_segment(data = lod_df, aes(x = 0.5, y = 0.35, xend = LOQ, yend = 0.35), linetype = 2, linewidth = 0.2) +
  geom_segment(data = lod_df, aes(x = LOQ, y = 0.35, xend = LOQ, yend = -0.1), linetype = 2, linewidth = 0.2) +
  stat_smooth(method = "nls", formula = y ~ SSasymp(x, Asym, R0, lrc), se = FALSE, linewidth = 0.35) +
  geom_point(data = lod_df, aes(x = LOQ, y = 0.35), color = "red3", shape = 8, size = 2) +
  annotate("text", y = 0.1, x= 2, label="LOD", color = "red", angle=90) +
  annotate("text", y = 0.1, x= 8.5, label="LOQ", color = "black", angle=90) +
  labs(title = "Exponential decay model\n(LOQ = 14 copies/reaction)", #"Coefficient of variation for copy number estimates",
       y = "CV for copy number estimates",
       x = "Standard DNA conc. (copies/reaction)") +
  coord_cartesian(xlim = c(1,2000000), ylim = c(0,0.7)) +
  theme(text = element_text(family = "Arial")) + theme_bw()
l3$layers[[3]]$aes_params$linewidth <- 0.35

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

# LOQ models
q1$layers$stat_smooth$aes_params$linewidth <- 0.3; q1 <- q1 + labs(y = "CV for copy number estimates")
q2$layers$stat_smooth$aes_params$linewidth <- 0.3; q2 <- q2 + labs(y = "CV for copy number estimates")
q3$layers$stat_smooth$aes_params$linewidth <- 0.3; q3 <- q3 + labs(y = "CV for copy number estimates")
q4$layers$stat_smooth$aes_params$linewidth <- 0.3; q4 <- q4 + labs(y = "CV for copy number estimates")
q5$layers$stat_smooth$aes_params$linewidth <- 0.3; q5 <- q5 + labs(y = "CV for copy number estimates")
q6$layers$stat_smooth$aes_params$linewidth <- 0.3; q6 <- q6 + labs(y = "CV for copy number estimates")
q7$layers$stat_smooth$aes_params$linewidth <- 0.3; q7 <- q7 + labs(y = "CV for copy number estimates")
fig_loq_models <- q1 + q2 + q3 + q4 + q5 + q6 + q7 + plot_spacer() + 
  plot_layout(ncol = 2) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))

# ------------------------------------------------ #
# Amplification of other related species
# ------------------------------------------------ #
# Tissue sample amplification
std_pattern_tissue <- o1[[2]] + labs(title = NULL, x = "Cycle", y = "\u0394Rn") +
  guides(color = guide_legend(ncol = 2)) +
  labs(title = "Amplification curve for standard DNAs") +
  theme_bw() +
  theme(axis.text.y = element_markdown()) +
  NULL
tissue <- o1[[10]] + labs(title = NULL, x = "Cycle", y = "\u0394Rn") +
  scale_color_manual(values = c("gray50", "pink2"), name = NULL,
                     labels = c("standard" = "Standard", "tissue" = "CWD tissue")) +
  labs(title = "Tissue DNA amplification") +
  theme_bw() +
  NULL
# Closely-related species DNA amplification
std_pattern_otherspp <- o1[[3]] + labs(title = NULL, x = "Cycle", y = "\u0394Rn") +
  labs(title = "Amplification curve for standard DNAs") +
  theme_bw() +
  theme(axis.text.y = element_markdown())
stenella10_2 <- o1[[4]] + labs(title = NULL, x = "Cycle", y = "\u0394Rn") +
  labs(title = "_Stenella_ (100 copies/µL)") +
  scale_color_manual(values = c("gray50","red3"), name = NULL,
                     labels = c("standard" = "Standard", "Stenella" = "_Stenella_")) +
  theme_bw() +
  theme(axis.text.y = element_markdown(), plot.title = element_markdown(),
        legend.text = element_markdown(), legend.position = "none")
delphinus10_2 <- o1[[5]] + labs(title = NULL, x = "Cycle", y = "\u0394Rn") +
  labs(title = "_Delphinus_ (100 copies/µL)") +
  scale_color_manual(values = c("gray50","red3"), name = NULL,
                     labels = c("standard" = "Standard", "Delphinus" = "_Delphinus_")) +
  theme_bw() +
  theme(axis.text.y = element_markdown(), plot.title = element_markdown(),
        legend.text = element_markdown(), legend.position = "none")
tursiops10_2 <- o1[[6]] + labs(title = NULL, x = "Cycle", y = "\u0394Rn") +
  labs(title = "_Tursiops_ (100 copies/µL)") +
  scale_color_manual(values = c("gray50","red3"), name = NULL,
                     labels = c("standard" = "Standard", "Tursiops" = "_Tursiops_")) +
  theme_bw() +
  theme(axis.text.y = element_markdown(), plot.title = element_markdown(),
        legend.text = element_markdown(), legend.position = "none")
stenella10_3 <- o1[[7]] + labs(title = NULL, x = "Cycle", y = "\u0394Rn") +
  labs(title = "_Stenella_ (1,000 copies/µL)") +
  scale_color_manual(values = c("gray50","red3"), name = NULL,
                     labels = c("standard" = "Standard", "Stenella" = "_Stenella_")) +
  theme_bw() +
  theme(axis.text.y = element_markdown(), plot.title = element_markdown(),
        legend.text = element_markdown(), legend.position = "right")
delphinus10_3 <- o1[[8]] + labs(title = NULL, x = "Cycle", y = "\u0394Rn") +
  labs(title = "_Delphinus_ (1,000 copies/µL)") +
  scale_color_manual(values = c("gray50","red3"), name = NULL,
                     labels = c("standard" = "Standard", "Delphinus" = "_Delphinus_")) +
  theme_bw() +
  theme(axis.text.y = element_markdown(), plot.title = element_markdown(),
        legend.text = element_markdown(), legend.position = "right")
tursiops10_3 <- o1[[9]] + labs(title = NULL, x = "Cycle", y = "\u0394Rn") +
  labs(title = "_Tursiops_ (1,000 copies/µL)") +
  scale_color_manual(values = c("gray50","red3"), name = NULL,
                     labels = c("standard" = "Standard", "Tursiops" = "_Tursiops_")) +
  theme_bw() +
  theme(axis.text.y = element_markdown(), plot.title = element_markdown(),
        legend.text = element_markdown(), legend.position = "right")

# Combine all the figures
fig_tissue <- std_pattern_tissue / tissue +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))

fig_otherspp <- (std_pattern_otherspp + plot_spacer() + plot_layout(widths = c(1,0.5))) /
  (stenella10_2 + stenella10_3) /
  (delphinus10_2 + delphinus10_3) /
  (tursiops10_2 + tursiops10_3) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))


# ------------------------------------------------ #
# Lantau results
# ------------------------------------------------ #
# Main qPCR results
## Surface and bottom combined
fig_qpcr_qual <- 
  ((r1[[3]]) + ggtitle(NULL) + theme(legend.position = "right")) +
  NULL

fig_qpcr_conc <-
  ((r1[[1]] + ggtitle("2024-07, Surface layer") + theme(legend.position = "right")) +
     r1[[2]] + ggtitle("2024-07, Bottom layer") ) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))

# CTD sensor
fig_ctd <- (c1[[1]] + ggtitle("2024-07, Summer survey") + theme_bw()
            + c1[[2]] + theme_bw() + c1[[3]] + theme_bw()) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold")) +
  NULL

# Environment v.s. qPCR
## Customize plots
e2_1 <- e2[[1]]; e2_1$layers[[2]] <- NULL
e2_1 <- e2_1 + geom_beeswarm(aes(x = sur_bot, y = as.numeric(detection1)), cex = 2.2) + theme_bw()
# Adjust text positions
e2_2 <- e2[[2]]; e2_2$layers[[3]]$data$x <- 27.2
e2_3 <- e2[[3]]; e2_3$layers[[3]]$data$x <- 22
e2_4 <- e2[[4]]; e2_4$layers[[3]]$data$x <- 4.5
# Change theme
e2_2 <- e2_2 + scale_color_manual(values = c("royalblue","red3"), name = NULL) + theme_bw()
e2_3 <- e2_3 + scale_color_manual(values = c("royalblue","red3"))  + theme_bw() + theme(legend.position = "none")
e2_4 <- e2_4 + scale_color_manual(values = c("royalblue","red3"), name = NULL) + theme_bw()
# Linetype
e2_2$layers[[2]]$aes_params$linetype <- "dashed"
e2_3$layers[[2]]$aes_params$linetype <- "dashed"
e2_4$layers[[2]]$aes_params$linetype <- "dashed"

fig_env_summer <- ((e2_1 + e2_2) / (e2_3 + e2_4)) +
  plot_annotation(tag_levels = "a") & theme(plot.tag = element_text(face = "bold"))


# ------------------------------------------------ #
# Save results
# ------------------------------------------------ #
# Save figures
## Main figures
ggsave("formatted_figs/Fig02_RelatedSppAmp.pdf", fig_otherspp, device = cairo_pdf, height = 11, width = 11)
ggsave("formatted_figs/Fig03_StdCurve.pdf", fig_lodloq, height = 8, width = 9)
ggsave("formatted_figs/Fig04_SurveyqPCR.pdf", fig_qpcr_qual, height = 6, width = 7)

## Supporting figures
ggsave("formatted_figs/FigS01_TissueAmp.pdf", fig_tissue, device = cairo_pdf, height = 8, width = 10)
ggsave("formatted_figs/FigS02_LOQAllModels.pdf", fig_loq_models, device = cairo_pdf, height = 13, width = 12)
ggsave("formatted_figs/FigS03_SurveyCTD.pdf", fig_ctd, height = 7, width = 11)
ggsave("formatted_figs/FigS04_SurveyqPCR_conc.pdf", fig_qpcr_conc, height = 8, width = 13)
ggsave("formatted_figs/FigS05_SurveyqPCRxEnv.pdf", fig_env_summer, device = cairo_pdf, height = 8, width = 9)

# Save sessioninfo
macam::save_session_info()

