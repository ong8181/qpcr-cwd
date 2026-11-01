####
#### Testing amplification of closely related species
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0
library(patchwork); packageVersion("patchwork") # 1.3.2, 2026.10.22
library(cols4all); packageVersion("cols4all") # 0.10, 2026.10.22

# Create output
(outdir <- macam::outdir_create())

# Read data
d1 <- read.csv("data_otherspp/data_R0.csv") # Original (tissue sample included)
d2 <- read.csv("data_otherspp/data_R1.csv") # Revision (Related species included)
d_amp1 <- read.csv("data_otherspp/data_amplification_R0.csv") # Original
d_amp2 <- read.csv("data_otherspp/data_amplification_R1.csv") # Revision
#head(d)

# Extract standard samples
d_std1 <- d1 %>% filter(task == "STANDARD")
d_std2 <- d2 %>% filter(task == "STANDARD")
d_tis <- d1 %>% filter(sample_name == "Tissue")
# Conc = 100 copies/µL
d_stn1 <- d2 %>% filter(sample_name == "Stenella" & quantity == 100)
d_del1 <- d2 %>% filter(sample_name == "Delphinus" & quantity == 100)
d_tur1 <- d2 %>% filter(sample_name == "Tursiops" & quantity == 100)
# Conc = 1000 copies/µL
d_stn2 <- d2 %>% filter(sample_name == "Stenella" & quantity == 1000)
d_del2 <- d2 %>% filter(sample_name == "Delphinus" & quantity == 1000)
d_tur2 <- d2 %>% filter(sample_name == "Tursiops" & quantity == 1000)


# ------------------------------------------------------- #
# Calculate DNA amount
# ------------------------------------------------------- #
# DNA copy * 2 ul / (Avogadro constant) * (Molecular Weight) * fg
d_std1$pg <- d_std1$quantity * 2 / (6.022 * 10^23) * (152 * 660) * 10^12
d_std1$fg <- d_std1$pg*1000 %>% round(digits = 3)
d_std2$pg <- d_std2$quantity * 2 / (6.022 * 10^23) * (152 * 660) * 10^12
d_std2$fg <- d_std2$pg*1000 %>% round(digits = 3)


# ------------------------------------------------------- #
# Extract standard sample data
# ------------------------------------------------------- #
# Original test
d_amp3 <- d_amp1 %>% filter(d_amp1$Well %in% d_std1$Well)
d_amp3$conc <- d_std1[match(d_amp1$Well, d_std1$Well), "quantity"] %>% na.omit %>% as.numeric
d_amp3$pg <- d_std1[match(d_amp1$Well, d_std1$Well), "pg"] %>% na.omit %>% as.numeric
d_amp3$conc <- factor(d_amp3$conc, levels = c("1","10","100","1000","10000","1e+05","1e+06","1e+07","1e+08","1e+09"))
# Test for revision
d_amp4 <- d_amp2 %>% filter(d_amp2$Well %in% d_std2$Well)
d_amp4$conc <- d_std2[match(d_amp2$Well, d_std2$Well), "quantity"] %>% na.omit %>% as.numeric
d_amp4$pg <- d_std2[match(d_amp2$Well, d_std2$Well), "pg"] %>% na.omit %>% as.numeric
d_amp4$conc <- factor(d_amp4$conc, levels = c("0","1","10","100","1000","10000"))


# ------------------------------------------------------- #
# Visualize
# ------------------------------------------------------- #
g1 <- d_amp3 %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = factor(conc))) +
  geom_line(data = d_amp4, aes(x = Cycle, y = Delta.Rn, group = Well), color = "black") +
  geom_line() + scale_color_manual(values = rev(c4a("carto.bold"))[2:11],
                                   name = "copies/µL",
                                   labels = c("1","10","100","1,000","10,000",
                                              "100,000","1,000,000","10,000,000",
                                              "100,000,000","1,000,000,000")) +
  xlim(0,55) + ylim(0,1.4) + ggtitle("Delta Rn (CWD Standard DNA)") +
  NULL
g2 <- d_amp3 %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = factor(conc))) +
  geom_line() + scale_color_manual(values = rev(c4a("carto.bold"))[2:11],
                                   name = "copies/µL",
                                   labels = c("1","10","100","1,000","10,000",
                                              "100,000","1,000,000","10,000,000",
                                              "100,000,000","1,000,000,000")) +
  ylim(0,1.4) + ggtitle("CWD Standard DNA for tissue DNA") +
  NULL
g3 <- d_amp4 %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = factor(conc))) +
  geom_line() + scale_color_manual(values = rev(c4a("carto.bold"))[1:6],
                                   name = "copies/µL",
                                   labels = c("0", "1","10","100","1,000","10,000")) +
  xlim(0,55) + ylim(0,1.4) + ggtitle("CWD Standard DNA for closely-related species") +
  NULL

(g1/g2/g3)


# ------------------------------------------------------- #
# Closely related species
# ------------------------------------------------------- #
# 100 copies/µL samples
# Choose target wells
d_amp_stn1 <- d_amp2 %>% filter((d_amp2$Well %in% d_std2$Well) | (d_amp2$Well %in% d_stn1$Well))
d_amp_del1 <- d_amp2 %>% filter((d_amp2$Well %in% d_std2$Well) | (d_amp2$Well %in% d_del1$Well))
d_amp_tur1 <- d_amp2 %>% filter((d_amp2$Well %in% d_std2$Well) | (d_amp2$Well %in% d_tur1$Well))
# Assign sample names
d_amp_stn1$sample <- d_amp_del1$sample <- d_amp_tur1$sample <- "standard"
d_amp_stn1$sample[d_amp_stn1$Well %in% d_stn1$Well] <- "Stenella"
d_amp_del1$sample[d_amp_del1$Well %in% d_del1$Well] <- "Delphinus"
d_amp_tur1$sample[d_amp_tur1$Well %in% d_tur1$Well] <- "Tursiops"
d_amp_stn1$sample <- factor(d_amp_stn1$sample, levels = c("standard", "Stenella"))
d_amp_del1$sample <- factor(d_amp_del1$sample, levels = c("standard", "Delphinus"))
d_amp_tur1$sample <- factor(d_amp_tur1$sample, levels = c("standard", "Tursiops"))

# 1000 copies/µL samples
# Choose target wells
d_amp_stn2 <- d_amp2 %>% filter((d_amp2$Well %in% d_std2$Well) | (d_amp2$Well %in% d_stn2$Well))
d_amp_del2 <- d_amp2 %>% filter((d_amp2$Well %in% d_std2$Well) | (d_amp2$Well %in% d_del2$Well))
d_amp_tur2 <- d_amp2 %>% filter((d_amp2$Well %in% d_std2$Well) | (d_amp2$Well %in% d_tur2$Well))
# Assign sample names
d_amp_stn2$sample <- d_amp_del2$sample <- d_amp_tur2$sample <- "standard"
d_amp_stn2$sample[d_amp_stn2$Well %in% d_stn2$Well] <- "Stenella"
d_amp_del2$sample[d_amp_del2$Well %in% d_del2$Well] <- "Delphinus"
d_amp_tur2$sample[d_amp_tur2$Well %in% d_tur2$Well] <- "Tursiops"
d_amp_stn2$sample <- factor(d_amp_stn2$sample, levels = c("standard", "Stenella"))
d_amp_del2$sample <- factor(d_amp_del2$sample, levels = c("standard", "Delphinus"))
d_amp_tur2$sample <- factor(d_amp_tur2$sample, levels = c("standard", "Tursiops"))

# 100 copies/µL samples
# Stenella
c1 <- d_amp_stn1 %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = sample)) +
  geom_line() + scale_color_manual(values = c("gray60", "red2")) +
  xlim(0,55) + ylim(0,1.4) +
  ggtitle("Stenella (100 copies/µL)")
# Delphinus
c2 <- d_amp_del1 %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = sample)) +
  geom_line() + scale_color_manual(values = c("gray60", "red2")) +
  xlim(0,55) + ylim(0,1.4) +
  ggtitle("Delphinus (100 copies/µL)")
# Tursiops
c3 <- d_amp_tur1 %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = sample)) +
  geom_line() + scale_color_manual(values = c("gray60", "red2")) +
  xlim(0,55) + ylim(0,1.4) +
  ggtitle("Tursiops (100 copies/µL)")

# 1000 copies/µL samples
# Stenella
c4 <- d_amp_stn2 %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = sample)) +
  geom_line() + scale_color_manual(values = c("gray60", "red2")) +
  xlim(0,55) + ylim(0,1.4) +
  ggtitle("Stenella (1000 copies/µL)")
# Delphinus
c5 <- d_amp_del2 %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = sample)) +
  geom_line() + scale_color_manual(values = c("gray60", "red2")) +
  xlim(0,55) + ylim(0,1.4) +
  ggtitle("Delphinus (1000 copies/µL)")
# Tursiops
c6 <- d_amp_tur2 %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = sample)) +
  geom_line() + scale_color_manual(values = c("gray60", "red2")) +
  xlim(0,55) + ylim(0,1.4) +
  ggtitle("Tursiops (1000 copies/µL)")


# ------------------------------------------------------- #
# Tissue samples and natural samples
# ------------------------------------------------------- #
# Choose target wells
d_amp_tis <- d_amp1 %>% filter((d_amp1$Well %in% d_std1$Well) | (d_amp1$Well %in% d_tis$Well))
d_amp_tis$sample <- "standard"
d_amp_tis$sample[d_amp_tis$Well %in% d_tis$Well] <- "tissue"

# Tissue + Natural
t1 <- d_amp_tis %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = sample)) +
  geom_line() + scale_color_manual(values = c("gray60", "pink")) +
  ylim(0,1.4) +
  ggtitle("Tissue DNA")


# ------------------------------------------------------- #
# Combine all plots
# ------------------------------------------------------- #
# Closely-related species
comb1 <- (g3 | g3) + plot_layout(guides = "collect")
comb2 <- (c1 | c4) + plot_layout(guides = "collect")
comb3 <- (c2 | c5) + plot_layout(guides = "collect")
comb4 <- (c3 | c6) + plot_layout(guides = "collect")
right_panel1 <- (comb1 / comb2 / comb3 / comb4)
right_panel2 <- (comb2 / comb3 / comb4)
# Tissue samples
left_panel1 <- (g2 / t1)
# Combine
all_panel1 <- ((left_panel1 / plot_spacer()/ plot_spacer()) | right_panel1) + plot_layout(widths = c(1,2))
all_panel2 <- ((g2 / t1 / g3) | right_panel2) + plot_layout(widths = c(1,2))


# ------------------------------------------------------- #
# Save results
# ------------------------------------------------------- #
# Save PDFs
ggsave(file.path(outdir, "/AmplificationPlot2.pdf"), plot = all_panel2, width = 16, height = 8)

# Save R objects
saveRDS(list(g1, g2, g3, c1, c2, c3, c4, c5, c6, t1), "../04_FigCode/data_robj/Fig_OtherSpp.obj")

# Save session information
macam::save_session_info()
