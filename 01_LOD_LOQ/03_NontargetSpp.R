####
#### Testing amplification of closely related species
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0
library(patchwork); packageVersion("patchwork") # 1.2.0
library(cols4all); packageVersion("cols4all") # 0.8

# Create output
(outdir <- macam::outdir_create())

# Read data
d <- read.csv("data_otherspp/data.csv")
d_amp <- read.csv("data_otherspp/data_amplification.csv")
head(d)

# Extract standard samples
d_std <- d %>% filter(task == "STANDARD")
d_stn <- d %>% filter(sample_name == "Stenella")
d_del <- d %>% filter(sample_name == "Delphinus")
d_tur <- d %>% filter(sample_name == "Tursiops")
d_tis <- d %>% filter(sample_name == "Tissue")

# ------------------------------------------------------- #
# Calculate DNA amount
# ------------------------------------------------------- #
# DNA copy * 2 ul / (Avogadro constant) * (Molecular Weight) * fg
d_std$pg <- d_std$quantity * 2 / (6.022 * 10^23) * (152 * 660) * 10^12
d_std$fg <- d_std$pg*1000 %>% round(digits = 3)


# ------------------------------------------------------- #
# Extract standard sample data
# ------------------------------------------------------- #
d_amp2 <- d_amp %>% filter(d_amp$Well %in% d_std$Well)
d_amp2$conc <- d_std[match(d_amp$Well, d_std$Well), "quantity"] %>% na.omit %>% as.numeric
d_amp2$pg <- d_std[match(d_amp$Well, d_std$Well), "pg"] %>% na.omit %>% as.numeric

# ------------------------------------------------------- #
# Visualize
# ------------------------------------------------------- #
g1 <- d_amp2 %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = factor(conc))) +
  geom_line() + scale_color_manual(values = rev(c4a("carto.bold")),
                                   name = "Copies/µL",
                                   labels = c("1","10","100","1000","10,000",
                                              "100,000","1,000,000","10,000,000",
                                              "100,000,000","1,000,000,000")) +
  ggtitle("Delta Rn")


# ------------------------------------------------------- #
# Closely related species
# ------------------------------------------------------- #
# Choose target wells
d_amp_stn <- d_amp %>% filter((d_amp$Well %in% d_std$Well) | (d_amp$Well %in% d_stn$Well))
d_amp_del <- d_amp %>% filter((d_amp$Well %in% d_std$Well) | (d_amp$Well %in% d_del$Well))
d_amp_tur <- d_amp %>% filter((d_amp$Well %in% d_std$Well) | (d_amp$Well %in% d_tur$Well))
# Assign sample names
d_amp_stn$sample <- d_amp_del$sample <- d_amp_tur$sample <- "standard"
d_amp_stn$sample[d_amp_stn$Well %in% d_stn$Well] <- "Stenella"
d_amp_del$sample[d_amp_del$Well %in% d_del$Well] <- "Delphinus"
d_amp_tur$sample[d_amp_tur$Well %in% d_tur$Well] <- "Tursiops"
d_amp_stn$sample <- factor(d_amp_stn$sample, levels = c("standard", "Stenella"))
d_amp_del$sample <- factor(d_amp_del$sample, levels = c("standard", "Delphinus"))
d_amp_tur$sample <- factor(d_amp_tur$sample, levels = c("standard", "Tursiops"))

# Stenella
g2 <- d_amp_stn %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = sample)) +
  geom_line() + scale_color_manual(values = c("gray60", "red2")) +
  ggtitle("Delta Rn (Stenella + Standard)")

# Delphinus
g3 <- d_amp_del %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = sample)) +
  geom_line() + scale_color_manual(values = c("gray60", "red2")) +
  ggtitle("Delta Rn (Delphinus + Standard)")

# Tursiops
g4 <- d_amp_tur %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = sample)) +
  geom_line() + scale_color_manual(values = c("gray60", "red2")) +
  ggtitle("Delta Rn (Tursiops + Standard)")



# ------------------------------------------------------- #
# Tissue samples and natural samples
# ------------------------------------------------------- #
# Choose target wells
d_amp_tis <- d_amp %>% filter((d_amp$Well %in% d_std$Well) | (d_amp$Well %in% d_tis$Well))
d_amp_tis$sample <- "standard"
d_amp_tis$sample[d_amp_tis$Well %in% d_tis$Well] <- "tissue"

# Tissue + Natural
g5 <- d_amp_tis %>% 
  ggplot(aes(x = Cycle, y = Delta.Rn, group = Well, color = sample)) +
  geom_line() + scale_color_manual(values = c("gray60", "pink")) +
  ggtitle("Delta Rn (Tissue + Standard)")



# ------------------------------------------------------- #
# Save results
# ------------------------------------------------------- #
# Save PDFs
ggsave(paste0(outdir, "/AmplificationPlot2.pdf"), plot = g1, width = 8, height = 6)

# Save R objects
saveRDS(list(g1, g2, g3, g4, g5), "../04_FigCode/data_robj/Fig_OtherSpp.obj")

# Save session information
macam::save_session_info()
