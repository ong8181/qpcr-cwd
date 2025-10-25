####
#### Lantau CWD qPCR v.s. Environment
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0
library(patchwork); packageVersion("patchwork") # 1.2.0
library(ggrepel); packageVersion("ggrepel") # 0.9.5

# Create output directory
dir.create("03_qPCRxEnvOut")

# Load data
d <- read.csv("data/data_qpcr_202501.csv")
d$date <- ymd(d$date)
llimit <- 2.195 # Max NC conc
d$detection1 <- (d$copies_ul > llimit)
d$detection2 <- d$copies_ul
d$detection2[d$copies_ul < llimit] <- 0 # Replace N.D. with 0
d_sur <- d %>% filter(depth == 0.2)
d_bot <- d %>% filter(depth > 0.2)


# ------------------------------------------ #
# qPCR data v.s. Environmental variables
# ------------------------------------------ #
# Extract subset
d_sub <- d %>% filter(date > "2025-01-16") %>%
  filter(sample_nc == "sample") %>% 
  select(sample_name, depth, sur_bot, temp, salinity, do, copies_ul)
d_long <- pivot_longer(d_sub, cols = -c(sample_name, depth, sur_bot))

# Visualize
## Surface v.s. Bottom
g1 <- d_sub %>%
  ggplot(aes(x = sur_bot, y = copies_ul*100)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(height = 0, width = 0.2) +
  xlab(NULL) +
  ylab("qPCR (copies/L)") +
  geom_hline(yintercept = llimit*100, linetype = 2) +
  geom_text(x = 0.8, y = llimit*100, label = "Detection limit", color = "black", size = 3.5)
## Temperature
g2 <- d_sub %>%
  ggplot(aes(x = temp, y = copies_ul*100, color = sur_bot)) +
  geom_point() + stat_smooth(method = "lm", se = FALSE) +
  geom_text_repel(aes(x = temp, y = copies_ul*100, label = sample_name), color = "black") +
  xlab(expression(paste("Temperature (", degree, "C)"))) +
  ylab("qPCR (copies/L)") +
  geom_hline(yintercept = llimit*100, linetype = 2) +
  geom_text(x = 0.8, y = llimit*100, label = "Detection limit", color = "black", size = 3.5)
## Salinity
g3 <- d_sub %>%
  ggplot(aes(x = salinity, y = copies_ul*100, color = sur_bot)) +
  geom_point() + stat_smooth(method = "lm", se = FALSE) +
  geom_text_repel(aes(x = salinity, y = copies_ul*100, label = sample_name), color = "black") +
  xlab(expression(paste("Salinity (\u2030)"))) +
  ylab("qPCR (copies/L)") +
  geom_hline(yintercept = llimit*100, linetype = 2) +
  NULL
## DO
g4 <- d_sub %>%
  ggplot(aes(x = do, y = copies_ul*100, color = sur_bot)) +
  geom_point() + stat_smooth(method = "lm", se = FALSE) +
  geom_text_repel(aes(x = do, y = copies_ul*100, label = sample_name), color = "black") +
  xlab("DO (mg/L)") +
  ylab("qPCR (copies/L)") +
  geom_hline(yintercept = llimit*100, linetype = 2) +
  NULL


# ------------------------------------------ #
# Save output
# ------------------------------------------ #
# Save figures
saveRDS(list(g1,g2,g3,g4), "../04_FigCode/data_robj/Fig_qPCRenv_winter.obj")

## cairo_pdf() for unicode
Cairo::CairoPDF("03_qPCRxEnvOut/qPCRxEnv_summary2.pdf", width = 10, height = 8)
(g1+g2)/(g3+g4)
dev.off()

# Save session information
macam::save_session_info()
