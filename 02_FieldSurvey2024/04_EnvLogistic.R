####
#### Lantau CWD qPCR v.s. Environment: Logistic regression
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0
library(patchwork); packageVersion("patchwork") # 1.2.0
library(ggrepel); packageVersion("ggrepel") # 0.9.5
library(ggforce); packageVersion("ggforce") # 0.4.2
library(ggbeeswarm); packageVersion("ggbeeswarm") # 0.7.2

# Create output directory
dir.create("04_EnvLogisticOut")

# Load data
d <- read.csv("data/data_qpcr_202407.csv")
d$date <- ymd(d$date)
llimit <- 2.875
d$detection1 <- (d$copies_ul > llimit)
d$detection2 <- d$copies_ul
d$detection2[d$copies_ul < llimit] <- 0 # Replace N.D. with 0
d_sur <- d %>% filter(depth == 0.2)
d_bot <- d %>% filter(depth > 0.2)


# ------------------------------------------ #
# qPCR data v.s. Environmental variables
# ------------------------------------------ #
# Extract subset
d_sub <- d %>% filter(date > "2024-07-28") %>%
  filter(sample_nc == "sample") %>% 
  select(sample_name, depth, sur_bot, temp, salinity, do, detection1)
d_long <- pivot_longer(d_sub, cols = -c(sample_name, depth, sur_bot))

# Logistic regression
## Temperature effect
m1 <- glm(detection1 ~ temp, data = d_sub, family = binomial)
m1_pred_df <- data.frame(temp = seq(25,30,by=0.5))
m1_pred_df$pred <- predict(m1, newdata = m1_pred_df, type = "response")
m1_pred_df$sur_bot <- NA
summary(m1)

# Visualize
g1 <- d_sub %>% 
  ggplot(aes(x = temp, y = as.numeric(detection1), color = sur_bot)) +
  geom_point() +
  geom_line(data = m1_pred_df, aes(x = temp, y = pred), color = "black") +
  annotate("text", x = 26, y = 0.5, label = "N.S. (P = 0.720)") +
  xlab(expression(paste("Temperature (", degree, "C)"))) +
  ylab("Detection probability") +
  coord_cartesian(xlim=c(25.5,29)) +
  NULL

## Salinity effect
m2 <- glm(detection1 ~ salinity, data = d_sub, family = binomial)
m2_pred_df <- data.frame(salinity = seq(12,33,by=0.5))
m2_pred_df$pred <- predict(m2, newdata = m2_pred_df, type = "response")
m2_pred_df$sur_bot <- NA
summary(m2)

# Visualize
g2 <- d_sub %>% 
  ggplot(aes(x = salinity, y = as.numeric(detection1), color = sur_bot)) +
  geom_point() +
  geom_line(data = m2_pred_df, aes(x = salinity, y = pred), color = "black") +
  annotate("text", x = 26, y = 0.5, label = "N.S. (P = 0.386)") +
  xlab(expression(paste("Salinity (\u2030)"))) +
  ylab("Detection probability") +
  coord_cartesian(xlim=c(12,33)) +
  NULL

## Salinity effect
m3 <- glm(detection1 ~ do, data = d_sub, family = binomial)
m3_pred_df <- data.frame(do = seq(2,7.5,by=0.2))
m3_pred_df$pred <- predict(m3, newdata = m3_pred_df, type = "response")
m3_pred_df$sur_bot <- NA
summary(m3)

# Visualize
g3 <- d_sub %>% 
  ggplot(aes(x = do, y = as.numeric(detection1), color = sur_bot)) +
  geom_point() +
  geom_line(data = m3_pred_df, aes(x = do, y = pred), color = "black") +
  annotate("text", x = 6, y = 0.5, label = "N.S. (P = 0.691)") +
  xlab("DO (mg/L)") +
  ylab("Detection probability") +
  coord_cartesian(xlim=c(2,7.5)) +
  NULL


# ------------------------------------------ #
# qPCR data v.s. surface or bottom
# cross-table (x-squared test)
# ------------------------------------------ #
chisq.test(table(d_sub$sur_bot, d_sub$detection1))

c1 <- d_sub %>% 
  ggplot(aes(x = sur_bot, y = as.numeric(detection1))) +
  geom_violin() +
  geom_beeswarm() +
  annotate("text", x = 1.5, y = 0.5, label = "N.S. (P = 1.0)") +
  xlab("Surface or bottom") +
  ylab("Detection probability") +
  NULL



# ------------------------------------------ #
# Save outputs
# ------------------------------------------ #
# Save figures
saveRDS(list(c1,g1,g2,g3), "../04_FigCode/data_robj/Fig_qPCRenv_summer_logis.obj")

# Save figures
Cairo::CairoPDF("04_EnvLogisticOut/qPCRxEnv_summary2.pdf", width = 10, height = 8)
(c1+g1)/(g2+g3)
dev.off()

# Save session information
macam::save_session_info()

