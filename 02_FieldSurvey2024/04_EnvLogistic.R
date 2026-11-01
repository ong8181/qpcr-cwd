####
#### Lantau CWD qPCR v.s. Environment: Logistic regression
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0
library(patchwork); packageVersion("patchwork") # 1.3.2
library(ggrepel); packageVersion("ggrepel") # 0.9.8
library(ggforce); packageVersion("ggforce") # 0.5.0
library(ggbeeswarm); packageVersion("ggbeeswarm") # 0.7.3
library(lme4); packageVersion("lme4") # 2.0.6

# Create output directory
(outdir <- macam::outdir_create())

# Load data
d <- read.csv("data/data_qpcr_202407.csv")
d$date <- ymd(d$date)
llimit <- 2.874 # Max of NC
d$detection1 <- ((d$copies_ul > llimit) & (d$n_pos_rep >= 2))
d$detection2 <- d$copies_ul
d$detection2[d$copies_ul < llimit] <- 0 # Replace N.D. with 0
d_sur <- d %>% filter(depth == 0.2)
d_bot <- d %>% filter(depth > 0.2)


# ------------------------------------------ #
# qPCR data v.s. Environmental variables
# ------------------------------------------ #
# Extract subset (add 'station' ID for paired structure)
d_sub <- d %>% filter(date > "2024-07-28") %>%
  filter(sample_nc == "sample") %>% 
  mutate(station = str_remove(sample_name, "-[SB]$")) %>%
  select(sample_name, station, depth, sur_bot, temp, salinity, do, detection1)
d_long <- pivot_longer(d_sub, cols = -c(sample_name, station, depth, sur_bot))

# Helper function to get GLMM predictions, 95% CI, and Odds Ratio label
get_glmm_pred <- function(model, new_df, var_name) {
  X <- model.matrix(as.formula(paste0("~ ", var_name)), data = new_df)
  beta <- fixef(model)
  vcov_mat <- as.matrix(vcov(model))
  fit_logit <- as.numeric(X %*% beta)
  se_logit  <- sqrt(diag(X %*% vcov_mat %*% t(X)))
  new_df$pred <- plogis(fit_logit)
  new_df$lwr  <- plogis(fit_logit - 1.96 * se_logit)
  new_df$upr  <- plogis(fit_logit + 1.96 * se_logit)
  new_df$sur_bot <- NA
  return(new_df)
}

get_glmm_stats <- function(model, var_name) {
  coefs <- summary(model)$coefficients
  b  <- coefs[var_name, "Estimate"]
  se <- coefs[var_name, "Std. Error"]
  p  <- coefs[var_name, "Pr(>|z|)"]
  or <- exp(b)
  or_lwr <- exp(b - 1.96 * se)
  or_upr <- exp(b + 1.96 * se)
  tibble(Variable = var_name,
         Slope_beta = b, SE = se,
         OR = or,
         OR_95CI_lower = or_lwr,
         OR_95CI_upper = or_upr,
         P_value = p,
         Label = sprintf("N.S. (P = %.3f)", p))
}

# Logistic regression (GLMM with station as random intercept)
## Temperature effect
m1 <- glmer(detection1 ~ temp + (1 | station), data = d_sub, family = binomial)
m1_pred_df <- get_glmm_pred(m1, data.frame(temp = seq(25, 30, by = 0.1)), "temp")
m1_stat <- get_glmm_stats(m1, "temp")
summary(m1)

# Visualize
g1 <- d_sub %>% 
  ggplot(aes(x = temp, y = as.numeric(detection1), color = sur_bot)) +
  geom_point() +
  geom_line(data = m1_pred_df, aes(x = temp, y = pred), color = "black") +
  annotate("text", x = 27.2, y = 0.5, label = m1_stat$Label, size = 3.2) +
  xlab(expression(paste("Temperature (", degree, "C)"))) +
  ylab("Detection probability") +
  coord_cartesian(xlim = c(25.5, 29), ylim = c(0, 1)) +
  NULL

## Salinity effect
m2 <- glmer(detection1 ~ salinity + (1 | station), data = d_sub, family = binomial)
m2_pred_df <- get_glmm_pred(m2, data.frame(salinity = seq(12, 33, by = 0.5)), "salinity")
m2_stat <- get_glmm_stats(m2, "salinity")
summary(m2)

# Visualize
g2 <- d_sub %>% 
  ggplot(aes(x = salinity, y = as.numeric(detection1), color = sur_bot)) +
  geom_point() +
  geom_line(data = m2_pred_df, aes(x = salinity, y = pred), color = "black") +
  annotate("text", x = 22.5, y = 0.5, label = m2_stat$Label, size = 3.2) +
  xlab(expression(paste("Salinity (\u2030)"))) +
  ylab("Detection probability") +
  coord_cartesian(xlim = c(12, 33), ylim = c(0, 1)) +
  NULL

## DO effect
m3 <- glmer(detection1 ~ do + (1 | station), data = d_sub, family = binomial)
m3_pred_df <- get_glmm_pred(m3, data.frame(do = seq(2, 7.5, by = 0.2)), "do")
m3_stat <- get_glmm_stats(m3, "do")
summary(m3)

# Visualize
g3 <- d_sub %>% 
  ggplot(aes(x = do, y = as.numeric(detection1), color = sur_bot)) +
  geom_point() +
  geom_line(data = m3_pred_df, aes(x = do, y = pred), color = "black") +
  annotate("text", x = 4.8, y = 0.5, label = m3_stat$Label, size = 3.2) +
  xlab("DO (mg/L)") +
  ylab("Detection probability") +
  coord_cartesian(xlim = c(2, 7.5), ylim = c(0, 1)) +
  NULL


# ------------------------------------------ #
# qPCR data v.s. surface or bottom
# Paired contingency table (McNemar's test)
# ------------------------------------------ #
# Reshape to station-level paired format (14 stations)
d_pair <- d_sub %>%
  select(station, sur_bot, detection1) %>%
  pivot_wider(names_from = sur_bot, values_from = detection1)

pair_tab <- table(Surface = factor(d_pair$surface, levels = c(TRUE, FALSE)),
                  Bottom  = factor(d_pair$bottom,  levels = c(TRUE, FALSE)))

# McNemar's test and exact binomial test on discordant pairs
mcnemar.test(pair_tab)
exact_p <- binom.test(pair_tab["TRUE", "FALSE"],
                      pair_tab["TRUE", "FALSE"] + pair_tab["FALSE", "TRUE"])$p.value

c1 <- d_sub %>% 
  ggplot(aes(x = sur_bot, y = as.numeric(detection1))) +
  geom_violin() +
  geom_beeswarm() +
  annotate("text", x = 1.5, y = 0.5,
           label = sprintf("N.S. (P = %.2f)", exact_p)) +
  xlab("Surface or bottom") +
  ylab("Detection probability") +
  NULL

# ------------------------------------------ #
# qPCR data v.s. surface or bottom
# cross-table (x-squared test; Old version)
# ------------------------------------------ #
# chisq.test(table(d_sub$sur_bot, d_sub$detection1))
# c1 <- d_sub %>% 
#   ggplot(aes(x = sur_bot, y = as.numeric(detection1))) +
#   geom_violin() +
#   geom_beeswarm() +
#   annotate("text", x = 1.5, y = 0.5, label = "N.S. (P = 1.0)") +
#   xlab("Surface or bottom") +
#   ylab("Detection probability") +
#   NULL


# ------------------------------------------ #
# Save outputs
# ------------------------------------------ #
# Save statistical summary table for the manuscript / Reviewer 2
env_stats_table <- bind_rows(m1_stat, m2_stat, m3_stat)
write.csv(env_stats_table, file.path(outdir, "qPCRxEnv_GLMM_summary.csv"), row.names = FALSE)

# Save figures
saveRDS(list(c1,g1,g2,g3), "../04_FigCode/data_robj/Fig_qPCRenv_summer_logis.obj")

# Save figures
## cairo_pdf() for unicode
Cairo::CairoPDF(file.path(outdir, "qPCRxEnv_summary2.pdf"), width = 10, height = 8)
(c1+g1)/(g2+g3)
dev.off()

# Save session information
macam::save_session_info()
