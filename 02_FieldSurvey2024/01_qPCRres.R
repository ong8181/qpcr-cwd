####
#### Summer Lantau CWD qPCR
####

# Load library
library(tidyverse); packageVersion("tidyverse") # 2.0.0
library(patchwork); packageVersion("patchwork") # 1.2.0
library(ggmap); packageVersion("ggmap") # 4.0.0
library(ggrepel); packageVersion("ggrepel") # 0.9.5
library(ggtext); packageVersion("ggtext") # 0.1.2
library(rphylopic); packageVersion("rphylopic") # 1.5.0

# Create output directory
dir.create("01_qPCRresOut")

# Load data
d <- read.csv("data/data_qpcr_202407.csv")
llimit <- 2.875 # Max of NC
d$detection1 <- (d$copies_ul > llimit) # Max of NC
d$detection2 <- d$copies_ul
d$detection2[d$copies_ul < llimit] <- 0 # Replace N.D. with 0
d$date <- ymd(d$date)

d_op <- d %>% filter(date < "2024-07-28", sample_nc == "sample") # OP sapmles
d_sur <- d %>% filter(date > "2024-07-28") %>% 
  filter(sur_bot == "surface") # Currently exclude the OP samples
d_bot <- d %>% filter(sur_bot == "bottom")


# ------------------------------------------ #
# Visualize sampling locations
# ------------------------------------------ #
# Register AIP
ggmap::register_stadiamaps("xxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx")
# Draw Hong Kong Map
hk_lonlat <- c(left = 113.8, right = 114.45, top = 22.6, bottom = 22.1)
lt_lonlat <- c(left = 113.75, right = 114.1, top = 22.4, bottom = 22.12)
lt_lonlat2 <- c(left = 113.81, right = 113.95, top = 22.3, bottom = 22.18)
# Maptype "alidade_smooth", "outdoors", "stamen_terrain"
hk_map <- get_stadiamap(hk_lonlat, zoom = 11, maptype = "stamen_terrain")
lt_map <- get_stadiamap(lt_lonlat, zoom = 12, maptype = "stamen_terrain")
lt_map2 <- get_stadiamap(lt_lonlat2, zoom = 14, maptype = "stamen_terrain")
ggmap(lt_map)

# Sampling locations
g1 <- ggmap(lt_map2) +
  geom_point(data = d, aes(x = lon, y = lat), color = "red3", shape = 4) +
  geom_text_repel(data = d, aes(x = lon, y = lat, label = sample_name), max.overlaps = 25) +
  xlab(expression(paste("Longitude (", degree, "E)"))) +
  ylab(expression(paste("Latitude (", degree, "N)")))
g1_2 <- ggmap(lt_map2) +
  geom_point(data = d, aes(x = lon, y = lat), color = "red3", shape = 4) +
  xlab(expression(paste("Longitude (", degree, "E)"))) +
  ylab(expression(paste("Latitude (", degree, "N)")))


# ------------------------------------------ #
# Visualize dolphin sighting
# ------------------------------------------ #
# Prepare Phylopic
img <- pick_phylopic(name = "Sousa chinensis", n = 1)
img2 <- recolor_phylopic(img, fill = "pink2")
get_attribution("7bed394a-90e1-4da0-9d07-5b080a2061f4")

# Qualitative visualization
g2 <- ggmap(lt_map2) +
  geom_phylopic(aes(x = 113.849, y = 22.192), img = img2, fill = "pink2", color = "gray20", alpha = 0.75, size = 0.005) +
  geom_point(data = d, aes(x = lon, y = lat), color = "red3", shape = 4) +
  xlab(expression(paste("Longitude (", degree, "E)"))) +
  ylab(expression(paste("Latitude (", degree, "N)"))) +
  geom_text_repel(data = d, aes(x = lon, y = lat, label = sample_name), max.overlaps = 25) +
  NULL


# ------------------------------------------ #
# Visualize qPCR data
# ------------------------------------------ #
# Qualitative visualization
g3 <- ggmap(lt_map2) +
  geom_phylopic(aes(x = 113.849, y = 22.188), fill = "pink2", color = "gray20", alpha = 0.75, img = img2, size = 0.005) +
  geom_point(data = d_sur, aes(x = lon, y = lat, fill = detection1), size = 4, shape = 21) +
  geom_text_repel(data = d_sur, aes(x = lon, y = lat, label = sample_name %>% str_sub(end = 4)), max.overlaps = 25) +
  scale_fill_manual(name = "eDNA detected?", values = c("gray", "red3")) +
  xlab(expression(paste("Longitude (", degree, "E)"))) +
  ylab(expression(paste("Latitude (", degree, "N)"))) +
  ggtitle("Surface samples") +
  NULL

g4 <- ggmap(lt_map2) +
  geom_phylopic(aes(x = 113.849, y = 22.188), fill = "pink2", color = "gray20", alpha = 0.75, img = img2, size = 0.005) +
  geom_point(data = d_bot, aes(x = lon, y = lat, fill = detection1), size = 4, shape = 21) +
  geom_text_repel(data = d_bot, aes(x = lon, y = lat, label = sample_name %>% str_sub(end = 4)), max.overlaps = 25) +
  scale_fill_manual(name = "eDNA detected?", values = c("gray", "red3")) +
  xlab(expression(paste("Longitude (", degree, "E)"))) +
  ylab(expression(paste("Latitude (", degree, "N)"))) +
  ggtitle("Bottom samples") +
  NULL

# Quantitative visualization
g5 <- ggmap(lt_map2) +
  geom_phylopic(aes(x = 113.849, y = 22.188), fill = "pink2", color = "gray20", alpha = 0.75, img = img2, size = 0.005) +
  geom_point(data = d_sur, aes(x = lon, y = lat, fill = detection2*100, size = detection2*100), shape = 21) +
  geom_text_repel(data = d_sur, aes(x = lon, y = lat, label = sample_name %>% str_sub(end = 4)), max.overlaps = 25) +
  scale_fill_gradient(name = "eDNA\n(copies/L seawater)", limits = c(0,10.0)*100, low = "gray", high = "red3") +
  scale_size_continuous(name = "eDNA\n(copies/L seawater)", limits = c(0,10.0)*100) +
  xlab(expression(paste("Longitude (", degree, "E)"))) +
  ylab(expression(paste("Latitude (", degree, "N)"))) +
  ggtitle("Surface samples") +
  NULL

g6 <- ggmap(lt_map2) +
  geom_phylopic(aes(x = 113.849, y = 22.188), fill = "pink2", color = "gray20", alpha = 0.75, img = img2, size = 0.005) +
  geom_point(data = d_bot, aes(x = lon, y = lat, fill = detection2*100, size = detection2*100), shape = 21) +
  geom_text_repel(data = d_bot, aes(x = lon, y = lat, label = sample_name %>% str_sub(end = 4)), max.overlaps = 25) +
  scale_fill_gradient(name = "eDNA\n(copies/L seawater)", limits = c(0,10.00)*100, low = "gray", high = "red3") +
  scale_size_continuous(name = "eDNA\n(copies/L seawater)", limits = c(0,10.00)*100) +
  xlab(expression(paste("Longitude (", degree, "E)"))) +
  ylab(expression(paste("Latitude (", degree, "N)"))) +
  ggtitle("Bottom samples") +
  NULL


# ------------------------------------------ #
# Combine surface and bottom samples
# (Qualitative visualization only)
# ------------------------------------------ #
# Qualitative visualization
g7 <- ggmap(lt_map2) +
  geom_phylopic(aes(x = 113.849, y = 22.188), fill = "pink2", color = "gray20", alpha = 0.75, img = img2, size = 0.005) +
  geom_rect(data = d_sur, aes(xmin = lon - 0.005/1.9, xmax = lon + 0.005/1.9,
                              ymin = lat, ymax = lat + 0.6 * 0.0027,
                              fill = detection1), linewidth = 0.3, color = "black", alpha = 1) +
  geom_rect(data = d_bot, aes(xmin = lon - 0.005/1.9, xmax = lon + 0.005/1.9,
                              ymin = lat, ymax = lat - 0.6 * 0.0027,
                              fill = detection1), linewidth = 0.3, color = "black", alpha = 1) +
  geom_text_repel(data = d_sur, aes(x = lon, y = lat, label = sample_name %>% str_sub(end = 4)),
                  nudge_x = -0.008, max.overlaps = 25, min.segment.length = 2) +
  scale_fill_manual(name = "eDNA detected?", values = c("gray", "red3")) +
  xlab(expression(paste("Longitude (", degree, "E)"))) +
  ylab(expression(paste("Latitude (", degree, "N)"))) +
  ggtitle("Summer survey in 2024.7") +
  NULL


# ------------------------------------------ #
# Save figures
# ------------------------------------------ #
# Save Robj
saveRDS(list(g5,g6,g7), "../04_FigCode/data_robj/Fig_Survey_summer.obj")

# Save figures
ggsave("01_qPCRresOut/LantauSummer_map.pdf", g1, width = 10, height = 8)
ggsave("01_qPCRresOut/LantauSummer_map_CWD.pdf", g2, width = 10, height = 8)
ggsave("01_qPCRresOut/qPCR_detection_summer.pdf", g3 + g4, width = 15, height = 8)
ggsave("01_qPCRresOut/qPCR_detection_summer_comb.pdf", g7, width = 9, height = 8)
ggsave("01_qPCRresOut/qPCR_concs_summer.pdf", g5 + g6, width = 15, height = 8)

# Save session information
macam::save_session_info()
