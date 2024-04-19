# Make map of samples 
# Read in metadata
# Make plot

# This script plot maps for the manuscript --------------------------#####
# Libraries ---------------------------------------------------------#####
library(ggplot2)
library(rgdal)
library(plyr)
library(tidyverse)
library(gridExtra)
library(cowplot)  
library(rnaturalearthdata)
# library(ggrepel)  # Not used in final version of making the plots as 
# labels were added manually
# Define inputs -----------------------------------------------------#####
# Paths to files may change depending on where it is stored in the computer
path2metadatafile <- "data/metadata/myna_WGS_meta_pop_env.csv"  # This path may change
path2shapefile <- "../shared_data/A_tristis_distribution_map/data_0.shp"  # This path will depend on where the shapefile is stored.

# Read metadata to get lat/lon --------------------------------------#####
meta_dt <- read.csv(path2metadatafile, header = T, na.strings = "n/a")
meta_dt$pop3_abb[meta_dt$pop3_abb == "Australia_unknown"] <- "Cairns"
# Calculate average lat/lon for plotting ----------------------------#####
meta_dt_pop <- meta_dt %>% 
  dplyr::group_by(pop3_abb) %>%
  dplyr::summarise(latitude = mean(as.numeric(latitude), na.rm = T),
                   longitude = mean(as.numeric(longitude), na.rm = T),
                   n = n())

meta_dt_pop$pop3_lab <- gsub(x = meta_dt_pop$pop3_abb, pattern = "_", " ", fixed = T)
meta_dt_pop$pop3_lab[meta_dt_pop$pop3_lab == "Maharashtra subpopulation A"] <- "Maharashtra\nsubpop. A"
meta_dt_pop$pop3_lab[meta_dt_pop$pop3_lab == "Great Barrier Island"] <- "Great Barrier\nIsland"

# Make map ----------------------------------------------------------#####
## Read shape file for Common myna distribution ---------------------#####
sf_basemap <- sf::read_sf("../shared_data/ne_10m_coastline/ne_10m_coastline.shp")

sf <- sf::read_sf(path2shapefile)
sf$range <- "Invasive"
sf$range[sf$ORIGIN == 1] <- "Native"
## Make base map ----------------------------------------------------#####
# Make base map of the continents and the distribution

worldmap <- ne_countries(scale = 'medium', type = 'map_units',
                         returnclass = 'sf')
mapWorld <- sf::st_as_sf(map_data('world'))
# map_dist <- geom_polygon(data = shapefile_df, aes(x = long.edit, y = lat, group = group, fill = ranges), alpha = 1) 
latbrks <- seq(-30, 30, by = 30)
lonbrks <- seq(60, 180, by = 60)

# ggplot() + 
  # geom_sf(data = sf, mapping = aes(group = ORIGIN))
crop_extent <- st_as_sfc(st_bbox(c(xmin = 10, xmax = 180, ymin = -50, ymax = 50), crs = st_crs(sf)))
st_points <- st_as_sf(meta_dt_pop, coords = c("longitude", "latitude"), crs = st_crs(sf))
# sfb <- st_intersection(sf_basemap, crop_extent)
sf_use_s2(FALSE)
mp <- ggplot() + 
  geom_sf(data = st_intersection(worldmap, crop_extent), fill = "lightgrey", colour = "lightgrey") +
  geom_sf(data = st_intersection(sf, crop_extent), mapping = aes(fill = range), alpha = 0.8) +
  geom_sf(data = st_points, size = 2, alpha = 0.8) +
  # geom_sf_label(data = st_points, aes(label = pop3_abb)) +
  # geom_text_repel(data = meta_dt_pop, aes(x = longitude, y = latitude, label = pop3_lab), box.padding = 0.5) +
  theme_bw() +
  theme(legend.justification = c(1,1), 
        legend.position = c(0.99,0.99),
        legend.direction = "horizontal",
        legend.text=element_text(size=15, face = "bold"),
        # legend.title=element_text(size=12, face = 'bold'),
        legend.text.align = 0,
        legend.background = element_rect(color = "black", size = 0.1, linetype = "solid"),
        axis.text = element_blank(),
        axis.line = element_line(colour = 'black'),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        plot.title = element_text(margin = margin(b = -25, unit = "pt"), size = 20),
        plot.margin = unit(c(0,20,0,5), 'pt'),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(size = 0.5, 
                                        linetype = "solid")) +
  coord_sf(xlim = c(10, 180), ylim = c(-50, 50), expand = F, ) +
  labs(fill = NULL) +
  scale_fill_manual(values = c("#d95f02","#377eb8")) 
# mp

png("results/sample_mapv1b.png", width = 8.3, height = 11.7/2, units = "in", res = 600)
mp
dev.off()



