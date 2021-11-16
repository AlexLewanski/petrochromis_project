##############################################################################################
### SCRIPT NAME: admixture_abundance_mod.R
### PURPOSE: analyze and visualize the relationship between admixture vs. abundance using
###          generalized additive models (GAMs)
### PRODUCTS:
###     abundance_admixture_gam_9_8_2021.png: main text figure displaying the GAM results
###                                           examining relationship between admixture and
###                                           abundance and the abundance of each species
###     polyodon_gam_mod_fit_9_8_2021.png: supplementary figure of several model fit plots
###                                        for the polyodon GAM
###     kazumbe_gam_mod_fit_9_8_2021.png: supplementary figure of several model fit plots
###                                       for the kazumbe GAM
###     gam_summary_table.txt: supplementary table containing the details for each GAM.
###     survey_info_table.txt: supplementary table containing visual surveys for polyodon
###                            and kazumbe
##############################################################################################


#####################
### SCRIPT SET-UP ###
#####################

### Loading libraries ###
library(tidyverse)
library(sf)
library(sp)
library(ggmap)
library(ggsn)
library(cowplot)
library(ggrepel)
library(rgeos)
library(here)


### Loading data/files ###
#shapefile_path <- '/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/non_vcf_datafiles/shapefiles/'
#output_path <- '/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/'
africa_lake_shp <- read_sf(here('non_vcf_datafiles', 'shapefiles', 'africawaterbody', 'Africa_waterbody.shp'))
tz_river_shp <- read_sf(here('non_vcf_datafiles', 'shapefiles', 'AFRICOVER_TZ_RIVERS-shapefile', 'AFRICOVER_TZ_RIVERS.shp'))

cichlid_metadata <- read.csv(here('working_metadata_file', 'monster_23jul20.csv'), stringsAsFactors = FALSE) #metadata


### Mapping resources ###
#https://geodata.lib.berkeley.edu/catalog/AFRICOVER_TZ_RIVERS
#https://geocompr.github.io/post/2019/ggplot2-inset-maps/
#https://gis.stackexchange.com/questions/206929/r-create-a-boundingbox-convert-to-polygon-class-and-plot/206952



#africa_lake_shp_project <- spTransform(africa_lake_shp, CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
#africa_lake_shp_project <- st_transform(africa_lake_shp, CRS("+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"))
#africa_lake_shp_project <- st_transform(africa_lake_shp, "+proj=longlat +ellps=WGS84 +datum=WGS84")



#######################
### PROCESSING DATA ###
#######################

sampling_locations <- cichlid_metadata %>% 
  dplyr::filter(sciname2 %in% c('Petrochromis kazumbe', 'Petrochromis polyodon') & location != 'Harembe') %>% 
  select(location, long, lat) %>%
  distinct(location, .keep_all = TRUE) %>% 
  mutate(Region = case_when(location == 'Ulombola' ~ 'South',
                            location %in% c('Hilltop', 'Jakob', 'Bangwe') ~ 'Mid',
                            location %in% c('Nondwa', 'Kalala', 'Ska', 'Katongwe_S', 'Katongwe_N', 'Gombe_South') ~ 'North')) %>% 
  arrange(desc(lat)) %>% 
  mutate(number_label = LETTERS[1:n()],
         #number_label = 1:n(),
         Region = factor(Region, levels = c('North', 'Mid', 'South')))

tz_river_points <- sf::st_point_on_surface(tz_river_shp)

# retrieve the coordinates
tz_river_coords <- as.data.frame(sf::st_coordinates(tz_river_points))
tz_river_coords$NAME <- tz_river_shp$AUTO_ID

# ggplot() +
#   geom_sf(data = africa_lake_shp, fill = '#ade8f4', color = '#ade8f4', size = 3) +
#   geom_sf(data = tz_river_shp, fill = '#ade8f4', color = '#ade8f4', size = 1.5) +
#   geom_text(data = tz_river_coords, aes(X, Y, label = NAME), colour = "black") +
#   xlim(29.55, 29.79) + ylim(-5.03, -4.72)
# 
# ggplot() +
#   geom_sf(data = africa_lake_shp, fill = '#ade8f4', color = '#ade8f4', size = 3) +
#   geom_sf(data = tz_river_shp[tz_river_shp$AUTO_ID %in% c('5153', '5070'),], fill = '#ade8f4', color = '#ade8f4', size = 1.5) +
#   xlim(29.55, 29.79) + ylim(-5.03, -4.72)



#########################
### CREATING MAP PLOT ###
#########################

water_color <- '#caf0f8'

#lake tanganyika label: 29.61, -4.97
#Luiche River label: 29.72, -4.85

kigoma_region_map <- ggplot() +
  geom_sf(data = africa_lake_shp, fill = water_color, color = water_color, size = 2) +
  geom_sf(data = tz_river_shp[tz_river_shp$AUTO_ID %in% c('5153', '5070'),], fill = water_color, color = water_color, size = 1.5) +
  ggrepel::geom_label_repel(data = sampling_locations,
                            aes(x = long, y = lat, label = number_label),
                            size = 6.25, alpha = 1,
                            label.size = NA,
                            fill = NA,
                            #segment.color = "black", segment.size = 1,
                            box.padding = 0.009, #0.01
                            point.padding = 0,
                            #min.segment.length = 0.01, nudge_x = 0.02,
                            #min.segment.length = 1, nudge_x = 0.011,
                            min.segment.length = 1, nudge_x = -0.0112,
                            seed = 1003) +
  geom_text(data = data.frame(label = c('Lake Tanganyika', 'Luiche River'),
                              lon = c(29.61, 29.72),
                              lat = c(-4.97, -4.85)),
            aes(x = lon, y = lat, label = label, family = "Times", fontface = "plain"), size = 5.5, color = '#0077b6') +
  geom_point(data = sampling_locations, aes(x = long, y = lat, fill = Region), size = 8, pch=21, color = '#e7e7e6') +
  scale_fill_manual(values = c('#602437', '#b9375e', '#ff9ebb')) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.72, 1),  #c(0.95, 0.45) c(0.63, 0.173) --> lower left, next to scale bar
        legend.justification = c("right", "top"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.key=element_blank(),
        legend.background=element_blank()) +
  #coord_cartesian(clip = 'off') +
  xlim(29.55, 29.79) + ylim(-5.03, -4.72) + 
  scalebar(x.min = 29.55, x.max = 29.79, y.min = -5.03, y.max = -4.72, 
           dist = 5, dist_unit = "km",
           st.dist = .02,
           st.size = 4.5,
           height = .013,
           border.size = 0.6,
           transform = TRUE, model = "WGS84", st.bottom = FALSE, location = "bottomleft")

kigoma_box <- st_as_sfc(rgeos::bbox2SP(n = -4.72, s = -5.03, w = 29.55, e = 29.79,
                         proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")))

kigoma_inset_map <- ggplot() + 
  geom_sf(data = africa_lake_shp[africa_lake_shp$NAME_OF_WA == "Tanganyika", ], fill = water_color, color = water_color, size = 0.1) +
  geom_sf(data = kigoma_box, fill = NA, color = "red", size = 1) +
  #xlim(28.6, 31.3) + 
  ylim(-8.8, -3.3) +
  scale_x_continuous(breaks = c(29, 30, 31)) + 
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) + 
  theme_minimal() +
  theme(#panel.grid.major = element_blank(),
        axis.text = element_text(size = 9),
        panel.grid.major = element_line(color = '#e7e7e6'))

ggdraw() +
  draw_plot(kigoma_region_map, width = 1, height = 1) +
  draw_plot(kigoma_inset_map, x = 0.65, y = 0.64, width = 0.37, height = 0.37)

ggsave2(paste0(output_path, 'kigoma_sampling_map_point_outline.png'), width = 7, height = 8.75)

#/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/


