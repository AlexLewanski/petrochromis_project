##############################################################################################
### SCRIPT NAME: studysystem_overview_fig.R
### PURPOSE: processing and visualizing pop gen statistics and creating study overview plot
###          for the main text
### PRODUCT:
###     prac_full_multipanel_5_12_2021.png: study overview figure that includes a map of the
###                                         study regions and sampling locations, plots of
###                                         pop gen stats (Fst, nucleotide diversity, dxy),
###                                         images of the two study species, and PCA plots
##############################################################################################


#####################
### SCRIPT SET-UP ###
#####################

### Loading libraries ###
library(tidyverse)
library(RColorBrewer)
library(egg)
library(sf)
library(sp)
library(ggmap)
library(ggsn)
library(cowplot)
library(ggrepel)
library(rgeos)
library(magick)
library(rgdal) #for adding bathymetry vis
library(maptools) #for adding bathymetry vis
library(here)

### Resources ###
#plotting
#https://stackoverflow.com/questions/27690729/greek-letters-symbols-and-line-breaks-inside-a-ggplot-legend-label
#https://stackoverflow.com/questions/42438450/make-legend-invisible-but-keep-figure-dimensions-and-margins-the-same
#https://stackoverflow.com/questions/47614314/how-to-put-plots-without-any-space-using-plot-grid
#https://stackoverflow.com/questions/57168938/ggplot-italicize-part-of-legend-across-multiple-lines

#mapping
#https://geodata.lib.berkeley.edu/catalog/AFRICOVER_TZ_RIVERS


### Loading data/files ###
cichlid_metadata <- read.csv(here('working_metadata_file', 'monster_23jul20.csv'), stringsAsFactors = FALSE) #metadata

#popgen stats
Fst_dataframe_updated <- read.csv(here('diversity_divergence_stats', 'processed_results', 'fst_7_2021.csv'))
pi_info <- read.csv(here('diversity_divergence_stats', 'processed_results', 'pi_7_2021.csv'))
dxy_info <- read.csv(here('diversity_divergence_stats', 'processed_results', 'dxy_7_2021.csv'))

#PCAs
both_species_pca <- read.csv(here('PCAs', 'kazumbe_polyodon_PCA_withmetadata_9_7_2021.csv'))
polyodon_pca <- read.csv(here('PCAs', 'polyodon_PCA_withmetadata_9_7_2021.csv'))
kazumbe_pca <- read.csv(here('PCAs', 'kazumbe_PCA_withmetadata_9_7_2021.csv'))
pca_importance_info <- read.csv(here('PCAs', 'pca_importance_info.csv'))

#map
tz_river_shp <- read_sf(here('non_vcf_datafiles', 'shapefiles', 'AFRICOVER_TZ_RIVERS-shapefile/AFRICOVER_TZ_RIVERS.shp'))
africa_lake_shp <- read_sf(here('non_vcf_datafiles', 'shapefiles', 'africawaterbody', 'Africa_waterbody.shp'))
contours <- readOGR(dsn=path.expand(here('non_vcf_datafiles', 'shapefiles', 'GIS_Package')),
                    layer="tcarta_tanganyika_contours",
                    stringsAsFactors=TRUE,
                    integer64="warn.loss")

#fish images
fish_image_list <- list(
  kazumbe = image_read(here('Fish_images', 'kazumbe_cropped.tif')),
  polyodon = image_read(here('Fish_images', 'polyodon_pic_cropped.tif'))
)


### Creating color palette ###
location_palette <- brewer.pal(11, "Set3")
names(location_palette) <- c("Hilltop", "Jakobsen's S", "Nondwa", "Kalalangabo", "Ulombola", "Katongwe N", "Gombe S", "Katongwe S", "Bangwe", "S Kagango", "Harembe")
location_palette[location_palette == '#FFFFB3'] <- '#F7DC6F' #switch the yellow out
location_palette[location_palette == "#D9D9D9"] <- "#AAB7B8" #switch the gray out
location_palette[names(location_palette) == "Hilltop"] <- "#996633"
location_palette[names(location_palette) == "Ulombola"] <- "#0000e6"
location_palette[names(location_palette) == "Nondwa"] <- "#4dd2ff"
location_palette[names(location_palette) == "Gombe S"] <- "#66ff66"
location_palette[names(location_palette) == "Harembe"] <- "#99ffeb"

#switch the colors of Katongwe N and Ulombola
Ulombola_col_switch <- location_palette[names(location_palette) == "Ulombola"]
Katongwe_N_col_switch <- location_palette[names(location_palette) == "Katongwe N"]

location_palette[names(location_palette) == "Ulombola"] <- Katongwe_N_col_switch
location_palette[names(location_palette) == "Katongwe N"] <- Ulombola_col_switch

#switch Jakobsen's S and Bangwe
Jakob_col_switch <- location_palette[names(location_palette) == "Jakobsen's S"]
Bangwe_col_switch <- location_palette[names(location_palette) == "Bangwe"]

location_palette[names(location_palette) == "Jakobsen's S"] <- Bangwe_col_switch
location_palette[names(location_palette) == "Bangwe"] <- Jakob_col_switch

location_palette_letter <- location_palette

location_letter_df <- data.frame(location = c('Gombe S', 'Katongwe N', 'Katongwe S', 'S Kagango', 'Kalalangabo', 'Nondwa', 'Hilltop', 'Bangwe', "Jakobsen's S", 'Ulombola', 'Harembe'),
                                 loc_letter = LETTERS[1:11])

for (i in seq_along(location_palette_letter)) {
  focal_loc <- names(location_palette_letter)[i]
  focal_let <- location_letter_df$loc_letter[location_letter_df$location == focal_loc]
  
  #replace name with associated letter
  names(location_palette_letter)[i] <- focal_let
}



################################
### VISUALIZING POPGEN STATS ###
################################

pi_plot <- ggplot(data = pi_info) +
  #geom_pointrange(aes(x = population, y = mean, ymin = lower, ymax = upper, group = species, color = species), size = 1.2) +
  geom_linerange(aes(x = population, ymin = lower, ymax = upper, group = species, color = species), size = 1.5) +
  geom_point(aes(x = population, y = mean, group = species, color = species), size = 5) +
  scale_color_manual(values = c("#DC7633", "#3498DB")) +
  #ylim(0, max(pi_info$upper)) +
  xlab("Population") + ylab("pi") +
  theme_cowplot() +
  theme(panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
        legend.position = "none") +
  coord_flip() +
  ggtitle('pi')

#combined plot of pi and dxy
divergence_diversity_plot <- ggplot(data = pi_info %>%
                                      mutate(Population = factor(Population, levels = (c('Full', LETTERS[10:1]) )) ) %>% 
                                      mutate(metric = case_when(species == 'kazumbe' ~ '\u03C0 (kazumbe)',
                                                                species == 'polyodon' ~ '\u03C0 (polyodon)')) %>% 
                                      select(-c(pop, species)) %>%
                                      rbind(., dxy_info %>% 
                                              select(-comparison_id) %>% 
                                              mutate(metric = 'dxy')) %>% 
                                      mutate(metric = factor(metric, levels = c('\u03C0 (kazumbe)', '\u03C0 (polyodon)', 'dxy')))) +
  geom_linerange(aes(x = Population, ymin = lower, ymax = upper, group = metric, color = metric), size = 2) +
  geom_point(aes(x = Population, y = mean, group = metric, color = metric), shape = 21, stroke = 2, size = 5.5) +
  #scale_color_manual(values = c("#3498DB", "#DC7633", "#6c757d"),
  #                   labels = c(expression(paste('\u03C0 (', italic('P. polyodon'), ')')), 
  #                              expression(paste('\u03C0 (', italic('P. kazumbe'), ')')), 'dxy' ) ) +
  scale_color_manual(values = c("#DC7633", "#3498DB", "#6c757d"),
                     labels = c(expression(paste('\u03C0 (', italic("P"), ". sp. 'kazumbe'", ')')),
                                expression(paste('\u03C0 (', italic("P"), ". cf. ",  italic("polyodon"), ')')),
                                expression(d[XY])) ) +
  #ylim(0, max(pi_info$upper)) +
  xlab("Population") + ylab("Average sequence divergence") +
  theme_cowplot() +
  theme(plot.margin = margin(5.5, 0, 5.5, 5.5),
        panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
        panel.grid.major.x = element_line(colour = 'gray', size = 0.4),
        legend.title = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        legend.position="bottom",
        legend.justification = "center",
        legend.spacing.x = unit(0.2, 'cm'),
        legend.text = element_text(size = 14, margin = margin(r = 4, unit = "pt"))) +
  #theme(panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
  #      legend.position = "none") +
  coord_flip()

Fst_plot <- ggplot(data = Fst_dataframe_updated %>% 
                     mutate(group = 'a',
                            Population = factor(Population, levels = (c('Full', LETTERS[10:1]) )))) +
  #geom_pointrange(aes(x = population, y = mean, ymin = lower, ymax = upper, group = species, color = species), size = 1.2) +
  geom_linerange(aes(x = Population, ymin = lower_CI, ymax = upper_CI, color = group), size = 2) +
  geom_point(aes(x = Population, y = Fst, color = group), shape = 21, stroke = 2, size = 5.5) +
  #ylim(0, max(pi_info$upper)) +
  xlab("Population") + 
  #ylab(expression(italic(F[ST]))) +
  ylab(expression(F[ST])) +
  theme_cowplot() +
  theme(plot.margin = margin(5.5, 5.5, 0, 5.5),
        panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
        panel.grid.major.x = element_line(colour = 'gray', size = 0.4),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_text(size = 20),
        axis.text.x = element_text(size = 15),
        legend.position = "bottom",
        legend.justification = "left",
        legend.text = element_text(size = 14, color = alpha("black", 0)),
        legend.title = element_text(color = "white"),
        legend.key = element_rect(fill = "white")) +
  scale_x_discrete(drop = FALSE) +
  scale_color_manual(values = c('#5e60ce'),
                     guide = guide_legend(override.aes = list(color = "white"))) +
  #scale_color_discrete(guide = guide_legend(override.aes = list(color = "white"))) +
  coord_flip()

#multipanel of pi/dxy plot and Fst plot
diversity_multipanel <- ggarrange(divergence_diversity_plot, Fst_plot,
                                  ncol = 2)



#################################
### CREATING STUDY REGION MAP ###
#################################

### Map set-up ###
#water_color <- '#caf0f8'
water_color <- '#ADE4FF'

#lake tanganyika label: 29.61, -4.97
#Luiche River label: 29.72, -4.85

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

#processing bathymetry data for inset map
bf.dataframe <- data.frame(id=rownames(contours@data),
                           values=sample(1:10,length(contours),replace=T),
                           contours@data, stringsAsFactors=F)
bf.fort   <- fortify(contours)
bf.merged <- plyr::join(bf.fort, bf.dataframe, by="id")

# retrieve the coordinates
tz_river_coords <- as.data.frame(sf::st_coordinates(tz_river_points))
tz_river_coords$NAME <- tz_river_shp$AUTO_ID

#kigoma_box <- st_as_sfc(rgeos::bbox2SP(n = -4.72, s = -5.03, w = 29.55, e = 29.79,
#                                       proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")))

kigoma_box <- st_as_sfc(rgeos::bbox2SP(n = -4.72, s = -5.03, w = 29.48, e = 29.79,
                                       proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")))


### Creating map ###
kigoma_region_map <- ggplot() +
  geom_sf(data = africa_lake_shp, fill = water_color, color = water_color, size = 2) +
  geom_sf(data = tz_river_shp[tz_river_shp$AUTO_ID %in% c('5153', '5070'),], fill = water_color, color = water_color, size = 1.5) +
  ggrepel::geom_label_repel(data = sampling_locations,
                            aes(x = long, y = lat, label = number_label),
                            size = 7, alpha = 1,
                            label.size = NA,
                            fill = NA,
                            #segment.color = "black", segment.size = 1,
                            box.padding = 0.009, #0.01
                            point.padding = 0,
                            #min.segment.length = 0.01, nudge_x = 0.02,
                            #min.segment.length = 1, nudge_x = 0.011,
                            #min.segment.length = 1, nudge_x = -0.0112,
                            min.segment.length = 1, nudge_x = -0.0132,
                            seed = 1002) + #1003
  geom_text(data = data.frame(label = c('Lake Tanganyika', 'Luiche River'),
                              lon = c(29.61, 29.705),
                              lat = c(-4.97, -4.847)),
            aes(x = lon, y = lat, label = label, family = "Times", fontface = "plain"), size = 8, color = '#0077b6') +
  geom_point(data = sampling_locations, aes(x = long, y = lat, fill = Region), size = 8.5, pch=21, color = '#e7e7e6') +
  scale_fill_manual(values = c('#602437', '#b9375e', '#ff9ebb')) +
  theme(panel.background = element_rect(fill = "white"),
        panel.grid.major = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = c(0.72, 1),  #c(0.95, 0.45) c(0.63, 0.173) --> lower left, next to scale bar
        legend.justification = c("right", "top"),
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        legend.key=element_blank(),
        legend.background=element_blank()) +
  #coord_cartesian(clip = 'off') +
  #xlim(29.55, 29.79) + ylim(-5.03, -4.72) + 
  xlim(29.55, 29.80) + ylim(-5.05, -4.68) + 
  scalebar(x.min = 29.55, x.max = 29.79, y.min = -5.03, y.max = -4.72, 
           dist = 5, dist_unit = "km",
           st.dist = .02,
           st.size = 4.5,
           height = .013,
           border.size = 0.6,
           transform = TRUE, model = "WGS84", st.bottom = FALSE, location = "bottomleft")

#inset map of kigoma
kigoma_inset_map <- ggplot() + 
  geom_sf(data = africa_lake_shp[africa_lake_shp$NAME_OF_WA == "Tanganyika",], fill = water_color, color = water_color, size = 0.1) +
  #xlim(28.6, 31.3) + 
  ylim(-8.8, -3.3) +
  scale_x_continuous(breaks = c(29, 30, 31)) + 
  #  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) + 
  theme_minimal() +
  theme(#panel.grid.major = element_blank(),
    axis.text = element_text(size = 9),
    panel.grid.major = element_line(color = '#e7e7e6'),
    axis.title = element_blank())  +
  geom_polygon(data=bf.merged[bf.merged$ELEVATION %in% c(-300),], size=0.4,
               aes(x=long,y=lat,group=group), fill='#5CC6FF') +
  geom_polygon(data=bf.merged[bf.merged$ELEVATION %in% c(-600),], size=0.4,
               aes(x=long,y=lat,group=group), fill='#0AA9FF') +
  geom_polygon(data=bf.merged[bf.merged$ELEVATION %in% c(-900),], size=0.4,
               aes(x=long,y=lat,group=group), fill='#0092E0') +
  geom_polygon(data=bf.merged[bf.merged$ELEVATION %in% c(-1200),], size=0.4,
            aes(x=long,y=lat,group=group), fill='#0077b6') +
  geom_sf(data = kigoma_box, fill = NA, color = "red", size = 1.1, alpha = 1)
  #geom_path(data=bf.merged[bf.merged$ELEVATION %in% c(-100, -400, -800, -1200),], size=0.4,
  #          aes(x=long,y=lat,group=group),color=scales::alpha('#0077b6',0.6))

tanganyika_sampling_map_final <- ggdraw() +
  draw_plot(kigoma_region_map, width = 1, height = 1) +
  #draw_plot(kigoma_inset_map, x = 0.65, y = 0.64, width = 0.4, height = 0.37) +
  draw_plot(kigoma_inset_map, x = 0.65, y = 0.56, width = 0.47, height = 0.47) +
  #theme(plot.margin = margin(5.5, -5, 5.5, 5.5))
  theme(plot.margin = margin(5.5, 10, 5.5, -2))

multipanel_diversity_toprow <- plot_grid(tanganyika_sampling_map_final, ggdraw(diversity_multipanel),
                                         ncol = 2, rel_widths = c(0.42,0.5), labels = c('(a)', '(b)'), label_size = 32, vjust = 0, hjust = c(-0.1, 0.1) ) + #rel_widths = c(0.45,0.5)
  theme(plot.margin = margin(32, 0, 8, 0))



###################
### FISH IMAGES ###
###################

kazumbe_plot <- ggplot() +
  ggtitle(expression(paste(italic("Petrochromis"), " sp. 'kazumbe'", sep = ""))) +
  draw_image(
    fish_image_list$kazumbe, scale = 1, x = 1,
    hjust = 1, halign = 0.5, valign = 0
  ) +
  theme_void() +
  theme(#plot.margin=unit(c(0.75, 0.75, 0.4, 0),"cm"),
    plot.margin=unit(c(0.4, 0, 0.2, 0),"cm"),
    plot.title = element_text(size = 21, face = 'italic', hjust = 0.5)) +
  coord_cartesian(clip = 'off')

polyodon_plot <- ggplot() +
  #ggtitle('Petrochromis polyodon') +
  ggtitle(expression(paste(italic("Petrochromis"), " cf. ",  italic("polyodon"), sep = ""))) +
  draw_image(
    fish_image_list$polyodon, scale = 1, x = 1,
    hjust = 1, halign = 0.5, valign = 0
  ) +
  theme_void() +
  theme(#plot.margin=unit(c(0.65, 0.75, 0.5, 0),"cm"),
    plot.margin=unit(c(0.4, 0, 0.2, 0),"cm"),
    plot.title = element_text(size = 21, face = 'italic', hjust = 0.5)) +
  coord_cartesian(clip = 'off')

#create multipanel of fish images
fish_multipanel <- plot_grid(kazumbe_plot, polyodon_plot, nrow = 2) +
  theme(plot.margin = unit(c(0.2, -0.2, 0.1, 0.5), "cm"))



#####################################################
### PLOTTING PRINCIPAL COMPONENTS ANALYSES (PCAs) ###
#####################################################

scaleFUN <- function(x) sprintf("%.2f", x)
empty_plot <- ggplot() + theme_void()

#PCA containing both species
species_PCA <- ggplot(data = both_species_pca, aes(x = PC1, y= PC2, color = species_ID_entropy)) +
  #geom_point(size = 5) +
  geom_point(size = 4.5, alpha = 0.5) +
  geom_point(size = 4.5, shape = 21, alpha = 0.7, stroke = 1) +
  scale_color_manual(name = "Species", values = c("#DC7633", "#3498DB"),
                     labels = c(expression(paste(italic("P"), ". sp. 'kazumbe'", sep = "")), 
                                expression(paste(italic("P"), ". cf. ",  italic("polyodon"), sep = "")))) +
  ###scale_color_manual(name = "Species", values = c("#ffc966", "#ccccff")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 13),
        axis.title.y = element_text(margin = margin(0, 0, 0, 0)), 
        #legend.title = element_text(size = 18),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.justification = "center",
        legend.margin = margin(-6, 0, 0, 0),
        legend.text = element_text(size = 14, margin = margin(r = 4, unit = "pt")),
        plot.margin = margin(5, 9, 5, 9)) +
  scale_x_continuous(labels = scaleFUN) +
  scale_y_continuous(labels = scaleFUN) +
  xlab(paste0("PC1 (", round(pca_importance_info[pca_importance_info$dataset == "combined_kazumbe_polyodon" & pca_importance_info$pc == 'PC1',]$prop_variance*100, 2),"%)")) +
  ylab(paste0("PC2 (", round(pca_importance_info[pca_importance_info$dataset == "combined_kazumbe_polyodon" & pca_importance_info$pc == 'PC2',]$prop_variance*100, 2),"%)")) +
  coord_cartesian(clip = 'off')

kazumbe_pca_plus_locletters <- kazumbe_pca %>% 
  mutate(location_factor_letter = factor(case_when(location == 'Gombe_South' ~ "A",
                                                   location == 'Katongwe_N' ~ "B",
                                                   location == 'Katongwe_S' ~ "C",
                                                   location == 'Ska' ~ "D",
                                                   location == 'Kalala' ~ "E",
                                                   location == 'Nondwa' ~ "F",
                                                   location == 'Hilltop' ~ "G",
                                                   location == 'Bangwe' ~ "H",
                                                   location == 'Jakob' ~ "I",
                                                   location == 'Ulombola' ~ "J"),
                                         levels= LETTERS[1:10]))

#PCA containing only kazumbe
kazumbe_PCA_plot <- ggplot(data = kazumbe_pca_plus_locletters, 
                           aes(x = PC1, y= PC2, color = location_factor_letter)) +
  geom_point(size = 4.5, alpha = 0.5) +
  geom_point(size = 4.5, shape = 1, alpha = 0.5) +
  scale_color_manual(name = "Location", values = location_palette_letter[names(location_palette_letter) %in% unique(kazumbe_pca_plus_locletters$location_factor_letter) ]) +
  ###scale_color_manual(name = "Location", values = gradient_palette[names(gradient_palette) %in% unique(kazumbe_PCA_combinedinfo$location_plotting) ]) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 13),
        axis.title.y = element_text(margin = margin(0, 0, 0, 0)), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 22, face="italic"),
        plot.margin = margin(5, 5, 5, 9)) +
  ggtitle(expression(paste(italic("P"), ". sp. 'kazumbe'", sep = ""))) +
  #ggtitle("P. kazumbe") + 
  scale_x_continuous(labels = scaleFUN) +
  scale_y_continuous(labels = scaleFUN) +
  xlab(paste0("PC1 (", round(pca_importance_info[pca_importance_info$dataset == "kazumbe" & pca_importance_info$pc == 'PC1',]$prop_variance*100, 2),"%)")) +
  ylab(paste0("PC2 (", round(pca_importance_info[pca_importance_info$dataset == "kazumbe" & pca_importance_info$pc == 'PC2',]$prop_variance*100, 2),"%)"))

polyodon_pca_plus_locletters <- polyodon_pca %>% 
  mutate(location_factor_letter = factor(case_when(location == 'Gombe_South' ~ "A",
                                                   location == 'Katongwe_N' ~ "B",
                                                   location == 'Katongwe_S' ~ "C",
                                                   location == 'Ska' ~ "D",
                                                   location == 'Kalala' ~ "E",
                                                   location == 'Nondwa' ~ "F",
                                                   location == 'Hilltop' ~ "G",
                                                   location == 'Bangwe' ~ "H",
                                                   location == 'Jakob' ~ "I",
                                                   location == 'Ulombola' ~ "J"),
                                         levels= LETTERS[1:10]))
#PCA containing only polyodon
polyodon_PCA_plot <- ggplot(data = polyodon_pca_plus_locletters, 
                            aes(x = PC1, y= PC2, color = location_factor_letter)) +
  geom_point(size = 4.5, alpha = 0.5) +
  geom_point(size = 4.5, shape = 1, alpha = 0.5) +
  scale_color_manual(name = "Location", values = location_palette_letter[names(location_palette_letter) %in% unique(polyodon_pca_plus_locletters$location_factor_letter) ]) +
  ###scale_color_manual(name = "Location", values = gradient_palette[names(gradient_palette) %in% unique(polyodon_PCA_combinedinfo$location_plotting) ]) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_text(size = 13),
        axis.title.y = element_text(margin = margin(0, 0, 0, 0)), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 22, face= c("italic") ),
        plot.margin = margin(5, 5, 5, 9)) +
  ggtitle(expression(paste(italic("P"), ". cf. ",  italic("polyodon"), sep = ""))) +
  #ggtitle("P. polyodon") + 
  scale_x_continuous(labels = scaleFUN) +
  scale_y_continuous(labels = scaleFUN) +
  xlab(paste0("PC1 (", format(round(pca_importance_info[pca_importance_info$dataset == "polyodon" & pca_importance_info$pc == 'PC1',]$prop_variance*100, 2), nsmall = 2) ,"%)")) +
  ylab(paste0("PC2 (", round(pca_importance_info[pca_importance_info$dataset == "polyodon" & pca_importance_info$pc == 'PC2',]$prop_variance*100, 2),"%)"))

both_species_pca_with_locletters <- both_species_pca %>% 
  mutate(Location = factor(case_when(location == 'Gombe_South' ~ "A",
                                     location == 'Katongwe_N' ~ "B",
                                     location == 'Katongwe_S' ~ "C",
                                     location == 'Ska' ~ "D",
                                     location == 'Kalala' ~ "E",
                                     location == 'Nondwa' ~ "F",
                                     location == 'Hilltop' ~ "G",
                                     location == 'Bangwe' ~ "H",
                                     location == 'Jakob' ~ "I",
                                     location == 'Ulombola' ~ "J"),
                           levels= LETTERS[1:10]))

#extracting legend
location_palette_letter_processed <- location_palette_letter[names(location_palette_letter) %in% unique(both_species_pca_with_locletters$Location) ]
location_palette_letter_processed <- location_palette_letter_processed[match(LETTERS[1:10], names(location_palette_letter_processed))]

location_legend <- get_legend(ggplot(data = both_species_pca_with_locletters,  
                                     aes(x = PC1, y= PC2, color = Location)) +
                                geom_point(size = 4.5, alpha = 0.5) +
                                geom_point(size = 4.5, shape = 1, alpha = 0.5) +
                                scale_color_manual(name = "Location", values = location_palette_letter_processed) +
                                ##scale_color_manual(name = "Location", values = gradient_palette[names(gradient_palette) %in% unique(petro_PCA_combinedinfo$location_plotting) ]) +
                                theme_bw() +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      axis.title = element_text(size = 18),
                                      plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
                                      legend.title = element_text(size = 18),
                                      legend.text = element_text(size = 16),
                                      #legend.text = element_text(size = 18, margin = margin(r = 22, unit = "pt")),
                                      legend.position = "right",
                                      #legend.spacing.x = unit(0.3, 'cm')
                                )
)

#multipanel of PCA plots
pca_multipanel <- plot_grid(species_PCA, 
                            kazumbe_PCA_plot, 
                            polyodon_PCA_plot,
                            ncol = 3, align = 'hv', axis = "bt")



###############################################
### CONSTRUCTING FULL STUDY OVERVIEW FIGURE ###
###############################################

multipanel_diversity_bottomrow <- plot_grid(fish_multipanel, pca_multipanel, plot_grid(location_legend, empty_plot, nrow = 2, rel_heights = c(0.9, 0.107)), 
                                            ncol = 3, rel_widths = c(0.35, 0.95, 0.125),
                                            labels = c('(c)', '(d)',''), label_size = 32, vjust = 1, hjust = c(-0.1, -0.3, -0.5) ) +
  theme(plot.margin = margin(15, 0, 0, 0))

plot_grid(multipanel_diversity_toprow,
          multipanel_diversity_bottomrow, nrow = 2, rel_heights = c(0.7, 0.3))

ggsave(here('figures', 'prac_full_multipanel_5_12_2021.png'),
        width = 19, height = 17, bg = 'white')
#width = 19, height = 16, bg = 'white')



#################################
### CODE CURRENTLY NOT IN USE ###
#################################

# diversity_calc <- function(input_data,
#                            calc_type = c('raw', 'per_site'),
#                            no_sites_col,
#                            count_diffs_col,
#                            count_comparisons_col) {
# 
#   if (calc_type == 'raw') {
#     return( sum(input_data[, count_diffs_col], na.rm = TRUE)/sum(input_data[, count_comparisons_col], na.rm = TRUE) )
#   } else if (calc_type == 'per_site') {
#     return( (sum(input_data[, count_diffs_col], na.rm = TRUE)/sum(input_data[, count_comparisons_col], na.rm = TRUE))/sum(input_data[, no_sites_col], na.rm = TRUE) )
#   }
# }
# 
# 
# 
# 
# divergence_diversity_plot_horizontal <- ggplot(data = pi_info %>%
#                                                  mutate(metric = case_when(species == 'kazumbe' ~ '\u03C0 (kazumbe)',
#                                                                            species == 'polyodon' ~ '\u03C0 (polyodon)')) %>% 
#                                                  select(-c(pop, species)) %>%
#                                                  rbind(., dxy_info %>% 
#                                                          select(-comparison_id) %>% 
#                                                          mutate(metric = 'dxy')) %>% 
#                                                  mutate(metric = factor(metric, levels = c('\u03C0 (polyodon)', '\u03C0 (kazumbe)', 'dxy')))) +
#   geom_linerange(aes(x = Population, ymin = lower, ymax = upper, group = metric, color = metric), size = 2) +
#   geom_point(aes(x = Population, y = mean, group = metric, color = metric), size = 5.5, shape = 21, stroke = 1.5) +
#   scale_shape(solid = FALSE) +
#   scale_color_manual(values = c("#3498DB", "#DC7633", "#6c757d"), 
#                      labels = c( (expression(paste("\u03C0 (", italic("P. polyodon"), ")"))),
#                                  (expression(paste("\u03C0 (", italic("P. kazumbe"), ")"))),
#                                  (("dxy")) )) +
#   #ylim(0, max(pi_info$upper)) +
#   xlab("Population") + ylab("Avg. seq. div.") +
#   theme_cowplot() +
#   theme(plot.margin = margin(0, 10, 0, 13),
#         panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
#         panel.grid.major.x = element_line(colour = 'gray', size = 0.4),
#         legend.title = element_blank(),
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(size = 22),
#         axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 11, b = 0, l = 0)),
#         axis.text.y = element_text(size = 18),
#         legend.text = element_text(size = 18, margin = margin(r = 22, unit = "pt")),
#         legend.position = "bottom",
#         legend.spacing.x = unit(0.3, 'cm')) 
# #theme(panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
# #      legend.position = "none") +
# #coord_flip()
# 
# 
# Fst_plot_horizontal <- ggplot(data = Fst_dataframe_updated %>% 
#                                 mutate(group = 'a')) +
#   #geom_pointrange(aes(x = population, y = mean, ymin = lower, ymax = upper, group = species, color = species), size = 1.2) +
#   geom_linerange(aes(x = Population, ymin = lower_CI, ymax = upper_CI, color = group), size = 2) +
#   geom_point(aes(x = Population, y = Fst, color = group), size = 5.5, shape = 21, stroke = 1.5) +
#   scale_shape(solid = FALSE) +
#   #ylim(0, max(pi_info$upper)) +
#   xlab("Population") + ylab("Fst") +
#   theme_cowplot() +
#   theme(plot.margin = margin(5.5, 10, 0, 13),
#         panel.grid.major.y = element_line(colour = 'gray', size = 0.4),
#         panel.grid.major.x = element_line(colour = 'gray', size = 0.4),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.text.y = element_text(size = 18),
#         axis.title.y = element_text(size = 22, margin = margin(t = 0, r = 11, b = 0, l = 0)),
#         legend.position = "none",
#         legend.text = element_text(color = "white"),
#         legend.title = element_text(color = "white"),
#         legend.key = element_rect(fill = "white")) +
#   scale_x_discrete(drop = FALSE) +
#   scale_color_manual(values = c('#5e60ce'),
#                      guide = guide_legend(override.aes = list(color = "white")))
# #scale_color_discrete(guide = guide_legend(override.aes = list(color = "white"))) +
# #coord_flip()
# 
# diversity_multipanel_horizontal <- ggarrange(Fst_plot_horizontal, divergence_diversity_plot_horizontal,
#                                              nrow = 2)
# 
# 
# 
# map_fish_multi <- plot_grid(tanganyika_sampling_map_final, fish_multipanel, ncol = 2, rel_widths = c(0.5, 0.5), labels = c('(a)', '(b)'), label_x = c(0, -0.03), label_size = 28)
# 
# plot_grid(map_fish_multi, ggdraw(diversity_multipanel_horizontal) + theme(plot.margin = unit(c(0, 0.5, 0.5, 0.8),"cm")), 
#           nrow = 2, rel_heights = c(0.65, 0.35), labels = c('', '(c)'), label_size = 28, align = "hv")
# 
# ggsave2('/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/map_divergencestats_fish_multipanel.png',
#         width = 18, height = 15)
# 
# 
# 
# map_fish_multi_alt <- plot_grid(fish_multipanel, tanganyika_sampling_map_final, ncol = 2, rel_widths = c(0.52, 0.48), labels = c('(a)', '(b)'), label_x = c(0, -0.02), label_size = 28)
# 
# plot_grid(map_fish_multi_alt, ggdraw(diversity_multipanel_horizontal) + theme(plot.margin = unit(c(0, 0.5, 0.5, 0.8),"cm")), 
#           nrow = 2, rel_heights = c(0.65, 0.35), labels = c('', '(c)'), label_size = 28)
# 
# ggsave2('/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/map_divergencestats_fish_multipanel_alt.png',
#         width = 18, height = 15)

