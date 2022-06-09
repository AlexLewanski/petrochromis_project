##############################################################################################
### SCRIPT NAME: popvae_vis.R
### PURPOSE: visualizing the popvae results, which is included as a supplementary figure
### PRODUCTS:
###     petro_multipanelpopvae_9_7_2021.png: Multipanel figure of of plots visualizing the
###                                          popvae results
##############################################################################################


#####################
### SCRIPT SET-UP ###
#####################

### Loading libraries ###
library(tidyverse)
library(cowplot)
library(rhdf5)
library(abind)
library(RColorBrewer)
library(here)


### Loading data ###
cichlid_metadata <- read.csv(here('working_metadata_file', 'monster_23jul20.csv'), stringsAsFactors = FALSE) #metadata
entropy_samples_order <- read.delim(here('working_metadata_file', 'SAMPLES_biallel_polykaz_minq20_maf1_minDP5_maxDP75_miss0.70.recode.txt'), header = FALSE)

popvae_upload <- read.table(here('popvae', 'results', 'results_summer2021', 'poly_kaz_fulldataset_latent_coords.txt'), header = TRUE)
popvae_kazumbe_upload <- read.table(here('popvae', 'results', 'results_summer2021', 'kazumbe_latent_coords.txt'), header = TRUE)
popvae_polyodon_upload <- read.table(here('popvae', 'results', 'results_summer2021', 'polyodon_latent_coords.txt'), header = TRUE)

#K = 2
mod_list <- list()
for (chain in 1:3) {
  mod_list[[paste0('chain_', chain)]] <- h5read(here('entropy', 'hdf5_files', paste0('kazumbe_polyodon_qmod_6_18_2021_', 2, 'rep', chain, '.hdf5')), "q")
}



#######################
### PROCESSING DATA ###
#######################

#changing column names of popvae dataframes
colnames(popvae_upload)[1:2] <- c("LD1","LD2")
colnames(popvae_kazumbe_upload)[1:2] <- c("LD1","LD2")
colnames(popvae_polyodon_upload)[1:2] <- c("LD1","LD2")


#metadata processing
metadata_subset_entropy <- cichlid_metadata[cichlid_metadata$sample_code %in% entropy_samples_order$V1, c('sample_code', 'location', 'sciname2', 'long', 'lat', 'species2')]
metadata_subset_entropy_reorder <- metadata_subset_entropy[match(as.character(entropy_samples_order$V1), metadata_subset_entropy$sample_code),]


### Processing entropy results and combining with metadata ###
k2_combined.q <- abind(mod_list$chain_1, mod_list$chain_2, mod_list$chain_3, along=1)
k2_q_estimate <- as.data.frame(t(apply(k2_combined.q, 2:3, mean)))

metadata_subset_entropy <- cichlid_metadata[cichlid_metadata$sample_code %in% entropy_samples_order$V1, c('sample_code', 'location', 'sciname2', 'long', 'lat', 'species2')]
metadata_subset_entropy_reorder <- metadata_subset_entropy[match(as.character(entropy_samples_order$V1), metadata_subset_entropy$sample_code),]

k2_q_estimate_with_metadata <- cbind(k2_q_estimate, metadata_subset_entropy_reorder)
#V2 greater than 0.5 --> majority kazumbe ancestry
#V2 less than 0.5 --> majority kazumbe ancestry
#there are no individuals with exactly 50/50 kazume/polyodon ancestry
k2_q_estimate_with_metadata$species_ID_entropy[k2_q_estimate_with_metadata$V2 > 0.5] <- 'P. kazumbe'
k2_q_estimate_with_metadata$species_ID_entropy[k2_q_estimate_with_metadata$V2 < 0.5] <- 'P. polyodon'

#adding location names for plotting
k2_q_estimate_with_metadata$location_plotting[k2_q_estimate_with_metadata$location == 'Gombe_South'] <- "Gombe S"
k2_q_estimate_with_metadata$location_plotting[k2_q_estimate_with_metadata$location == 'Katongwe_N'] <- "Katongwe N"
k2_q_estimate_with_metadata$location_plotting[k2_q_estimate_with_metadata$location == 'Katongwe_S'] <- "Katongwe S"
k2_q_estimate_with_metadata$location_plotting[k2_q_estimate_with_metadata$location == 'Ska'] <- "S Kagango"
k2_q_estimate_with_metadata$location_plotting[k2_q_estimate_with_metadata$location == 'Kalala'] <- "Kalalangabo"
k2_q_estimate_with_metadata$location_plotting[k2_q_estimate_with_metadata$location == 'Nondwa'] <- "Nondwa"
k2_q_estimate_with_metadata$location_plotting[k2_q_estimate_with_metadata$location == 'Hilltop'] <- "Hilltop"
k2_q_estimate_with_metadata$location_plotting[k2_q_estimate_with_metadata$location == 'Bangwe'] <- "Bangwe"
k2_q_estimate_with_metadata$location_plotting[k2_q_estimate_with_metadata$location == 'Jakob'] <- "Jakobsen's S"
k2_q_estimate_with_metadata$location_plotting[k2_q_estimate_with_metadata$location == 'Ulombola'] <- "Ulombola"
#k2_q_estimate_with_metadata$location_plotting[k2_q_estimate_with_metadata$location == 'Harembe'] <- "Harembe"

k2_q_estimate_with_metadata$location_plotting <- factor(k2_q_estimate_with_metadata$location_plotting, 
                                                        levels= c("Gombe S", "Katongwe N", "Katongwe S", "S Kagango", "Kalalangabo", "Nondwa", "Hilltop", "Bangwe", "Jakobsen's S", "Ulombola"))


### Creating qualitative color palette ###
location_palette <- brewer.pal(length(unique(k2_q_estimate_with_metadata$location_plotting)),"Set3")
names(location_palette) <- unique(k2_q_estimate_with_metadata$location_plotting)
location_palette[location_palette == '#FFFFB3'] <- '#F7DC6F' #switch the yellow out
location_palette[location_palette == "#D9D9D9"] <- "#AAB7B8" #switch the gray out
location_palette[names(location_palette) == "Hilltop"] <- "#996633"
location_palette[names(location_palette) == "Ulombola"] <- "#0000e6"
location_palette[names(location_palette) == "Nondwa"] <- "#4dd2ff"
location_palette[names(location_palette) == "Gombe S"] <- "#66ff66"
#location_palette[names(location_palette) == "Harembe"] <- "#99ffeb"

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

location_palette_df <- stack(location_palette) %>% 
  mutate(letter = case_when(ind == 'Gombe S' ~ 'A',
                            ind == 'Katongwe N' ~ 'B',
                            ind == 'Katongwe S' ~ 'C',
                            ind == 'S Kagango' ~ 'D',
                            ind == 'Kalalangabo' ~ 'E',
                            ind == 'Nondwa' ~ 'F',
                            ind == 'Hilltop' ~ 'G',
                            ind == 'Bangwe' ~ 'H',
                            ind == "Jakobsen's S" ~ 'I',
                            ind == 'Ulombola' ~ 'J'))

location_palette_letter <- setNames(location_palette_df$values, location_palette_df$letter)


### merge metadata/entropy results dataframe with popvae dataframe ###
popvae_metadata <- merge(popvae_upload, k2_q_estimate_with_metadata, 
                         by.x = 'sampleID', by.y = 'sample_code', all.x = TRUE, all.y = FALSE)

popvae_kazumbe_metadata <- merge(popvae_kazumbe_upload, k2_q_estimate_with_metadata, 
                                 by.x = 'sampleID', by.y = 'sample_code', all.x = TRUE, all.y = FALSE)

popvae_polyodon_metadata <- merge(popvae_polyodon_upload, k2_q_estimate_with_metadata, 
                                  by.x = 'sampleID', by.y = 'sample_code', all.x = TRUE, all.y = FALSE)



##################################
### VISUALIZING POPVAE RESULTS ###
##################################

### Popvae plots ###
#plot of results from popvae inlcuding both kazumbe and polyodon
species_popvae_plot <- ggplot(data = popvae_metadata, aes(x = LD1, y = LD2, color = species_ID_entropy)) +
  geom_point(size = 6, alpha = 0.5) +
  geom_point(size = 6, shape = 21, alpha = 0.7, stroke = 1) +
  scale_color_manual(name = "Species", 
                     values = c("#DC7633", "#3498DB"),
                     labels = c(expression(paste(italic("P"), ". sp. 'kazumbe'", sep = "")), 
                                expression(paste(italic("P"), ". cf. ",  italic("polyodon"), sep = "")))) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 13),
        legend.title = element_text(size = 21),
        legend.position = "bottom",
        legend.spacing.x = unit(0.3, 'cm'),
        legend.text = element_text(size = 18, face = "italic", margin = margin(r = 12, unit = "pt")))
  # theme(panel.grid.major = element_blank(),
  #       panel.grid.minor = element_blank(),
  #       axis.title = element_text(size = 22),
  #       axis.text = element_text(size = 13),
  #       legend.title = element_text(size = 21),
  #       legend.position="bottom",
  #       legend.text = element_text(size = 18, face = "italic"))

#plot of kazumbe popvae results
kazumbe_popvae_plot <- ggplot(data = popvae_kazumbe_metadata, aes(x = LD1, y= LD2, color = location_plotting)) +
  geom_point(size = 4.5, alpha = 0.5) +
  geom_point(size = 4.5, shape = 1, alpha = 0.5) +
  scale_color_manual(name = "Location", values = location_palette[names(location_palette) %in% unique(popvae_kazumbe_metadata$location_plotting) ]) +
  #scale_color_manual(name = "Location", values = gradient_palette[names(gradient_palette) %in% unique(polyodon_PCA_combinedinfo$location_plotting) ]) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 13),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 27, face= c("italic") )) +
  #ggtitle("P. kazumbe")
  ggtitle(expression(paste(italic("Petrochromis"), " sp. 'kazumbe'", sep = "")))

#\textit{Petrochromis} sp. `kazumbe'
#\textit{P}. cf. \textit{polyodon}

#plot of polyodon popvae results
polyodon_popvae_plot <- ggplot(data = popvae_polyodon_metadata, aes(x = LD1, y= LD2, color = location_plotting)) +
  geom_point(size = 4.5, alpha = 0.5) +
  geom_point(size = 4.5, shape = 1, alpha = 0.5) +
  scale_color_manual(name = "Location", values = location_palette[names(location_palette) %in% unique(popvae_polyodon_metadata$location_plotting) ]) +
  #scale_color_manual(name = "Location", values = gradient_palette[names(gradient_palette) %in% unique(polyodon_PCA_combinedinfo$location_plotting) ]) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = 22),
        axis.text = element_text(size = 13),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 27, face= c("italic") )) +
  #ggtitle("P. polyodon")
  ggtitle(expression(paste(italic("Petrochromis"), " cf. ",  italic("polyodon"), sep = "")))


### Create legend ###
popvae_metadata_updated <- popvae_metadata %>% 
  mutate(location_updated = factor(case_when(location == 'Gombe_South' ~ 'A',
                                     location == 'Katongwe_N' ~ 'B',
                                     location == 'Katongwe_S' ~ 'C',
                                     location == 'Ska' ~ 'D',
                                     location == 'Kalala' ~ 'E',
                                     location == 'Nondwa' ~ 'F',
                                     location == 'Hilltop' ~ 'G',
                                     location == 'Bangwe' ~ 'H',
                                     location == 'Jakob' ~ 'I',
                                     location == 'Ulombola' ~ 'J'),
                           levels = LETTERS[1:10]))

location_palette_color_final <- location_palette_letter[names(location_palette_letter) %in% unique(popvae_metadata_updated$location_updated) ]
location_palette_color_final <- location_palette_color_final[match(levels(popvae_metadata_updated$location_updated), names(location_palette_color_final))]

location_legend <- get_legend(ggplot(data = popvae_metadata_updated, 
                                     aes(x = LD1, y = LD2, color = location_updated)) +
                                geom_point(size = 4.5, alpha = 0.5) +
                                geom_point(size = 4.5, shape = 1, alpha = 0.5) +
                                scale_color_manual(name = "Location", values = location_palette_color_final) +
                                #scale_color_manual(name = "Location", values = gradient_palette[names(gradient_palette) %in% unique(petro_PCA_combinedinfo$location_plotting) ]) +
                                theme_bw() +
                                theme(panel.grid.major = element_blank(),
                                      panel.grid.minor = element_blank(),
                                      axis.title = element_text(size = 18),
                                      plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
                                      legend.title = element_text(size = 21),
                                      legend.text=element_text(size = 18)))


### Create the multipanel plot ###
multipanel_kazumbe_polyodon_popvae <- plot_grid(kazumbe_popvae_plot, polyodon_popvae_plot, nrow = 2, labels = c("(b)", "(c)"), label_size = 26.5)
multipanel_kazumbe_polyodon_popvae_withlegend <- plot_grid(multipanel_kazumbe_polyodon_popvae, 
                                                           location_legend, ncol = 2, rel_widths = c(0.9, 0.25))

plot_grid(species_popvae_plot + theme(plot.margin = unit(c(0.75, 0.5, 0, 0.75), "cm")), 
          multipanel_kazumbe_polyodon_popvae_withlegend, 
          rel_widths = c(0.6, 0.4), labels = c('(a)', ''), label_size = 26.5)

#Export the final plot
ggsave2(here('figures', 'petro_multipanelpopvae_9_7_2021.png'), width = 22, height = 12, bg = "white")



#################################
### CODE NOT CURRENTLY IN USE ###
#################################

# ggplot() +
#   geom_point(data = popvae_metadata[which(popvae_metadata$species2 == 'kazumbe'),], 
#              size = 6, alpha = 0.5, aes(x = LD1, y = LD2)) +
#   geom_point(data = popvae_metadata[which(popvae_metadata$species2 == 'kazumbe'),], 
#              size = 6, shape = 21, alpha = 0.7, stroke = 1, aes(x = LD1, y = LD2)) +
#   geom_point(data = popvae_metadata[which(popvae_metadata$species2 == 'polyodon'),], 
#              size = 6, alpha = 0.5, aes(x = LD1, y = LD2)) +
#   geom_point(data = popvae_metadata[which(popvae_metadata$species2 == 'polyodon'),], 
#              size = 6, shape = 21, alpha = 0.7, stroke = 1, aes(x = LD1, y = LD2)) +
#   scale_color_manual(name = "Species", values = c("#DC7633", "#3498DB")) +
#   theme_bw() +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.title = element_text(size = 22),
#         axis.text = element_text(size = 13),
#         legend.title = element_text(size = 21),
#         legend.position="bottom",
#         legend.text=element_text(size = 18, face = "italic"))

# out_path <- '/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/species_lists_and_metadata/'
# 
# kazumbe_ids <- k2_q_estimate_with_metadata[which(k2_q_estimate_with_metadata$species_ID_entropy == 'P. kazumbe'),]$sample_code
# 
# polyodon_ids <- k2_q_estimate_with_metadata[which(k2_q_estimate_with_metadata$species_ID_entropy == 'P. polyodon'),]$sample_code
# 
# write.table(kazumbe_ids, paste0(out_path, 'KAZUMBE_IDs_kazumbe_polyodon_pundamilia_8_19_2020_maf0.01_maxmissing0.70_minDP5_post_4.recode.txt'), 
#             sep="\t", col.names = FALSE, quote = FALSE, row.names = FALSE)
# 
# write.table(polyodon_ids, paste0(out_path, 'POLYODON_IDs_kazumbe_polyodon_pundamilia_8_19_2020_maf0.01_maxmissing0.70_minDP5_post_4.recode.txt'), 
#             sep="\t", col.names = FALSE, quote = FALSE, row.names = FALSE)

# library(readxl)
# sample_info_2007 <- read_excel('/Users/alexlewanski/Downloads/MASTEREXTRACTIONLIST2007.xls', sheet = "All")
# 
# popvae_kazumbe_with_sexinfo <- sample_info_2007 %>% 
#   mutate(sampleID = gsub("\\.", "_", `Specimen #`),
#          sex = if_else(is.na(sex), 'unknown', sex)) %>% 
#   select(sampleID, sex) %>% 
#   left_join(popvae_kazumbe_metadata, ., by = 'sampleID') %>% 
#   mutate(sex = if_else(is.na(sex), 'unknown', sex))
# 
# sexinfo_samples <- sample_info_2007 %>% 
#   mutate(sampleID = gsub("\\.", "_", `Specimen #`),
#          sex = if_else(is.na(sex), 'unknown', sex)) %>% 
#   filter(Species == 'Petrochromis kazumbe' & sex != 'unknown') %>% 
#   pull(sampleID)
# 
# sexinfo_samples[sexinfo_samples %in% popvae_kazumbe_with_sexinfo$sampleID]
# 
# table(popvae_kazumbe_with_sexinfo$sex)
# 
# 
# popvae_kazumbe_with_sexinfo %>% 
#   ggplot() +
#   geom_point(aes(x = LD1, y = LD2, color = sex, size = sex, alpha = sex)) +
#   scale_size_manual(values = c(5, 5, 3)) +
#   scale_color_manual(values = c('#e56b6f', '#00afb9', 'gray')) +
#   scale_alpha_manual(values = c(0.8, 0.8, 0.5)) +
#   theme_bw()

