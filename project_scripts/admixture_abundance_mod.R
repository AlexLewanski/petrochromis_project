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
library(rhdf5)
library(abind)
library(tidyverse)
library(cowplot)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(mgcv)
library(gratia)
library(kableExtra)
library(here)


### Loading data ###
cichlid_metadata <- read.csv(here('working_metadata_file', 'monster_23jul20.csv'), stringsAsFactors = FALSE) #metadata
entropy_samples_order <- read.delim(here('working_metadata_file', 'SAMPLES_biallel_polykaz_minq20_maf1_minDP5_maxDP75_miss0.70.recode.txt'), header = FALSE) #order of samples in entropy results
polyodon_kazumbe_survey_data <- read.csv(here('non_vcf_datafiles', 'ppol_pkaz.csv'), stringsAsFactors = FALSE) #survey data

#K = 2
mod_list <- list()
for (chain in 1:3) {
  mod_list[[paste0('chain_', chain)]] <- h5read(here('entropy', 'hdf5_files', paste0('kazumbe_polyodon_qmod_6_18_2021_', 2, 'rep', chain, '.hdf5')), "q")
}

k2_combined.q <- abind(mod_list$chain_1, mod_list$chain_2, mod_list$chain_3, along=1)
k2_q_estimate <- as.data.frame(t(apply(k2_combined.q, 2:3, mean)))


### Custom functions ###
construct_gam_confint <- function(mod, number_predicted_vals = 50) {
  #https://stats.stackexchange.com/questions/33327/confidence-interval-for-gam-model
  
  predictor_var <- colnames(mod$model)[2]
  
  new_data_df <- data.frame(x = seq(min(mod$model[, 2]), max(mod$model[, 2]), length.out = number_predicted_vals)) 
  colnames(new_data_df) <- predictor_var
  
  #point estimate
  estimate_df <- data.frame(predictor = new_data_df[, 1],
                            estimate = predict(mod, new_data_df, type = "response", se.fit = FALSE))
  
  #confidence intervals
  mod_predict <- predict(mod, new_data_df, type = "link", se.fit = TRUE)
  
  lower_confint <- mod_predict$fit - (2 * mod_predict$se.fit)
  upper_confint <- mod_predict$fit + (2 * mod_predict$se.fit)
  
  estimate_df$lower_confint <- as.vector(mod$family$linkinv(lower_confint))
  estimate_df$upper_confint <- as.vector(mod$family$linkinv(upper_confint))
  
  return(estimate_df)
}

theme_Q_entropy <- function (base_size = 12) { 
  theme(plot.background = element_rect(colour = NA),
        axis.ticks = element_line(colour = "black", size = 1.25),
        axis.ticks.length=unit(1.2, "mm"),
        #panel.grid.major = element_line(colour = "#E8E8E8", size = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        axis.line = element_line(colour = "black", size = 1.4, lineend = 'square'),
        axis.text=element_text(size=rel(1.5)),
        #axis.title.x = element_text(face = "bold", size=rel(1.95), margin = unit(c(3, 0, 0, 0), "mm")),
        #axis.title.y = element_text(face = "bold", size=rel(1.95), margin = unit(c(0, 3, 0, 0), "mm")),
        axis.text.x=element_text(margin = margin(t = 6)),
        axis.text.y=element_text(margin = margin(r = 6)),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        legend.key = element_rect(fill = "white"),
        legend.key.size= unit(8, "mm"),
        #legend.title = element_text(size = rel(1.5)),
        legend.title = element_text(size = 22),
        legend.text=element_text(size=rel(1.5)))
}

### Script resources ###
#https://stackoverflow.com/questions/32555531/how-to-italicize-part-one-or-two-words-of-an-axis-title
#https://stats.stackexchange.com/questions/33327/confidence-interval-for-gam-model
#https://stats.stackexchange.com/questions/356067/which-method-is-correct-generalized-additive-model-mgcv
#https://fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/
#https://stackoverflow.com/questions/39282293/r-ggplot2-using-italics-and-non-italics-in-the-same-category-label (for italicizing part of axis text)


#######################
### PROCESSING DATA ###
#######################

### Metadata processing ###
#updating location names
cichlid_metadata <- cichlid_metadata %>%
  mutate(location_updated = case_when(location == 'Ska' ~ 'Skagongo',
                                      location == 'Kalala' ~ 'Kalalangabo',
                                      !(location %in% c('Ska', 'Kalala')) ~ location))

metadata_subset <- cichlid_metadata[cichlid_metadata$sample_code %in% entropy_samples_order$V1, c('sample_code', 'location', 'sciname2', 'long', 'lat', 'species2', 'location_updated')]
metadata_subset_reorder <- metadata_subset[match(as.character(entropy_samples_order$V1), metadata_subset$sample_code),]

#add a column of "nice" location names to be plotted
metadata_subset_reorder$location_plotting <- NA
metadata_subset_reorder$location_plotting[metadata_subset_reorder$location == 'Gombe_South'] <- "Gombe S"
metadata_subset_reorder$location_plotting[metadata_subset_reorder$location == 'Katongwe_N'] <- "Katongwe N"
metadata_subset_reorder$location_plotting[metadata_subset_reorder$location == 'Katongwe_S'] <- "Katongwe S"
metadata_subset_reorder$location_plotting[metadata_subset_reorder$location == 'Ska'] <- "S Kagango"
metadata_subset_reorder$location_plotting[metadata_subset_reorder$location == 'Kalala'] <- "Kalalangabo"
metadata_subset_reorder$location_plotting[metadata_subset_reorder$location == 'Nondwa'] <- "Nondwa"
metadata_subset_reorder$location_plotting[metadata_subset_reorder$location == 'Hilltop'] <- "Hilltop"
metadata_subset_reorder$location_plotting[metadata_subset_reorder$location == 'Bangwe'] <- "Bangwe"
metadata_subset_reorder$location_plotting[metadata_subset_reorder$location == 'Jakob'] <- "Jakobsen's S"
metadata_subset_reorder$location_plotting[metadata_subset_reorder$location == 'Ulombola'] <- "Ulombola"
metadata_subset_reorder$location_plotting[metadata_subset_reorder$location == 'Harembe'] <- "Harembe"
metadata_subset_reorder$location_plotting <- factor(metadata_subset_reorder$location_plotting, 
                                                    levels= c("Gombe S", "Katongwe N", "Katongwe S", "S Kagango", "Kalalangabo", "Nondwa", "Hilltop", "Bangwe", "Jakobsen's S", "Ulombola", "Harembe"))

#add a column of abbreviated names for kazumbe and polyodon
metadata_subset_reorder$reduced_species_names <- ifelse(metadata_subset_reorder$sciname2 == "Petrochromis kazumbe",
                                                        "P. kazumbe",
                                                        "P. polyodon")

### survey data processing ###
#summarizing survey data by species and location
SUMMARIZED_polyodon_kazumbe_survey_data <- polyodon_kazumbe_survey_data %>%
  group_by(location, abbreviation) %>%
  summarise(total_sample_area = sum(samp_area), n_total = sum(n)) %>% 
  ungroup() %>% 
  mutate(total_density = n_total/total_sample_area,
         species_location_id = paste(abbreviation, location, sep = "_")) %>% 
  as.data.frame()

#summarizing survey data by species, year and location and then summarizing by location and species
SUMMARIZED_BY_YEAR_polyodon_kazumbe_survey_data <- polyodon_kazumbe_survey_data %>%
  mutate(survey_year = paste0("20", gsub(" \\d*:\\d*", "",  gsub("\\d*/\\d*/" , "", date_time) ) )) %>% 
  group_by(location, abbreviation, survey_year) %>%
  summarise(total_sample_area_year = sum(samp_area), n_total_per_year = sum(n)) %>% 
  ungroup() %>%
  mutate(density_per_year = n_total_per_year/total_sample_area_year,
         species_location_id = paste(abbreviation, location , sep = "_")) %>%
  group_by(species_location_id) %>% 
  summarise(mean_density = mean(density_per_year),
            location = first(location),
            abbreviation = first(abbreviation)) %>% 
  ungroup() %>% 
  as.data.frame()


#merge the two survey dataframes together
SUMMARIZED_polyodon_kazumbe_survey_data_FULL <- merge(SUMMARIZED_polyodon_kazumbe_survey_data, 
                                                      SUMMARIZED_BY_YEAR_polyodon_kazumbe_survey_data[,c('species_location_id', 'mean_density')], by = "species_location_id")


### Color palettes ###
location_palette <- setNames(brewer.pal(11,"Set3"), c("Hilltop", "Jakobsen's S", "Nondwa", "Kalalangabo", "Ulombola", "Katongwe N", "Gombe S", "Katongwe S", "Bangwe", "S Kagango", "Harembe"))
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


k2_q_estimate_with_metadata <- cbind(k2_q_estimate, metadata_subset_reorder)

k2_q_long <- k2_q_estimate_with_metadata %>%
  gather(key = cluster, value = proportion, c('V1', 'V2'))



### Creating dataset used in analysis ###
#POLYODON
polyodon_individuals <- k2_q_estimate_with_metadata[k2_q_estimate_with_metadata$V1 > 0.5,] #extracting polyodon individuals based on entropy ancestry estimates
merged_dataset_polyodon <- merge(polyodon_individuals, 
                                 SUMMARIZED_polyodon_kazumbe_survey_data_FULL[SUMMARIZED_polyodon_kazumbe_survey_data_FULL$abbreviation == 'Ppolyodon', c('location', 'total_density', 'mean_density')], 
                                 by.x = 'location_updated', by.y = 'location', all.x = FALSE, all.y = FALSE)
colnames(merged_dataset_polyodon)[colnames(merged_dataset_polyodon) == "total_density"] <- "total_density_polyodon"
colnames(merged_dataset_polyodon)[colnames(merged_dataset_polyodon) == "mean_density"] <- "mean_density_polyodon"

merged_dataset_polyodon <- merge(merged_dataset_polyodon, 
                                 SUMMARIZED_polyodon_kazumbe_survey_data_FULL[SUMMARIZED_polyodon_kazumbe_survey_data_FULL$abbreviation == 'Pkazumbe', c('location', 'total_density', 'mean_density')], 
                                 by.x = 'location_updated', by.y = 'location', all.x = FALSE, all.y = FALSE)
colnames(merged_dataset_polyodon)[colnames(merged_dataset_polyodon) == "total_density"] <- "total_density_kazumbe"
colnames(merged_dataset_polyodon)[colnames(merged_dataset_polyodon) == "mean_density"] <- "mean_density_kazumbe"


#KAZUMBE
kazumbe_individuals <- k2_q_estimate_with_metadata[which(k2_q_estimate_with_metadata$V2 > 0.5),]
merged_dataset_kazumbe <- merge(kazumbe_individuals, SUMMARIZED_polyodon_kazumbe_survey_data_FULL[which(SUMMARIZED_polyodon_kazumbe_survey_data_FULL$abbreviation == 'Pkazumbe'), c('location', 'total_density', 'mean_density')], 
                                by.x = 'location_updated', by.y = 'location', all.x = FALSE, all.y = FALSE)
#colnames(merged_dataset_kazumbe)[colnames(merged_dataset_kazumbe) == "density"] <- "density_kazumbe"
colnames(merged_dataset_kazumbe)[colnames(merged_dataset_kazumbe) == "total_density"] <- "total_density_kazumbe"
colnames(merged_dataset_kazumbe)[colnames(merged_dataset_kazumbe) == "mean_density"] <- "mean_density_kazumbe"

merged_dataset_kazumbe <- merge(merged_dataset_kazumbe, SUMMARIZED_polyodon_kazumbe_survey_data_FULL[which(SUMMARIZED_polyodon_kazumbe_survey_data_FULL$abbreviation == 'Ppolyodon'), c('location', 'total_density', 'mean_density')], 
                                by.x = 'location_updated', by.y = 'location', all.x = FALSE, all.y = FALSE)
#colnames(merged_dataset_kazumbe)[colnames(merged_dataset_kazumbe) == "density"] <- "density_polyodon"
colnames(merged_dataset_kazumbe)[colnames(merged_dataset_kazumbe) == "total_density"] <- "total_density_polyodon"
colnames(merged_dataset_kazumbe)[colnames(merged_dataset_kazumbe) == "mean_density"] <- "mean_density_polyodon"



#############################################
### GAM MODELS OF ADMIXTURE VS. ABUNDANCE ###
#############################################
v1_kazumbedens_gam_beta <- gam(V1 ~ s(mean_density_kazumbe, k = 4, bs = "tp"),
                               family = betar(link = "logit"),
                               data = merged_dataset_kazumbe, 
                               method = "REML")

v1_polyodondens_gam_beta <- gam(V2 ~ s(mean_density_polyodon, k = 4, bs = "tp"),
                                family = betar(link = "logit"),
                                data = merged_dataset_polyodon, 
                                method = "REML")



############################
### GAM MAIN TEXT FIGURE ###
############################

#create dataframe of predictor value, estimate, and lower and upper CIs for each GAM
kazumbe_mod_info <- construct_gam_confint(mod = v1_kazumbedens_gam_beta)
polyodon_mod_info <- construct_gam_confint(mod = v1_polyodondens_gam_beta)

#polydon dataframe used in plotting
polyodon_plotting_df <- merged_dataset_polyodon %>%
  mutate(Location = case_when(location == 'Gombe_South' ~ "A",
                              location == 'Katongwe_N' ~ "B",
                              location == 'Katongwe_S' ~ "C",
                              location == 'Ska' ~ "D",
                              location == 'Kalala' ~ "E",
                              location == 'Nondwa' ~ "F",
                              location == 'Hilltop' ~ "G",
                              location == 'Bangwe' ~ "H",
                              location == 'Jakob' ~ "I",
                              location == 'Ulombola' ~ "J"))

#visualizaton of polyodon GAM
gam_polyodon_plot <- ggplot(data = polyodon_mod_info, aes(predictor, estimate)) +
  geom_point(data = polyodon_plotting_df, mapping = aes(mean_density_polyodon, V2, color = Location), size = 6.5) +
  geom_ribbon(aes(predictor, ymin = lower_confint, ymax = upper_confint), color = NA, fill = 'gray', alpha = 0.65) +
  geom_line(aes(predictor, estimate), size = 3.5, lineend = "round", color = 'black', linetype = "dashed") +
  scale_color_manual(name = "Location", values = location_palette_letter[names(location_palette_letter) %in% unique(polyodon_plotting_df$Location) ]) +
  #geom_smooth_ci(size = 2, lineend = "round", color = 'black') +
  #xlab(expression(italic('P. polyodon') ~ 'abundance (standardized)' )) +
  xlab(expression(paste(italic("P"), ". cf. ",  italic("polyodon"), ' abundance (standardized)', sep = ""))) +
  #ylab(expression(italic('P. kazumbe') ~ 'ancestry' )) +
  ylab(expression(paste(italic("P"), ". sp. 'kazumbe' ancestry", sep = ""))) +
  #ggtitle(expression(italic('P. polyodon'))) +
  ggtitle(expression(paste(italic("Petrochromis"), " cf. ",  italic("polyodon"), sep = ""))) +
  theme_cowplot() +
  theme(legend.position = "none",
        plot.margin=unit(c(0.8, 0.25, 0, 0.25),"cm"),
        axis.line = element_line(size = 1.5),
        axis.ticks.length=unit(.15, "cm"),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 19),
        plot.title = element_text(hjust = 0.5, size = 26))

#kazumbe dataframe used in plotting
kazumbe_plotting_df <- merged_dataset_kazumbe %>%
  mutate(Location = case_when(location == 'Gombe_South' ~ "A",
                                     location == 'Katongwe_N' ~ "B",
                                     location == 'Katongwe_S' ~ "C",
                                     location == 'Ska' ~ "D",
                                     location == 'Kalala' ~ "E",
                                     location == 'Nondwa' ~ "F",
                                     location == 'Hilltop' ~ "G",
                                     location == 'Bangwe' ~ "H",
                                     location == 'Jakob' ~ "I",
                                     location == 'Ulombola' ~ "J"))

#visualizaton of kazumbe GAM
gam_kazumbe_plot <- ggplot(data = kazumbe_mod_info, aes(predictor, estimate)) +
  geom_point(data = kazumbe_plotting_df,
             mapping = aes(mean_density_kazumbe, V1, color = Location), size = 6.5) +
  geom_ribbon(aes(predictor, ymin = lower_confint, ymax = upper_confint), fill = 'gray', color = NA, size = 1.5, alpha = 0.65) +
  geom_line(aes(predictor, estimate), size = 3.5, lineend = "round", color = 'black', linetype = "dashed") +
  scale_color_manual(name = "Location", values = location_palette_letter[names(location_palette_letter) %in% unique(kazumbe_plotting_df$Location) ]) +
  #geom_smooth_ci(size = 2, lineend = "round", color = 'black') +
  #xlab(expression(italic('P. kazumbe') ~ 'abundance (standardized)' )) +
  xlab(expression(paste(italic("P"), ". sp. 'kazumbe' abundance (standardized)", sep = ""))) +
  #ylab(expression(italic('P. polyodon') ~ 'ancestry' )) +
  ylab(expression(paste(italic("P"), ". cf. ",  italic("polyodon"), ' ancestry', sep = ""))) +
  #ggtitle(expression(italic('P. kazumbe'))) +
  ggtitle(expression(paste(italic("Petrochromis"), " sp. 'kazumbe'", sep = ""))) +
  theme_cowplot()  +
  theme(legend.position = "none",
        plot.margin=unit(c(0.8, 0.25, 0, 0.25),"cm"),
        axis.line = element_line(size = 1.5),
        axis.ticks.length=unit(.15, "cm"),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 19),
        plot.title = element_text(hjust = 0.5, size = 26))


#processing abundance information for visualization
kaz_pol_abundance_df <- merged_dataset_polyodon %>%
  distinct(location_updated, .keep_all = TRUE) %>%
  select(location_updated, mean_density_polyodon, mean_density_kazumbe) %>%
  rename(polyodon = mean_density_polyodon, kazumbe = mean_density_kazumbe) %>%
  pivot_longer(!location_updated, names_to = "species", values_to = "density") %>%
  mutate(species_plotting = ifelse(species == 'polyodon', "P. polyodon", "P. kazumbe"),
         Location = factor(case_when(location_updated == 'Gombe_South' ~ "A",
                                     location_updated == 'Katongwe_N' ~ "B",
                                     location_updated == 'Katongwe_S' ~ "C",
                                     location_updated == 'Skagongo' ~ "D",
                                     location_updated == 'Kalalangabo' ~ "E",
                                     location_updated == 'Nondwa' ~ "F",
                                     location_updated == 'Hilltop' ~ "G",
                                     location_updated == 'Bangwe' ~ "H",
                                     location_updated == 'Jakob' ~ "I",
                                     location_updated == 'Ulombola' ~ "J"), levels =  LETTERS[10:1]))

#visualizing abundance information
poly_kaz_abundance_plot <- ggplot(data = kaz_pol_abundance_df,
       aes(x = species_plotting, y = density)) +
  geom_line(aes(group = Location, color = Location), size = 3) +
  geom_point(aes(color = Location), size = 6.5) +
  scale_x_discrete(expand = c(0.15, 0.15),
                   labels = c(expression(paste(italic("P"), ". sp. 'kazumbe'", sep = "")), 
                              expression(paste(italic("P"), ". cf. ",  italic("polyodon"), sep = "")))) +
  scale_color_manual(name = "Location", 
                     values = location_palette_letter[names(location_palette_letter) %in% unique(kaz_pol_abundance_df$Location) ]) +
  xlab('Species') + ylab('Standardized abundance') +
  theme_cowplot() +
  theme(legend.position = "none",
        plot.margin = unit(c(0.8, 0.25, 0, 0.25),"cm"),
        axis.line = element_line(size = 1.5),
        axis.ticks.length=unit(.15, "cm"),
        axis.ticks = element_line(colour = "black", size = 1.5),
        axis.text = element_text(size = 16),
        axis.title = element_text(size = 20),
        axis.text.x = element_text(face = 'italic'))

#creating legend for multipanel plot
location_palette_letter_processed <- location_palette_letter[names(location_palette_letter) %in% unique(kaz_pol_abundance_df$Location) ]
location_palette_letter_processed <- location_palette_letter_processed[match(sort(names(location_palette_letter_processed)), names(location_palette_letter_processed))]
location_legend <- get_legend(ggplot(data = kaz_pol_abundance_df %>%
                                       mutate(Location = factor(Location, levels = c('D', 'E', 'F', 'G', 'I'))),
                                     aes(x = species_plotting, y = density)) +
                                geom_line(aes(group = Location, color = Location), size = 3) +
                                geom_point(aes(color = Location), size = 6.5) +
                                scale_color_manual(name = "Location", values = location_palette_letter_processed) +
                                xlab('Species') + ylab('Standardized abundance') +
                                theme(legend.position = "bottom",
                                      legend.spacing.x = unit(0.3, 'cm'),
                                      legend.background = element_blank(),
                                      legend.key = element_rect(fill = "white"),
                                      legend.title = element_text(size = 20),
                                      legend.text = element_text(size = 18, margin = margin(r = 10, unit = "pt")),
                                      plot.margin = unit(c(0.8, 0.25, 0, 0.25),"cm"),
                                      axis.line = element_line(size = 1.5),
                                      axis.ticks.length=unit(.15, "cm"),
                                      axis.ticks = element_line(colour = "black", size = 1.5),
                                      axis.text = element_text(size = 15),
                                      axis.title = element_text(size = 18),
                                      axis.text.x = element_text(face = 'italic')) +
                                guides(fill = guide_legend(label.position = "bottom")) )


#creating multipanel plot of GAM results and abundance information
plot_grid(plot_grid(gam_kazumbe_plot, gam_polyodon_plot, poly_kaz_abundance_plot,
                    nrow = 1, labels = c('(a)', '(b)', '(c)'), label_size = 27, vjust = 2, align = 'hv', axis = "bt"),
          location_legend,
          nrow = 2, rel_heights = c(1, 0.1))

#saving multipanel figure
ggsave2(here('figures', 'abundance_admixture_gam_9_8_2021.png'), width = 18, height = 6.5, bg = "white")



#############################
### GAM MODEL FIT VISUALS ###
#############################
polyodon_mod_fit <- cowplot::plot_grid(qq_plot(v1_polyodondens_gam_beta) + theme_bw(), residuals_linpred_plot(v1_polyodondens_gam_beta) + theme_bw(),
                                       residuals_hist_plot(v1_polyodondens_gam_beta) + theme_bw(), observed_fitted_plot(v1_polyodondens_gam_beta) + theme_bw(), 
                                       nrow = 2)
ggsave2(here('figures', 'polyodon_gam_mod_fit_9_8_2021.png'), width = 14, height = 10, bg = "white")


kazumbe_mod_fit <- cowplot::plot_grid(qq_plot(v1_kazumbedens_gam_beta) + theme_bw(), residuals_linpred_plot(v1_kazumbedens_gam_beta) + theme_bw(),
                                      residuals_hist_plot(v1_kazumbedens_gam_beta) + theme_bw(), observed_fitted_plot(v1_kazumbedens_gam_beta) + theme_bw(), nrow = 2)

ggsave2(here('figures', 'kazumbe_gam_mod_fit_9_8_2021.png'), width = 14, height = 10, bg = "white")



###########################################################################
### SUMMARY TABLES (GAM MODEL RESULTS AND POLYODON/KAZUMBE SURVEY INFO) ###
###########################################################################

### Model results table ###
polyodon_gam_summary <- summary(v1_polyodondens_gam_beta)
kazumbe_gam_summary <- summary(v1_kazumbedens_gam_beta)

mod_summary <- data.frame("Model formula" = c("polyodon ancestry $\\sim$ s(kazumbe abundance)", "kazumbe ancestry $\\sim$ s(polyodon abundance)"),
                          "$R^{2}_{adj}$" = c(round(kazumbe_gam_summary$r.sq, 2), round(polyodon_gam_summary$r.sq, 2)),
                          "edf" = c(round(kazumbe_gam_summary$edf, 2), round(polyodon_gam_summary$edf, 2)),
                          "chi.sq" = c(round(kazumbe_gam_summary$chi.sq, 2), round(polyodon_gam_summary$chi.sq, 2)),
                          "$P$" = c(ifelse(kazumbe_gam_summary$s.pv < 0.001, "\\textless0.001", round(kazumbe_gam_summary$s.pv, 3)),
                                    ifelse(polyodon_gam_summary$s.pv < 0.001, "\\textless0.001", round(polyodon_gam_summary$s.pv, 3))),
                          check.names = FALSE)
rownames(mod_summary) <- NULL


gam_mod_summary_table <- kbl(mod_summary, "latex", booktabs = TRUE, escape = FALSE, align = "c") %>% 
  add_header_above(c("", "", "Smooth term" = 3)) %>% 
  kable_styling(latex_options = c("hold_position"))

cat(gam_mod_summary_table, file = here('tables', 'gam_summary_table.txt'), append = FALSE)


### Survey information ###
survey_info_table_df <- SUMMARIZED_polyodon_kazumbe_survey_data_FULL %>%
  filter(location %in% k2_q_estimate_with_metadata$location_updated) %>% 
  mutate("Location name" = case_when(location == 'Gombe_South' ~ 'Gombe S',
                                     location == 'Katongwe_N' ~ 'Katongwe N',
                                     location == 'Katongwe_S' ~ 'Katongwe S',
                                     location == 'Skagongo' ~ 'S Kagongo',
                                     location == 'Kalalangabo' ~ 'Kalalangabo',
                                     location == 'Nondwa' ~ 'Nondwa',
                                     location == 'Hilltop' ~ 'Hilltop',
                                     location == 'Bangwe' ~ 'Bangwe',
                                     location == 'Jakob' ~ "Jakobsen's Beach",
                                     location == 'Ulombola' ~ 'Ulombola'),
         "Location letter" = case_when(location == 'Gombe_South' ~ "A",
                                       location == 'Katongwe_N' ~ 'B',
                                       location == 'Katongwe_S' ~ 'C',
                                       location == 'Skagongo' ~ 'D',
                                       location == 'Kalalangabo' ~ 'E',
                                       location == 'Nondwa' ~ 'F',
                                       location == 'Hilltop' ~ 'G',
                                       location == 'Bangwe' ~ 'H',
                                       location == 'Jakob' ~ 'I',
                                       location == 'Ulombola' ~ 'J'),
         abbreviation = case_when(abbreviation == 'Pkazumbe' ~ 'kazumbe',
                                  abbreviation == 'Ppolyodon' ~ 'polyodon'),
         total_density = round(total_density, 4),
         mean_density = round(mean_density, 4)) %>%
  select(abbreviation, `Location letter`, `Location name`, n_total, total_sample_area, total_density, mean_density) %>% 
  rename(Species = abbreviation,
         N = n_total, 
         `Total area (m^2)` = total_sample_area,
         `Density (total)` = total_density,
         `Density (avg per year)` = mean_density) %>% 
  arrange(Species, `Location letter`)

survey_info_table <- kbl(survey_info_table_df, "latex", booktabs = TRUE, align = "c") %>% 
  kable_styling(latex_options = c("scale_down", "hold_position")) %>% 
  collapse_rows(columns = 1, valign = "middle", latex_hline = "major")

cat(survey_info_table, file = here('tables', 'survey_info_table.txt'), append = FALSE)



#################################
### CODE CURRENTLY NOT IN USE ###
#################################

# glm_confint_dataframe <- function(mod, data, predictor) {
#   #creating proper confidence intervals for glms
#   #https://fromthebottomoftheheap.net/2018/12/10/confidence-intervals-for-glms/
#   
#   fam <- family(mod)
#   ilink <- fam$linkinv
#   
#   ndata <- with(data, data_frame(var = seq(min(data[, predictor]), max(data[, predictor]),
#                                            length = 100)))
#   colnames(ndata) <- predictor
#   ndata <- add_column(ndata, fit = predict(mod, newdata = ndata, type = 'response'))
#   ndata <- bind_cols(ndata, setNames(as_tibble(predict(mod, newdata = as.data.frame(ndata), se.fit = TRUE)[1:2]),
#                                      c('fit_link','se_link')))
#   
#   ndata <- mutate(ndata,
#                   fit_resp  = ilink(fit_link),
#                   right_upr = ilink(fit_link + (2 * se_link)),
#                   right_lwr = ilink(fit_link - (2 * se_link)))
#   
#   return(ndata)
#   
# }

