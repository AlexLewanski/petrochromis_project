###################################################################################################
### SCRIPT NAME: entropy_results_visualization.R
### PURPOSE: Process and visualize the entropy model results.
### PRODUCTS:
###     k2_Q_q_multipanel_FULL_barchart_marginal_qbarplot_9_7_2021.png: main text entropy figure
###                                                                     including Q vs. q bivariate
###                                                                     plots, summary of hybrid 
###                                                                     classes, and barplots of q
###     entropy_k2through5.png: barplots of q for K = 2 through K = 5. This plot is included in the
###                             supplementary material
###     entropy_trace_plots.png: examples of trace plots for the K = 2 and K = 3 entropy models.
###                              his plot is included in the supplementary material.
###     The "RESULTS SECTION SUMMARIES" section includes summaries of the entropy analyses that
###     are reported in the text of the Results section.
###################################################################################################


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
library(ggExtra)
library(kableExtra)
library(here)

#load entropy visualization functions
source(here('project_scripts', 'entropy_vis_functions.R'))


#table_save_dir <- '/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/tables/'


### Loading data ###
#paths to data
#hdf5_file_path_q <- '/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/entropy/hdf5_files/'
#hdf5_file_path_Q <- '/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/entropy/hdf5_files/'

#metadata
cichlid_metadata <- read.csv(here('working_metadata_file', 'monster_23jul20.csv'), stringsAsFactors = FALSE) #metadata
entropy_samples_order_q <- read.delim(here('working_metadata_file', 'SAMPLES_biallel_polykaz_minq20_maf1_minDP5_maxDP75_miss0.70.recode.txt'), header = FALSE)
entropy_samples_order_Q_north <- read.delim(here('entropy/sample_lists', 'SAMPLES_north_biallel_polykaz_minq20_maf1_minDP5_maxDP75_miss0.70.txt'), header = FALSE)
entropy_samples_order_Q_mid <- read.delim(here('entropy', 'sample_lists', 'SAMPLES_mid_biallel_polykaz_minq20_maf1_minDP5_maxDP75_miss0.70.txt'), header = FALSE)


#full dataset q results
entropy_upload_list <- lapply(split(2:5, paste0('K_', 2:5)), function(k) {
  mod_list <- list()
  
  for (chain in 1:3) {
    #mod_list[[paste0('chain_', chain)]] <- h5read(paste0(hdf5_file_path_q, paste0('kazumbe_polyodon_qmod_6_18_2021_', k, 'rep', chain, '.hdf5' )), "q")
    mod_list[[paste0('chain_', chain)]] <- h5read(here('entropy', 'hdf5_files', paste0('kazumbe_polyodon_qmod_6_18_2021_', k, 'rep', chain, '.hdf5')), "q")
  }
  
  message('finished the following upload: k = ', k)
  return(mod_list)
} )

#region specific Q/q results
k2_Q_q_upload_list <- list()
for (region in c('north', 'mid')) {
  for (chain in 1:3) {
    #k2_Q_q_upload_list[[paste0(region, '_q')]][[paste0('chain_', chain)]] <- h5read(paste0(hdf5_file_path_Q, paste0(region, '_kazumbe_polyodon_ancestrycomplement_6_18_2021_2rep', chain, '.hdf5' )), "q")
    #k2_Q_q_upload_list[[paste0(region, '_Q')]][[paste0('chain_', chain)]] <- h5read(paste0(hdf5_file_path_Q, paste0(region, '_kazumbe_polyodon_ancestrycomplement_6_18_2021_2rep', chain, '.hdf5' )), "Q")
    
    k2_Q_q_upload_list[[paste0(region, '_q')]][[paste0('chain_', chain)]] <- h5read(here('entropy', 'hdf5_files', paste0(region, '_kazumbe_polyodon_ancestrycomplement_6_18_2021_2rep', chain, '.hdf5' )), "q")
    k2_Q_q_upload_list[[paste0(region, '_Q')]][[paste0('chain_', chain)]] <- h5read(here('entropy', 'hdf5_files', paste0(region, '_kazumbe_polyodon_ancestrycomplement_6_18_2021_2rep', chain, '.hdf5' )), "Q")
  }
}


### ggplot aesthetics ###
#ggplot plotting theme for Q vs. q plots
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

theme_Q_entropy_alt <- function (base_size = 12) { 
  theme(plot.background = element_rect(colour = NA),
        axis.ticks = element_line(colour = "black", size = 1.25),
        axis.ticks.length=unit(1.2, "mm"),
        #panel.grid.major = element_line(colour = "#E8E8E8", size = 0.5), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white", colour = NA),
        axis.line = element_line(colour = "black", size = 1.4, lineend = 'square'),
        axis.text=element_text(size=rel(1)),
        #axis.title.x = element_text(face = "bold", size=rel(1.95), margin = unit(c(3, 0, 0, 0), "mm")),
        #axis.title.y = element_text(face = "bold", size=rel(1.95), margin = unit(c(0, 3, 0, 0), "mm")),
        axis.text.x=element_text(size = rel(1.2), margin = margin(t = 6)),
        axis.text.y=element_text(size = rel(1.2), margin = margin(r = 6)),
        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
        legend.key = element_rect(fill = "white"),
        legend.key.size= unit(5, "mm"),
        #legend.title = element_text(size = rel(1.5)),
        legend.title = element_text(size = 16),
        legend.text=element_text(size=rel(1.2)))
}


### color palette for q plots ###
#col_pal <- c("#80CEEF", "#B5D572", "#E86C4A", "#F5A3F3", "#FCE573", "#F0901A")
q_color_palette <- brewer.pal(12,"Paired")[seq(1, 12, by = 2)] #extract the light colors

#color palette for Q vs. q plots (colored based on location)
location_palette <- brewer.pal(11,"Set3")
#names(location_palette) <- unique(metadata_subset_reorder$location_plotting)
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


location_conversion_df <- data.frame(location = c("Gombe S", "Katongwe N", "Katongwe S", "S Kagango", "Kalalangabo", "Nondwa", "Hilltop", "Bangwe", "Jakobsen's S", "Ulombola"),
                                     location_plotting = LETTERS[1:10])


location_palette_letters <- location_palette
for (i in names(location_palette_letters)) {
  
  if (i %in% location_conversion_df$location) {
    names(location_palette_letters)[names(location_palette_letters) == i] <- location_conversion_df[location_conversion_df$location == i,]$location_plotting
  }
}



#######################
### Processing data ###
#######################

### Processing metadata ###

# location location_plotting
#       Gombe S                 A
#    Katongwe N                 B
#    Katongwe S                 C
#     S Kagango                 D
#   Kalalangabo                 E
#        Nondwa                 F
#       Hilltop                 G
#        Bangwe                 H
#  Jakobsen's S                 I
#      Ulombola                 J

processed_metadata <- lapply(list(full = entropy_samples_order_q,
                                  north = entropy_samples_order_Q_north,
                                  mid = entropy_samples_order_Q_mid),
                             function(samples_list, metadat) {
                               metadata_subset <- metadat[metadat$sample_code %in% samples_list$V1, c('sample_code', 'location', 'sciname2', 'long', 'lat', 'species2')]
                               metadata_subset_reorder <- metadata_subset[match(as.character(samples_list$V1), metadata_subset$sample_code),]
                               
                               metadata_subset_reorder %>% 
                                 mutate(location_clean = factor(case_when(location == 'Gombe_South'~ "Gombe S",
                                                                          location == 'Katongwe_N' ~ "Katongwe N",
                                                                          location == 'Katongwe_S' ~ "Katongwe S",
                                                                          location == 'Ska' ~ "S Kagango",
                                                                          location == 'Kalala' ~ "Kalalangabo",
                                                                          location == 'Nondwa' ~ "Nondwa",
                                                                          location == 'Hilltop' ~ "Hilltop",
                                                                          location == 'Bangwe' ~ "Bangwe",
                                                                          location == 'Jakob' ~ "Jakobsen's S",
                                                                          location == 'Ulombola' ~"Ulombola"),
                                                                levels= c("Gombe S", "Katongwe N", "Katongwe S", "S Kagango", "Kalalangabo", "Nondwa", "Hilltop", "Bangwe", "Jakobsen's S", "Ulombola", "Harembe")),
                                        location_plotting = factor(case_when(location == 'Gombe_South'~ "A",
                                                                             location == 'Katongwe_N' ~ "B",
                                                                             location == 'Katongwe_S' ~ "C",
                                                                             location == 'Ska' ~ "D",
                                                                             location == 'Kalala' ~ "E",
                                                                             location == 'Nondwa' ~ "F",
                                                                             location == 'Hilltop' ~ "G",
                                                                             location == 'Bangwe' ~ "H",
                                                                             location == 'Jakob' ~ "I",
                                                                             location == 'Ulombola' ~"J"),
                                                                   levels= LETTERS[1:10]),
                                        reduced_species_names = ifelse(sciname2 == "Petrochromis kazumbe",
                                                                       "P. kazumbe",
                                                                       "P. polyodon")
                                 )
                               
                             }, metadat = cichlid_metadata)


### Processing entropy results ###

#entropy_upload_list k2_Q_q_upload_list
#full dataset q results
entropy_list <- lapply(split(names(entropy_upload_list), names(entropy_upload_list)), 
                       function(k, metadata_subset_reorder, entropy_list) {
                         combined.q <- abind(entropy_upload_list[[k]][['chain_1']], entropy_upload_list[[k]][['chain_2']], entropy_upload_list[[k]][['chain_3']], along=1)
                         q_estimate <- as.data.frame(t(apply(combined.q, 2:3, mean)))
                         q_estimate_with_metadata <- cbind(q_estimate, metadata_subset_reorder)
  
                         q_long <- q_estimate_with_metadata %>%
                           gather(key = cluster, value = proportion, paste0('V', 1:as.numeric(gsub('K_', '', k))) )
                        
                         list(combined.q = combined.q,
                              q_estimate = q_estimate,
                              q_estimate_with_metadata = q_estimate_with_metadata,
                              q_long = q_long)
                         }, 
                       metadata_subset_reorder = processed_metadata$full,
                       entropy_list = entropy_upload_list)


## Get credible intervals and calculate mean

# entropy_list$K_2$combined.q
# 
# 
# k2_q_Q_list$north$Q$combined.Q
# 
# q.ci <- apply(entropy_list$K_3$combined.q, 2:3, quantile, probs = c(0.025,0.975))
# q.ci.width <- q.ci[2,2,] - q.ci[1,2,]
# 
# 
# Q.ci <- apply(allQ, 2:3, quantile, probs=c(0.025,0.975))
# q.ci.width <- q.ci[2,2,] - q.ci[1,2,]
# Q.ci.width <- Q.ci[2,2,] - Q.ci[1,2,]
# 
# mean(q.ci.width)


#focal_var is the ancestry parameter that I will focus on --> it is the proportion of polyodon ancestry
focal_var <- "V2"

k2_q_Q_list <- lapply(list(mid = 'mid', north = 'north'), function(region, entropy_results, meta_data, focal_var) {
  
  col <- switch(focal_var, V1 = {1}, V2 = {2})
  
  k2_Qmod_combined.q <- abind(entropy_results[[paste0(region, '_q')]]$chain_1, entropy_results[[paste0(region, '_q')]]$chain_2, entropy_results[[paste0(region, '_q')]]$chain_3, along = 1)
  k2_Qmod_q_estimate <- as.data.frame(t(apply(k2_Qmod_combined.q, 2:3, mean)))
  k2_Qmod_q_estimate_with_metadata <- cbind(k2_Qmod_q_estimate, meta_data[[region]])
  
  k2_combined.Q <- abind(entropy_results[[paste0(region, '_Q')]]$chain_1, entropy_results[[paste0(region, '_Q')]]$chain_2, entropy_results[[paste0(region, '_Q')]]$chain_3, along = 1)
  
  q_ci <- apply(k2_Qmod_combined.q, 2:3, quantile, probs = c(0.025,0.975)) #95 CI% for q
  Q_ci <- apply(k2_combined.Q, 2:3, quantile, probs = c(0.025,0.975)) #95 CI% for Q
  
  q_collapse <- do.call(rbind, lapply(seq(dim(q_ci)[3]), function(x) q_ci[ , col, x]))
  
  ci_df <- do.call(rbind, lapply(seq(dim(Q_ci)[3]), function(x) Q_ci[ , , x])) %>% 
    as.data.frame() %>% 
    rownames_to_column() %>%
    mutate(Q_ci = paste0('Qci_', sub(pattern = 'X([0-9]+\\.[0-9]).*', replacement = '\\1', x = rowname)),
           index = rep(1:(n()/2), each = 2)) %>%
    select(!c(rowname, V1, V3) ) %>% 
    pivot_wider(names_from = 'Q_ci', values_from = 'V2') %>%
    cbind(., q_collapse %>%
            as.data.frame() %>% 
            rename(qci_2.5 = `2.5%`,
                   qci_97.5 = `97.5%`))
  
  k2_Qmod_q_estimate_with_metadata_andci <- cbind(k2_Qmod_q_estimate_with_metadata, ci_df)
  
  return(
    list(q = list(combined.q = k2_Qmod_combined.q,
                  q_estimate = k2_Qmod_q_estimate,
                  q_estimate_with_metadata = k2_Qmod_q_estimate_with_metadata,
                  q_estimate_with_metadata_and_ci = k2_Qmod_q_estimate_with_metadata_andci,
                  q_long = k2_Qmod_q_estimate_with_metadata %>%
                    gather(key = cluster, value = proportion, c('V1', 'V2')) ),
         Q = list(combined.Q = k2_combined.Q,
                  Q_estimate = data.frame(t(apply(k2_combined.Q, 2:3, mean)))) )
  )
}, entropy_results = k2_Q_q_upload_list, meta_data = processed_metadata, focal_var = focal_var)




# q.ci.width <- q.ci[2,2,] - q.ci[1,2,]


#checking for evidence of label switching
#V1 should be the dominant ancestry for the same species in region; if this isn't the case, label switching occurred
k2_q_Q_list$north$q$q_estimate_with_metadata %>% 
  group_by(reduced_species_names) %>% 
  summarise(mean_V1 = mean(V1), 
            mean_V2 = mean(V2))

cbind(k2_q_Q_list$north$q$q_estimate_with_metadata, k2_q_Q_list$north$Q$Q_estimate) %>% 
  group_by(reduced_species_names) %>% 
  summarise(mean_X1 = mean(X1), 
            mean_X2 = mean(X2),
            mean_X3 = mean(X3))
  
k2_q_Q_list$mid$q$q_estimate_with_metadata %>% 
  group_by(reduced_species_names) %>% 
  summarise(mean_V1 = mean(V1), 
            mean_V2 = mean(V2))

cbind(k2_q_Q_list$mid$q$q_estimate_with_metadata, k2_q_Q_list$mid$Q$Q_estimate) %>% 
  group_by(reduced_species_names) %>% 
  summarise(mean_X1 = mean(X1), 
            mean_X2 = mean(X2),
            mean_X3 = mean(X3))

#is there evidence of label switching?
label_switch <- 'no'

if (label_switch == 'yes') {
  #switch the labels for V1 <-> V2 and X1 <-> X2 for the mid region to harmonize the labels between the regions
  k2_q_Q_combined <- rbind(cbind(k2_q_Q_list$north$q$q_estimate_with_metadata_and_ci, k2_q_Q_list$north$Q$Q_estimate),
                           cbind(k2_q_Q_list$mid$q$q_estimate_with_metadata_and_ci, k2_q_Q_list$mid$Q$Q_estimate) %>%
                                mutate(V1_new = V2,
                                       V2_new = V1,
                                       X1_new = X3,
                                       X3_new = X1) %>%
                                select(!c(X1, X3, V1, V2)) %>%
                                rename(V1 = V1_new,
                                       V2 = V2_new,
                                       X1 = X1_new,
                                       X3 = X3_new) %>%
                                select(c(V1, V2, sample_code, location, sciname2, long, lat, species2, location_clean, location_plotting, reduced_species_names, X1, X2, X3))
  )

} else {
  ## IF THERE IS NO LABEL SWITCHING, YOU CAN JUST BIND THE DATA FOR THE REGIONS TOGETHER
  k2_q_Q_combined <- rbind(cbind(k2_q_Q_list$north$q$q_estimate_with_metadata_and_ci, k2_q_Q_list$north$Q$Q_estimate),
                           cbind(k2_q_Q_list$mid$q$q_estimate_with_metadata_and_ci, k2_q_Q_list$mid$Q$Q_estimate))

}



#focal_var is the ancestry parameter that I will focus on --> it is the proportion of polyodon ancestry
#focal_var <- "V2"
#X2 --> prop of interspecific hybridization
#V2 --> prop of polyodon genome


### Classifying into hybrid classes based on Q12 and q ###
#From Mandeville et al. (2019):
#equation quantifying the proximity of an individual to the expected line between parental species: Q - 2q; Q + 2q
#kazumbe backcross: −0.1 < Q − 2q < 0.1
#polyodon backcross: 1.9 < Q + 2q < 2.1. 

#full polyodon: q > 0.9 and Q ≤ 0.25. 
#full kazumbe: q < 0.1 and Q ≤ 0.25. 
#F1 hybrids: q 0.4–0.6 and Q > 0.8. 
#F2 hybrids: q 0.4–0.6 (the same as F1 hybrids) and Q 0.4–0.6

#Other: did not meet any of these conditions for specific hybrid classes, but were intermediate between parental species.
k2_q_Q_combined$Q_minus_2q <- k2_q_Q_combined$X2 - (2*k2_q_Q_combined[, focal_var])
k2_q_Q_combined$Q_plus_2q <- k2_q_Q_combined$X2 + (2*k2_q_Q_combined[, focal_var])

k2_q_Q_combined$classification <- NA

k2_q_Q_combined$classification[ (k2_q_Q_combined[, focal_var] > 0.9) & (k2_q_Q_combined$X2 <= 0.25)] <- "polyodon"
k2_q_Q_combined$classification[ (k2_q_Q_combined[, focal_var] < 0.1) & (k2_q_Q_combined$X2 <= 0.25)] <- "kazumbe"
k2_q_Q_combined$classification[ (k2_q_Q_combined[, focal_var] > 0.4) & (k2_q_Q_combined[, focal_var] < 0.6) & (k2_q_Q_combined$X2 > 0.8)] <- "F1"
k2_q_Q_combined$classification[ (k2_q_Q_combined[, focal_var] > 0.4) & (k2_q_Q_combined[, focal_var] < 0.6) & (k2_q_Q_combined$X2 > 0.4) & (k2_q_Q_combined$X2 < 0.6)] <- "F2"

k2_q_Q_combined$classification[is.na(k2_q_Q_combined$classification) & (k2_q_Q_combined$Q_minus_2q > -0.1) & (k2_q_Q_combined$Q_minus_2q < 0.1)] <- "bc.kazumbe"
k2_q_Q_combined$classification[is.na(k2_q_Q_combined$classification) & (k2_q_Q_combined$Q_plus_2q > 1.9) & (k2_q_Q_combined$Q_plus_2q < 2.1)] <- "bc.polyodon"
k2_q_Q_combined$classification[ is.na(k2_q_Q_combined$classification) ] <- "other"
k2_q_Q_combined$classification <- factor(k2_q_Q_combined$classification, levels = c("kazumbe", "polyodon", "F1", "F2", "bc.kazumbe", "bc.polyodon", "other"))


#number of individuals per sampling location
location_sample_size <- as.data.frame(table(droplevels(k2_q_Q_combined$location_plotting)))



##############################
### Creating Q vs. q plots ###
##############################

### Aesthetic decisions ###
title_size <- 25 #title size of multipanel plot (for entropy_bar_plot_multi)
vertical_adjust_title <- -1 #vertical adjustment of title in multipanel plot (for entropy_bar_plot_multi)

### Information and set-up ###
#X2 --> prop of interspecific hybridization
#V1 --> prop of kazumbe genome
#V2 --> prop of polyodon genome

#triangle for Q vs. q plots
triangle_dataframe <- data.frame(x = c(0, 0.5, 1), y = c(0, 1, 0))

#creating common y and x axis titles for multipanel Q vs. q plots
#https://stackoverflow.com/questions/33114380/centered-x-axis-label-for-muliplot-using-cowplot-package
#https://stackoverflow.com/questions/44152058/wider-margins-for-grid-arrange-function
y.grob <- textGrob(expression(paste(Q[12], ' (interspecific ancestry)')), 
                   gp=gpar(fontface="plain", col="black", fontsize=20), rot=90)

x.grob <- textGrob("q (proportion of ancestry)", 
                   gp=gpar(fontface="plain", col="black", fontsize=20))


### Creating entropy plots ###
k2_plot_list <- list() #list to store K2 Q vs. q plots

k2_plot_list[["k2_full_dataset"]] <- ggplot(data = k2_q_Q_combined, aes_string(x = focal_var, y = "X2", color = "location_plotting") ) +
  geom_line(data = triangle_dataframe, aes(x = x, y = y), color = '#e5e5e5', size = 1.5, lineend = "round") +
  geom_point(size = 6.5, alpha = 0.4) +
  geom_point(size = 6.5, shape = 1, alpha = 0.4) +
  scale_color_manual(name = "Location", values = location_palette_letters[names(location_palette_letters) %in% unique(k2_q_Q_combined$location_plotting) ]) +
  scale_x_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
  scale_y_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
  xlab("q (proportion of ancestry)") +
  ylab(expression(Q[12]~"(interspecific ancestry)") ) +
  theme_Q_entropy_alt() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        axis.title = element_text(size = 20, face = "plain"),
        axis.text = element_text(size = 16),
        plot.margin=unit(c(5, 1.5, 1, 2),"mm"),
        legend.position = 'none')


#loop through each location with at least two individuals of each species, create a Q vs. q plot, and store the plot in k2_plot_list
for (i in unique(k2_q_Q_combined$location_plotting)) {
  
  #subset down to location i
  location_dataframe <- k2_q_Q_combined %>%
    filter(location_plotting == i)
  
  #create plot and store in k2_plot_list in the entry "k2_i" where i is the location (e.g. k2_Hilltop)
  k2_plot_list[[paste0("k2_", gsub(" ", "_", i) )]] <- ggplot(data = location_dataframe, aes_string(x = focal_var, y = "X2", color = "location_plotting") ) +
    geom_line(data = triangle_dataframe, aes(x = x, y = y), color = '#e5e5e5', size = 1.5, lineend = "round") +
    geom_point(size = 6.5, alpha = 0.4) +
    geom_point(size = 6.5, shape = 1, alpha = 0.4) +
    scale_color_manual(name = "Location", values = location_palette_letters[names(location_palette_letters) %in% unique(location_dataframe$location_plotting) ]) +
    scale_x_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
    scale_y_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
    labs(title = i) +
    theme_Q_entropy() +
    theme(axis.title = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 24),
          legend.position = "none") +
    geom_segment(aes(x = qci_2.5, xend = qci_97.5, y = X2, yend = X2), size = 1) +
    geom_segment(aes(y = Qci_2.5, yend = Qci_97.5, x = V2, xend = V2), size = 1)
  
  #print message saying that the plot for location i is finished
  message(paste0("completed plot for ", i))
}


### barchart of hybrid class frequencies broken down by sampling location ###
#to add custom text to (i.e. the sample sizes of each location placed above each bar), you need to put the data and aes in the geom_bar() function, not the ggplot() function
hybrid_distribution_barchart <- ggplot() + 
  geom_bar(data = k2_q_Q_combined, aes(x = location_plotting, fill = classification), position = "fill") +
  scale_fill_manual(values = c("#DC7633", "#3498DB", "#edba99", "#99cbed", "#66b266")) +
  theme(#legend.position = "none", 
    plot.background = element_rect(fill = "white"), 
    panel.background = element_rect(fill = "white"),
    axis.text.y=element_text(size = 18, color = 'black'),
    #axis.text.x=element_text(angle = 45, margin=margin(25,0, -40,0), hjust=1, size = 14.4, color = 'black', vjust = 1.6),
    axis.text.x=element_text(hjust=0.5, size = 18, color = 'black', vjust = 3.5),
    axis.ticks.x=element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 20, face = "plain"),
    legend.text = element_text(size = 15),
    legend.title = element_blank(),
    plot.margin=unit(c(3,0,0,0),"mm")) +
  ylab("Proportion of individuals") +
  geom_text(aes(x = Var1, y = 1, label = Freq), size = 6, data = location_sample_size, vjust = -0.4) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
  coord_cartesian(clip = "off")


### Creating full multipanel plot ###

#ROW 1 --> (left) FULL DATA Q12 vs. q plot and (right) hybrid frequency barchart
location_palette_letters_processed <- location_palette_letters[names(location_palette_letters) %in% unique(k2_q_Q_combined$location_plotting) ]
location_palette_letters_processed <- location_palette_letters_processed[match(sort(names(location_palette_letters_processed)), names(location_palette_letters_processed))] 
legend_Q_q_full_plot_alt <- get_legend(ggplot(data = k2_q_Q_combined, aes_string(x = focal_var, y = "X2", color = "location_plotting") ) +
                                         geom_line(data = triangle_dataframe, aes(x = x, y = y), color = '#bbbaba', size = 1.5, lineend = "round") +
                                         geom_point(size = 6.5, alpha = 0.5) +
                                         geom_point(size = 6.5, shape = 1, alpha = 0.5) +
                                         scale_color_manual(name = "Location", values = location_palette_letters_processed) +
                                         scale_x_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
                                         scale_y_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
                                         xlab("q (proportion of ancestry)") +
                                         ylab(expression("Q[12] (interspecific ancestry)") ) +
                                         theme_Q_entropy_alt() +
                                         theme(plot.title = element_text(hjust = 0.5, size = 16),
                                               axis.title = element_text(size = 16, face = "plain"),
                                               legend.text = element_text(size = 18),
                                               legend.title = element_text(size = 20),
                                               plot.margin=unit(c(5, 0, 1, 4), "mm")))

Q_q_full_plot_alt_marg_hist <- ggMarginal(k2_plot_list[["k2_full_dataset"]], type = "histogram",
                                          xparams = list(binwidth = 0.05, fill = '#bbbaba', col = 'black'),
                                          yparams = list(binwidth = 0.05, fill = '#bbbaba', col = 'black'))

Q_q_full_plot_alt_marg_hist_withlegend <- cowplot::plot_grid(Q_q_full_plot_alt_marg_hist, ggdraw(legend_Q_q_full_plot_alt) + theme(plot.margin=unit(c(0, 0, 0, 1),"cm")), 
                                                             rel_widths = c(2.6, 0.5))

upper_plot_Qq_bar_withmarg <- cowplot::plot_grid(Q_q_full_plot_alt_marg_hist_withlegend, hybrid_distribution_barchart, 
                                                 nrow = 1, labels = c("(a)", "(b)"), label_size = 25, scale = 0.9, hjust = c(-0.5, 0))



#ROWS 2 AND 3 --> LOCATION-SPECIFIC Q12 vs q PLOTS 
#arrange each location's Q vs. q plot into a single multi-panel plot

k2_Q_q_models <- plot_grid(k2_plot_list$k2_A, k2_plot_list$k2_B,
                           k2_plot_list$k2_C, k2_plot_list$k2_D,
                           k2_plot_list$k2_E, k2_plot_list$k2_F,
                           k2_plot_list$k2_G, k2_plot_list$k2_I, nrow = 2, align = 'vh', vjust=1, scale = 1) +
  theme(plot.margin = margin(0, 0, 0, 5))

#add the common x and y axis titles for the multipanel plot of location-specific Q vs. q plots
k2_Q_q_multipanel_lower <- grid.arrange(arrangeGrob(k2_Q_q_models, left = y.grob, bottom = x.grob),
                                        vp = viewport(width = 0.97, height = 0.97))

Qq_multipanel_FULL_withbarchart_withmarg <- cowplot::plot_grid(upper_plot_Qq_bar_withmarg,
                                                               ggdraw(k2_Q_q_multipanel_lower), nrow = 2, rel_heights = c(2, 2.6),
                                                               labels = c("", "(c)"),
                                                               label_size = 25)


#ROWS 4 --> BARPLOTS OF q for K = 2 and K = 3
qinfo_extract <- lapply(entropy_list, function(x) x$q_long)

#ADDRESSING KATIE'S COMMENT (11/14/2021): switch ancestry group labels for K = 3 so that the polyodon
#individuals have the same color (blue) in both the K = 2 and K = 3 visualizations
qinfo_extract$K_3 <- qinfo_extract$K_3 %>% 
  mutate(cluster = case_when(cluster == 'V1' ~ 'V2',
                             cluster == 'V2' ~ 'V1',
                             cluster == 'V3' ~ 'V3') )

qinfo_extract$K_4 <- qinfo_extract$K_4 %>% 
  mutate(cluster = case_when(cluster == 'V1' ~ 'V3',
                             cluster == 'V3' ~ 'V1',
                             cluster == 'V2' ~ 'V2',
                             cluster == 'V4' ~ 'V4') )

qinfo_extract$K_5 <- qinfo_extract$K_5 %>% 
  mutate(cluster = case_when(cluster == 'V1' ~ 'V3',
                             cluster == 'V3' ~ 'V1',
                             cluster == 'V2' ~ 'V2',
                             cluster == 'V4' ~ 'V4',
                             cluster == 'V5' ~ 'V5') )

V2_plots <- create_V2_plots(qinfo_extract, q_col_pal = q_color_palette)

entropy_bar_plot_multi <- plot_grid(V2_plots$plot_list$K_2location_genetic + ggtitle("") +
                                      #ggtitle(paste0("K = 2; DIC = ", formatC(DIC_dataframe$DIC[2], format="e", digits = 2)) ) +
                                      ggtitle(expression(italic('K') ~ " = 2")) +
                                      theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size),
                                            plot.margin=unit(c(0.1, 0.35, 0.5, 0),"cm")) , 
                                    V2_plots$plot_list$K_3location_genetic + ggtitle("") +
                                      #ggtitle(paste0("K = 3; DIC = ", formatC(DIC_dataframe$DIC[3], format="e", digits = 2)) ) +
                                      ggtitle(expression(italic('K') ~ " = 3")) +
                                      theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size),
                                            plot.margin=unit(c(0.1, 0, 0.5, 0.35),"cm")), ncol = 2)



### CREATING AND EXPORTING FINAL ENTROPY PLOT ###
plot_grid(Qq_multipanel_FULL_withbarchart_withmarg, entropy_bar_plot_multi + theme(plot.margin = unit(c(0.85, 0.7, 0.3, 0.8),"cm")), 
          labels = c("", "(d)"), label_y = c(1, 1.01), nrow = 2, label_size = 25, rel_heights = c(0.8, 0.25))

#ggsave2(here('figures', 'k2_Q_q_multipanel_FULL_barchart_marginal_qbarplot_7_18_2022.png'), 
#        width = 19, height = 14, bg = "white")

ggsave2(here('figures', 'k2_Q_q_multipanel_FULL_barchart_marginal_qbarplot_7_18_2022.pdf'), 
        device = 'pdf', width = 19, height = 14, bg = "white", dpi = 1000)

# ggsave2(here('figures', 'k2_Q_q_multipanel_FULL_barchart_marginal_qbarplot_9_7_2021.png'), 
#         width = 19, height = 14, bg = "white")



##############################################################################
### PLOTS OF q VALUES FOR K = 2 to K = 5 (INCLUDED SUPPLEMENTARY MATERIAL) ###
##############################################################################

cowplot::plot_grid(V2_plots$plot_list$K_2location_genetic +
                     theme(plot.margin = unit(c(0.1, 0.25, 0.45, 0.25), "cm"),
                           plot.title = element_text(size = 24.5, vjust = -1)) +
                     ggtitle(expression(italic('K') ~ " = 2")),
                   V2_plots$plot_list$K_3location_genetic+
                     theme(plot.margin = unit(c(0.1, 0.25, 0.45, 0.25), "cm"),
                           plot.title = element_text(size = 24.5, vjust = -1)) +
                     ggtitle(expression(italic('K') ~ " = 3")),
                   V2_plots$plot_list$K_4location_genetic+
                     theme(plot.margin = unit(c(0.1, 0.25, 0.45, 0.25), "cm"),
                           plot.title = element_text(size = 24.5, vjust = -1)) +
                     ggtitle(expression(italic('K') ~ " = 4")),
                   V2_plots$plot_list$K_5location_genetic+
                     theme(plot.margin = unit(c(0.1, 0.25, 0.45, 0.25), "cm"),
                           plot.title = element_text(size = 24.5, vjust = -1)) +
                     ggtitle(expression(italic('K') ~ " = 5")), nrow = 4) +
  theme(plot.margin = unit(c(0, 0, 0.25, 0), "cm"))

#ggsave(here('figures', 'entropy_k2through5.png'),
#       width = 12, height = 10, bg = "white")

ggsave(here('figures', 'entropy_k2through5.pdf'), device = 'pdf',
       width = 12, height = 10, bg = "white", dpi = 1000)



#######################################
### TRACEPLOTS OF CHAINS FOR EACH K ###
#######################################

#For the K = 2 and K = 3 models, which are visualized in the main text, I plotted traces of 4
#randomly chosen parameters

#I set a random seed so that the plot included in the supplementary material is replicable. If
#you want to randomly select a different set of parameters, either don't set a random seed or
#change the random seed number.

set.seed(59021)

trace_multipanel_list <- list()

#loop through K = 2 and K = 3
for (i in 1:2) {
  
  trace_plot_list <- list()
  k <- i + 1
  
  #loop through 4 random parameters of the model
  for (x in 1:4) {
    parameter_choice <- sample(seq(k), size = 1) #choose a random parameter
    individual_choice <- sample(ncol(entropy_upload_list[[i]][[1]][1,,]), size = 1) #choose a random individual
    
    trace_title <- ggdraw() + 
      draw_label(
        paste0('K = ', k),
        fontface = 'bold',
        x = 0,
        hjust = 0, 
        size = 30
      ) +
      theme(
        # add margin on the left of the drawing canvas,
        # so title is aligned with left edge of first plot
        plot.margin = margin(0, 0, 0, 7)
      )
    
    trace_plot_list[[x]] <- ggplot() +
      geom_line(data = data.frame(Iteration = 1:length(entropy_upload_list[[i]][[1]][, parameter_choice, individual_choice]),
                                  Parameter = entropy_upload_list[[i]][[1]][, parameter_choice, individual_choice]),
                aes(x = Iteration, y = Parameter), color = "#ef476f", size = 1, alpha = 1) +
      geom_line(data = data.frame(Iteration = 1:length(entropy_upload_list[[i]][[2]][, parameter_choice, individual_choice]),
                                  Parameter = entropy_upload_list[[i]][[2]][, parameter_choice, individual_choice]),
                aes(x = Iteration, y = Parameter), color = "#06d6a0", size = 1, alpha = 0.85) +
      geom_line(data = data.frame(Iteration = 1:length(entropy_upload_list[[i]][[3]][, parameter_choice, individual_choice]),
                                  Parameter = entropy_upload_list[[i]][[3]][, parameter_choice, individual_choice]),
                aes(x = Iteration, y = Parameter), color = "#118ab2", size = 1, alpha = 0.7) +
      #ggtitle(paste0('k = ', k, '; param = ', parameter_choice, '; indiv = ', individual_choice)) +
      ylim(0, 1) +
      theme_cowplot() +
      theme(axis.title = element_text(size = 14),
            axis.text = element_text(size = 14))
  }
  
  trace_multipanel <- plot_grid(trace_plot_list[[1]], trace_plot_list[[2]],
                                trace_plot_list[[3]], trace_plot_list[[4]], nrow = 2)
  
  trace_multipanel_list[[paste0('k_', k)]] <- plot_grid(trace_title, trace_multipanel,
                                                        ncol = 1,
                                                        rel_heights = c(0.15, 1))
  
  message('completed the following k: ', k)
}

#create final multipanel plot of all traces
example_trace_plot <- plot_grid(
  trace_multipanel_list$k_2,
  trace_multipanel_list$k_3,
  nrow = 2
)

example_trace_plot
#ggsave2(here('figures', 'entropy_trace_plots.png'), 
#        width = 0.8*19, height = 0.8*14, bg = "white")
ggsave2(here('figures', 'entropy_trace_plots.pdf'), device = 'pdf', 
        width = 0.8*19, height = 0.8*14, bg = "white")




#################################
### RESULTS SECTION SUMMARIES ###
#################################

#X2 --> prop of interspecific hybridization
#V1 --> prop of kazumbe genome
#V2 --> prop of polyodon genome

summary_new_classify_df <- k2_q_Q_combined %>% 
  mutate(updated_id = case_when(V1 > 0.5 ~ 'kazumbe',
                                V1 < 0.5 ~ 'polyodon'))

#percentage of individuals with majority kazumbe ancestry or polyodon ancestry
#that are not categorized into a hybrid class
summary_new_classify_df %>%
  filter(!(classification %in% c('bc.polyodon', 'bc.kazumbe', 'other')) ) %>% 
  group_by(updated_id) %>%
  summarize(n=n()) %>% 
  ungroup() %>% 
  left_join(., summary_new_classify_df %>% group_by(updated_id) %>% summarize(n=n()),
            by = 'updated_id') %>% 
  rename(species_id = updated_id,
         hybrid = n.x,
         total = n.y) %>% 
  mutate(nonhybrid_prop_percentage = (hybrid/total)*100)

#percentage of individuals in dataset that are classified as having complete
#ancestry of one of the parental species (i.e either full polyodon ancestry or
#full kazumbe ancestry)
summary_new_classify_df %>%
  filter(!(classification %in% c('bc.polyodon', 'bc.kazumbe', 'other'))) %>%
  summarize(n=n()) %>%
  pull(n)/(summary_new_classify_df %>% summarize(n=n()) %>% pull(n))*100

#counts and proportion of each individuals that have hybrid ancestry, broken
#down by majority parental species ancestry
summary_new_classify_df %>%
  filter(classification %in% c('bc.polyodon', 'bc.kazumbe', 'other') ) %>% 
  group_by(updated_id) %>%
  summarize(n=n()) %>% 
  ungroup() %>% 
  left_join(., summary_new_classify_df %>% group_by(updated_id) %>% summarize(n=n()),
            by = 'updated_id') %>% 
  rename(species_id = updated_id,
         hybrid = n.x,
         total = n.y) %>% 
  mutate(hybrid_prop_percentage = (hybrid/total)*100)

#counts and proportion of each individuals that have hybrid ancestry, broken
#down by majority parental species ancestry and broken down by sampling
#location
summary_new_classify_df %>%
  filter(classification %in% c('bc.polyodon', 'bc.kazumbe', 'other')) %>% 
  group_by(location_plotting) %>%
  summarize(n=n()) %>% 
  ungroup() %>% 
  full_join(., summary_new_classify_df %>% group_by(location_plotting) %>% summarize(n=n()),
            by = 'location_plotting') %>% 
  rename(location = location_plotting,
         hybrid = n.x,
         total = n.y) %>%
  mutate(hybrid = replace_na(hybrid, 0)) %>% 
  mutate(hybrid_percentage = (hybrid/total)*100)



#################################
### CODE CURRENTLY NOT IN USE ###
#################################

### TABLE OF DIC VALUES ###
# DIC_dataframe <- read.table(here('entropy', 'DIC_values.txt'), header = TRUE)
# DIC_dataframe_processed <- DIC_dataframe %>% 
#   mutate(DIC = format(DIC, scientific = TRUE, big.mark = ",", digits = 3),
#          delta_DIC = format(delta_DIC, scientific = TRUE, big.mark = ",", digits = 3)) %>% 
#   rename("$\\Delta$DIC" = delta_DIC)
# 
# entropy_dic_table <- kable(DIC_dataframe_processed, "latex", booktabs = TRUE, align = "c", escape = FALSE) %>% 
#   kable_styling(latex_options = c("hold_position"))
# cat(entropy_dic_table, file = paste0(table_save_dir, 'entropy_dic_table.txt'), append = FALSE)
# 
# 


## Get credible intervals and calculate mean

# entropy_list$K_2$combined.q
# 
# 
# k2_q_Q_list$north$Q$combined.Q
# 
# q.ci <- apply(entropy_list$K_3$combined.q, 2:3, quantile, probs = c(0.025,0.975))
# q.ci.width <- q.ci[2,2,] - q.ci[1,2,]
# 
# 
# Q.ci <- apply(allQ, 2:3, quantile, probs=c(0.025,0.975))
# q.ci.width <- q.ci[2,2,] - q.ci[1,2,]
# Q.ci.width <- Q.ci[2,2,] - Q.ci[1,2,]
# 
# mean(q.ci.width)


# k2_Q_q_multipanel_FULL_withmarg <- cowplot::plot_grid(cowplot::plot_grid(NULL, Q_q_full_plot_alt_marg_hist_withlegend, NULL, rel_widths = c(1.3, 3, 0.7), nrow = 1, labels = c("", "(a)", ""),
#                                                                          label_size = 25, hjust = 1.2),
#                                                       ggdraw(k2_Q_q_multipanel_lower), nrow = 2, rel_heights = c(2.5, 3.5),
#                                                       labels = c("", "(b)"),
#                                                       label_size = 25)


#X2 --> prop of interspecific hybridization
#V1 --> prop of kazumbe genome
#V2 --> prop of polyodon genome



# ### K2 ###
# #Q vs. q plot of the full dataset (all individuals of polyodon and kazumbe)
# Q_q_full_plot <- ggplot(data = k2_q_Q_combined, aes(x = V1, y = X2, color = location_plotting) ) +
#   geom_line(data = triangle_dataframe, aes(x = x, y = y), color = '#bbbaba', size = 1.5, lineend = "round") +
#   geom_point(size = 6.5, alpha = 0.5) +
#   geom_point(size = 6.5, shape = 1, alpha = 0.5) +
#   scale_color_manual(name = "Location", values = location_palette_letters[names(location_palette_letters) %in% unique(k2_q_Q_combined$location_plotting) ]) +
#   scale_x_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
#   scale_y_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
#   xlab("q (proportion of ancestry)") +
#   ylab("Q (interspecific ancestry)") +
#   theme_Q_entropy() +
#   theme(plot.title = element_text(hjust = 0.5, size = 16),
#         axis.title = element_text(size = 20, face = "plain"),
#         axis.text = element_text(size = 18, color = 'black'))
# 
# 
# k2_plot_list <- list() #list to store K2 Q vs. q plots
# k2_plot_list[["k2_full_dataset"]] <- Q_q_full_plot #first entry is the full dataset Q vs. q plot
# 
# #loop through each location with at least two individuals of each species, create a Q vs. q plot, and store the plot in k2_plot_list
# for (i in unique(k2_q_Q_combined$location_plotting)) {
#   
#   #subset down to location i
#   location_dataframe <- k2_q_Q_combined %>%
#     filter(location_plotting == i)
#   
#   #create plot and store in k2_plot_list in the entry "k2_i" where i is the location (e.g. k2_Hilltop)
#   k2_plot_list[[paste0("k2_", gsub(" ", "_", i) )]] <- ggplot(data = location_dataframe, aes(x = V1, y = X2, color = location_plotting) ) +
#     geom_line(data = triangle_dataframe, aes(x = x, y = y), color = '#bbbaba', size = 1.5, lineend = "round") +
#     geom_point(size = 6.5, alpha = 0.5) +
#     geom_point(size = 6.5, shape = 1, alpha = 0.5) +
#     scale_color_manual(name = "Location", values = location_palette_letters[names(location_palette_letters) %in% unique(location_dataframe$location_plotting) ]) +
#     scale_x_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
#     scale_y_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
#     labs(title = i) +
#     theme_Q_entropy() +
#     theme(axis.title = element_blank(),
#           plot.title = element_text(hjust = 0.5, size = 24),
#           legend.position = "none")
#   
#   #print message saying that the plot for location i is finished
#   message(paste0("completed plot for ", i))
#   
# }
# 
# 
# 
# #arrange each location's Q vs. q plot into a single multi-panel plot
# k2_Q_q_models <-plot_grid(k2_plot_list$k2_A, k2_plot_list$k2_B,
#                           k2_plot_list$k2_C, k2_plot_list$k2_D,
#                           k2_plot_list$k2_E, k2_plot_list$k2_F,
#                           k2_plot_list$k2_G, k2_plot_list$k2_I, nrow = 2, align='vh', vjust=1, scale = 1) +
#   theme(plot.margin = margin(0, 0, 0, 5))
# 
# #add the common x and y axis titles for the multipanel plot of location-specific Q vs. q plots
# k2_Q_q_multipanel_lower <- grid.arrange(arrangeGrob(k2_Q_q_models, left = y.grob, bottom = x.grob), vp=viewport(width=0.97, height=0.97))
# 
# #combine the full dataset Q vs q plot as (a) and the location-specific multi-plot panel of Q vs. q plots as (b)
# k2_Q_q_multipanel_FULL <- cowplot::plot_grid(cowplot::plot_grid(NULL, k2_plot_list$k2_full_dataset, NULL, rel_widths = c(1.6, 3, 1.6), nrow = 1, labels = c("", "(a)", ""), 
#                                                                 label_size = 25, hjust = 1.2), 
#                                              ggdraw(k2_Q_q_multipanel_lower), nrow = 2, rel_heights = c(2,3), 
#                                              labels = c("", "(b)"), 
#                                              label_size = 25)
# 
# #Export the completed K2 Q vs. q plot
# k2_Q_q_multipanel_FULL
# ggsave2("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/entropy/results/results_5_19_2021/plots/k2_Q_q_multipanel_FULL_5_19_2021.png", width = 20, height = 12)
# ##############################

#X2 --> prop of interspecific hybridization
#V1 --> prop of polyodon genome


#From Mandeville et al. (2019):
#equation quantifying the proximity of an individual to the expected line between parental species: Q - 2q; Q + 2q
#kazumbe backcross: −0.1 < Q − 2q < 0.1
#polyodon backcross: 1.9 < Q + 2q < 2.1. 

#full polyodon: q > 0.9 and Q ≤ 0.25. 
#full kazumbe: q < 0.1 and Q ≤ 0.25. 
#F1 hybrids: q 0.4–0.6 and Q > 0.8. 
#F2 hybrids: q 0.4–0.6 (the same as F1 hybrids) and Q 0.4–0.6

#Other: did not meet any of these conditions for specific hybrid classes, but were intermediate between parental species.
# k2_q_Q_combined$Q_minus_2q <- k2_q_Q_combined$X2 - (2*k2_q_Q_combined$V1)
# k2_q_Q_combined$Q_plus_2q <- k2_q_Q_combined$X2 + (2*k2_q_Q_combined$V1)
# 
# k2_q_Q_combined$classification <- NA
# 
# k2_q_Q_combined$classification[ (k2_q_Q_combined$V1 > 0.9) & (k2_q_Q_combined$X2 <= 0.25)] <- "polyodon"
# k2_q_Q_combined$classification[ (k2_q_Q_combined$V1 < 0.1) & (k2_q_Q_combined$X2 <= 0.25)] <- "kazumbe"
# k2_q_Q_combined$classification[ (k2_q_Q_combined$V1 > 0.4) & (k2_q_Q_combined$V1 < 0.6) & (k2_q_Q_combined$X2 > 0.8)] <- "F1"
# k2_q_Q_combined$classification[ (k2_q_Q_combined$V1 > 0.4) & (k2_q_Q_combined$V1 < 0.6) & (k2_q_Q_combined$X2 > 0.4) & (k2_q_Q_combined$X2 < 0.6)] <- "F2"
# 
# k2_q_Q_combined$classification[is.na(k2_q_Q_combined$classification) & (k2_q_Q_combined$Q_minus_2q > -0.1) & (k2_q_Q_combined$Q_minus_2q < 0.1)] <- "bc.kazumbe"
# k2_q_Q_combined$classification[is.na(k2_q_Q_combined$classification) & (k2_q_Q_combined$Q_plus_2q > 1.9) & (k2_q_Q_combined$Q_plus_2q < 2.1)] <- "bc.polyodon"
# k2_q_Q_combined$classification[ is.na(k2_q_Q_combined$classification) ] <- "other"
# k2_q_Q_combined$classification <- factor(k2_q_Q_combined$classification, levels = c("kazumbe", "polyodon", "F1", "F2", "bc.kazumbe", "bc.polyodon", "other"))
# 
# 
# location_sample_size <- as.data.frame(table(droplevels(k2_q_Q_combined$location_plotting)))


# #to add custom text to (i.e. the sample sizes of each location placed above each bar), you need to put the data and aes in the geom_bar() function, not the ggplot() function
# hybrid_distribution_barchart <- ggplot() + 
#   geom_bar(data = k2_q_Q_combined, aes(x = location_plotting, fill = classification), position = "fill") +
#   scale_fill_manual(values = c("#DC7633", "#3498DB", "#edba99", "#99cbed", "#66b266")) +
#   theme(#legend.position = "none", 
#     plot.background = element_rect(fill = "white"), 
#     panel.background = element_rect(fill = "white"),
#     axis.text.y=element_text(size = 18, color = 'black'),
#     #axis.text.x=element_text(angle = 45, margin=margin(25,0, -40,0), hjust=1, size = 14.4, color = 'black', vjust = 1.6),
#     axis.text.x=element_text(hjust=0.5, size = 18, color = 'black', vjust = 3.5),
#     axis.ticks.x=element_blank(),
#     axis.title.x = element_blank(),
#     axis.title.y = element_text(size = 20, face = "plain"),
#     legend.text = element_text(size = 15),
#     legend.title = element_blank(),
#     plot.margin=unit(c(3,0,0,0),"mm")) +
#   ylab("Proportion of individuals") +
#   geom_text(aes(x = Var1, y = 1, label = Freq), size = 6, data = location_sample_size, vjust = -0.4) +
#   scale_y_continuous(breaks = seq(0, 1, 0.25), labels = c("0", "0.25", "0.5", "0.75", "1")) +
#   coord_cartesian(clip = "off")
# 
# 
# Q_q_full_plot <- ggplot(data = k2_q_Q_combined, aes(x = V1, y = X2, color = location_plotting) ) +
#   geom_line(data = triangle_dataframe, aes(x = x, y = y), color = '#bbbaba', size = 1.5, lineend = "round") +
#   geom_point(size = 6.5, alpha = 0.5) +
#   geom_point(size = 6.5, shape = 1, alpha = 0.5) +
#   scale_color_manual(name = "Location", values = location_palette_letters[names(location_palette_letters) %in% unique(k2_q_Q_combined$location_plotting) ]) +
#   scale_x_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
#   scale_y_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
#   xlab("q (proportion of ancestry)") +
#   ylab(expression(Q[12]~"(interspecific ancestry)") ) +
#   theme_Q_entropy_alt() +
#   theme(plot.title = element_text(hjust = 0.5, size = 16),
#         axis.title = element_text(size = 16, face = "plain"),
#         axis.text = element_text(size = 18),
#         plot.margin=unit(c(5,2,1,2),"mm"))
# 
# 
# upper_plot_Qq_bar<- cowplot::plot_grid(Q_q_full_plot, hybrid_distribution_barchart, nrow = 1, labels = c("(a)", "(b)"), label_size = 22, scale = 0.9, hjust = c(-0.5, 0))
# 
# Qq_multipanel_FULL_withbarchart <- cowplot::plot_grid(upper_plot_Qq_bar,
#                                                       ggdraw(k2_Q_q_multipanel_lower), nrow = 2, rel_heights = c(2,2.5), 
#                                                       labels = c("", "(c)"), 
#                                                       label_size = 22)
# 
# Qq_multipanel_FULL_withbarchart
# ggsave2("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/entropy/results/results_5_19_2021/plots/k2_Q_q_multipanel_FULL_withbarchart_5_19_2021.png", width = 20, height = 12)




# Q_q_full_plot_alt <- ggplot(data = k2_q_Q_combined, aes(x = V1, y = X2, color = location_plotting) ) +
#   geom_line(data = triangle_dataframe, aes(x = x, y = y), color = '#bbbaba', size = 1.5, lineend = "round") +
#   geom_point(size = 6.5, alpha = 0.5) +
#   geom_point(size = 6.5, shape = 1, alpha = 0.5) +
#   scale_color_manual(name = "Location", values = location_palette_letters[names(location_palette_letters) %in% unique(k2_q_Q_combined$location_plotting) ]) +
#   scale_x_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
#   scale_y_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
#   xlab("q (proportion of ancestry)") +
#   ylab(expression(Q[12]~"(interspecific ancestry)") ) +
#   theme_Q_entropy_alt() +
#   theme(plot.title = element_text(hjust = 0.5, size = 16),
#         axis.title = element_text(size = 20, face = "plain"),
#         axis.text = element_text(size = 16),
#         plot.margin=unit(c(5, 1.5, 1, 2),"mm"),
#         legend.position = 'none')
# 
# legend_Q_q_full_plot_alt <- get_legend(ggplot(data = k2_q_Q_combined, aes(x = V1, y = X2, color = location_plotting) ) +
#                                          geom_line(data = triangle_dataframe, aes(x = x, y = y), color = '#bbbaba', size = 1.5, lineend = "round") +
#                                          geom_point(size = 6.5, alpha = 0.5) +
#                                          geom_point(size = 6.5, shape = 1, alpha = 0.5) +
#                                          scale_color_manual(name = "Location", values = location_palette_letters[names(location_palette_letters) %in% unique(k2_q_Q_combined$location_plotting) ]) +
#                                          scale_x_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
#                                          scale_y_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
#                                          xlab("q (proportion of ancestry)") +
#                                          ylab(expression("Q[12] (interspecific ancestry)") ) +
#                                          theme_Q_entropy_alt() +
#                                          theme(plot.title = element_text(hjust = 0.5, size = 16),
#                                                axis.title = element_text(size = 16, face = "plain"),
#                                                legend.text = element_text(size = 18),
#                                                legend.title = element_text(size = 20),
#                                                plot.margin=unit(c(5, 0, 1, 4),"mm")))
# 
# Q_q_full_plot_alt_marg_hist <- ggMarginal(Q_q_full_plot_alt, type = "histogram",
#                                           xparams = list(binwidth = 0.05, fill = '#bbbaba', col = 'black'),
#                                           yparams = list(binwidth = 0.05, fill = '#bbbaba', col = 'black'))
# 
# Q_q_full_plot_alt_marg_hist_withlegend <- cowplot::plot_grid(Q_q_full_plot_alt_marg_hist, ggdraw(legend_Q_q_full_plot_alt) + theme(plot.margin=unit(c(0, 0, 0, 1),"cm")), 
#                                                              rel_widths = c(2.6, 0.5))
# 
# 
# k2_Q_q_multipanel_FULL_withmarg <- cowplot::plot_grid(cowplot::plot_grid(NULL, Q_q_full_plot_alt_marg_hist_withlegend, NULL, rel_widths = c(1.3, 3, 0.7), nrow = 1, labels = c("", "(a)", ""), 
#                                                                          label_size = 25, hjust = 1.2), 
#                                                       ggdraw(k2_Q_q_multipanel_lower), nrow = 2, rel_heights = c(2.5, 3.5), 
#                                                       labels = c("", "(b)"), 
#                                                       label_size = 25)
# 
# upper_plot_Qq_bar_withmarg <- cowplot::plot_grid(Q_q_full_plot_alt_marg_hist_withlegend, hybrid_distribution_barchart, 
#                                                  nrow = 1, labels = c("(a)", "(b)"), label_size = 25, scale = 0.9, hjust = c(-0.5, 0))
# 
# 
# Qq_multipanel_FULL_withbarchart_withmarg <- cowplot::plot_grid(upper_plot_Qq_bar_withmarg,
#                                                                ggdraw(k2_Q_q_multipanel_lower), nrow = 2, rel_heights = c(2, 2.6), 
#                                                                labels = c("", "(c)"), 
#                                                                label_size = 25)
# 
# 
# k2_Q_q_multipanel_FULL_withmarg
# ggsave2("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/entropy/results/results_5_19_2021/plots/k2_Q_q_multipanel_FULL_marginal_5_19_2021.png", width = 20, height = 12)
# 
# Qq_multipanel_FULL_withbarchart_withmarg
# ggsave2("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/entropy/results/results_5_19_2021/plots/k2_Q_q_multipanel_FULL_barchart_marginal_5_19_2021.png", width = 20, height = 12)
# 
# 
# entropy_list$K_2$q_long
# 
# entropy_list <- list(K2 = k2_q_long,
#                      K3 = k3_q_long,
#                      K4 = k4_q_long,
#                      K5 = k5_q_long,
#                      K6 = k6_q_long)
# 
# V2_plots <- create_V2_plots(lapply(entropy_list, function(x) x$q_long) , q_col_pal = q_color_palette)
# V2_plots$plot_list$K_2location_genetic
# 
# title_size <- 25
# vertical_adjust_title <- -1
# 
# 
# entropy_bar_plot_multi <- plot_grid(V2_plots$plot_list$K_2location_genetic + ggtitle("") +
#                                       #ggtitle(paste0("K = 2; DIC = ", formatC(DIC_dataframe$DIC[2], format="e", digits = 2)) ) +
#                                       ggtitle(expression(italic('K') ~ " = 2")) +
#                                       theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size),
#                                             plot.margin=unit(c(0.1, 0.35, 0.5, 0),"cm")) , 
#                                     V2_plots$plot_list$K_3location_genetic + ggtitle("") +
#                                       #ggtitle(paste0("K = 3; DIC = ", formatC(DIC_dataframe$DIC[3], format="e", digits = 2)) ) +
#                                       ggtitle(expression(italic('K') ~ " = 3")) +
#                                       theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size),
#                                             plot.margin=unit(c(0.1, 0, 0.5, 0.35),"cm")), ncol = 2)
# 
# plot_grid(Qq_multipanel_FULL_withbarchart_withmarg, entropy_bar_plot_multi + theme(plot.margin=unit(c(0.85, 0.7, 0.3, 0.8),"cm")), 
#           labels = c("", "(d)"), label_y = c(1, 1.01) ,nrow = 2, label_size = 25, rel_heights = c(0.8, 0.25))
# 
# ggsave2("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/entropy/results/results_2_12_2021/plots/k2_Q_q_multipanel_FULL_barchart_marginal_qbarplot_5_25_2021.png", 
#         width = 19, height = 14)


### q plots ###

# #list of entropy results dataframe with sequentially increasing K vals
# #dataframes contain metadata (most importantly location, sample code)
# entropy_list <- list(K2 = k2_q_long,
#                      K3 = k3_q_long,
#                      K4 = k4_q_long,
#                      K5 = k5_q_long,
#                      K6 = k6_q_long)
# 
# 
# ### Create list of V1 and V2 plots
# #show_progress = FALSE if you don't want functions to print plot creation progress
# V1_plots <- create_V1_plots(entropy_list, q_col_pal = q_color_palette)
# 
# V2_plots <- create_V2_plots(entropy_list, q_col_pal = q_color_palette)
# 
# 
# ### Create multipanel plots ###
# #Visualization parameters for multipanel plots
# title_size <- 23
# vertical_adjust_title <- -1
# 
# #V1 multipanel plot
# V1_multipanel <- plot_grid(
#   plot_grid(V1_plots$plot_list$K_2genetic +
#               ggtitle(paste0("k = 2; DIC = ", formatC(DIC_dataframe$DIC[2], format="e", digits = 2)) ) +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)) , 
#             V1_plots$plot_list$K_2location_genetic + ggtitle("") +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)), ncol = 2),
#   plot_grid(V1_plots$plot_list$K_3genetic +
#               ggtitle(paste0("k = 3; DIC = ", formatC(DIC_dataframe$DIC[3], format="e", digits = 2)) ) +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)) , 
#             V1_plots$plot_list$K_3location_genetic + ggtitle("") +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)), ncol = 2),
#   plot_grid(V1_plots$plot_list$K_4genetic +
#               ggtitle(paste0("k = 4; DIC = ", formatC(DIC_dataframe$DIC[4], format="e", digits = 2)) ) +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)) , 
#             V1_plots$plot_list$K_4location_genetic + ggtitle("") +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)), ncol = 2),
#   plot_grid(V1_plots$plot_list$K_5genetic +
#               ggtitle(paste0("k = 5; DIC = ", formatC(DIC_dataframe$DIC[5], format="e", digits = 2)) ) +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)) , 
#             V1_plots$plot_list$K_5location_genetic + ggtitle("") +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)), ncol = 2),
#   plot_grid(V1_plots$plot_list$K_6genetic +
#               ggtitle(paste0("k = 6; DIC = ", formatC(DIC_dataframe$DIC[6], format="e", digits = 2)) ) +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)) , 
#             V1_plots$plot_list$K_6location_genetic + ggtitle("") +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)), ncol = 2),
#   nrow = 5
# )
# 
# #V2 multipanel plot
# V2_multipanel <- plot_grid(
#   plot_grid(V2_plots$plot_list$K_2genetic +
#               ggtitle(paste0("k = 2; DIC = ", formatC(DIC_dataframe$DIC[2], format="e", digits = 2)) ) +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)) , 
#             V2_plots$plot_list$K_2location_genetic + ggtitle("") +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)), ncol = 2),
#   plot_grid(V2_plots$plot_list$K_3genetic +
#               ggtitle(paste0("k = 3; DIC = ", formatC(DIC_dataframe$DIC[3], format="e", digits = 2)) ) +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)) , 
#             V2_plots$plot_list$K_3location_genetic + ggtitle("") +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)), ncol = 2),
#   plot_grid(V2_plots$plot_list$K_4genetic +
#               ggtitle(paste0("k = 4; DIC = ", formatC(DIC_dataframe$DIC[4], format="e", digits = 2)) ) +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)) , 
#             V2_plots$plot_list$K_4location_genetic + ggtitle("") +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)), ncol = 2),
#   plot_grid(V2_plots$plot_list$K_5genetic +
#               ggtitle(paste0("k = 5; DIC = ", formatC(DIC_dataframe$DIC[5], format="e", digits = 2)) ) +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)) , 
#             V2_plots$plot_list$K_5location_genetic + ggtitle("") +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)), ncol = 2),
#   plot_grid(V2_plots$plot_list$K_6genetic +
#               ggtitle(paste0("k = 6; DIC = ", formatC(DIC_dataframe$DIC[6], format="e", digits = 2)) ) +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)) , 
#             V2_plots$plot_list$K_6location_genetic + ggtitle("") +
#               theme(plot.title = element_text(vjust = vertical_adjust_title, size = title_size)), ncol = 2),
#   nrow = 5
# )
# 
# #create the y axis title for the multipanel plots
# qplot_y.grob <- textGrob("Proportion of ancestry", 
#                          gp=gpar(fontface="plain", col="black", fontsize=30), rot=90)
# 
# 
# ### V1 plot ###
# #add y axis title to multipanel plot
# V1_multipanel_full <- ggdraw(grid.arrange(arrangeGrob(V1_multipanel, left = qplot_y.grob), vp=viewport(width=0.99, height=0.99)))
# 
# #export V1 multipanel plot
# V1_multipanel_full
# ggsave2("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/entropy/results/results_8_27_2020/plots/q_plot_multipanel_V1_8_31_2020.png", width = 21, height = 16)
# #ggsave2("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/entropy/results/results_8_27_2020/plots/k2_Q_q_multipanel_FULL_withbarchart_8_27_2020.png", width = 20, height = 12)
# 
# ### V2 plot ###
# #add y axis title to multipanel plot
# V2_multipanel_full <- ggdraw(grid.arrange(arrangeGrob(V2_multipanel, left = qplot_y.grob), vp=viewport(width=0.99, height=0.99)))
# 
# #export V2 multipanel plot
# V2_multipanel_full
# ggsave2("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/entropy/results/results_8_27_2020/plots/q_plot_multipanel_V2_8_31_2020.png", width = 21, height = 16)

# 
# mid_q_ci <- t(as.data.frame(apply(k2_q_Q_list$mid$q$combined.q, 2:3, quantile, probs = c(0.025,0.975))))
# 
# apply(k2_q_Q_list$mid$q$combined.q, 2:3, quantile, probs = c(0.025,0.975))
# 
# 
# mid_Q12_ci <- apply(k2_q_Q_list$mid$Q$combined.Q, 2:3, quantile, probs = c(0.025,0.975))[,,14]
# 
# 
# t(apply(k2_q_Q_list$mid$Q$combined.Q, 2:3, mean))
# 
# do.call(mid_Q12_ci, rbind)
# 
# 
# do.call(rbind, lapply(seq(dim(mid_Q12_ci)[3]), function(x) mid_Q12_ci[ , , x])) %>% 
#   as.data.frame() %>% 
#   rownames_to_column() %>%
#   mutate(Q_ci = paste0('Qci_', sub(pattern = 'X([0-9]+\\.[0-9]).*', replacement = '\\1', x = rowname)),
#          index = rep(1:(n()/2), each = 2)) %>%
#   select(!c(rowname, V1, V3) ) %>% 
#   pivot_wider(names_from = 'Q_ci', values_from = 'V2') %>%
#   cbind(., t(as.data.frame(mid_q_ci)) %>% 
#           as.data.frame() %>% 
#           rename(qci_2.5 = `2.5%`,
#                  qci_97.5 = `97.5%`)) %>% 
#   ggplot() +
#   geom_segment(aes(x = qci_2.5, xend = qci_97.5, y = Qci_2.5, yend = Qci_97.5))
# 
# sub(pattern = 'X([0-9]+\\.[0-9]).*', replacement = '\\1', x = 'X97.5..75')
# 
# data.frame(t(mid_Q12_ci))
# 
# data.frame(mean(mid_Q12_ci[2,2,] - mid_Q12_ci[1,2,]))
# 
# 
# mid_q_ci <- apply(k2_q_Q_list$mid$q$combined.q, 2:3, quantile, probs = c(0.025,0.975))
# 
# mean(mid_q_ci[2,2,] - mid_q_ci[1,2,])


# rbind(cbind(k2_q_Q_list$north$q$q_estimate_with_metadata_and_ci, k2_q_Q_list$north$Q$Q_estimate),
#       cbind(k2_q_Q_list$mid$q$q_estimate_with_metadata_and_ci, k2_q_Q_list$mid$Q$Q_estimate)) %>% 
#   ggplot(aes_string(x = focal_var, y = "X2", color = "location_plotting") ) +
#   geom_line(data = triangle_dataframe, aes(x = x, y = y), color = '#bbbaba', size = 1.5, lineend = "round") +
#   geom_point(size = 6.5, alpha = 0.4) +
#   geom_point(size = 6.5, shape = 1, alpha = 0.4) +
#   scale_color_manual(name = "Location", values = location_palette_letters[names(location_palette_letters) %in% unique(k2_q_Q_combined$location_plotting) ]) +
#   scale_x_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
#   scale_y_continuous(breaks = seq(0, 1, 0.5), labels = c("0", "0.5", "1")) +
#   xlab("q (proportion of ancestry)") +
#   ylab(expression(Q[12]~"(interspecific ancestry)") ) +
#   theme_Q_entropy_alt() +
#   theme(plot.title = element_text(hjust = 0.5, size = 16),
#         axis.title = element_text(size = 20, face = "plain"),
#         axis.text = element_text(size = 16),
#         plot.margin=unit(c(5, 1.5, 1, 2),"mm"),
#         legend.position = 'none') +
#   geom_segment(data = rbind(cbind(k2_q_Q_list$north$q$q_estimate_with_metadata_and_ci, k2_q_Q_list$north$Q$Q_estimate),
#                             cbind(k2_q_Q_list$mid$q$q_estimate_with_metadata_and_ci, k2_q_Q_list$mid$Q$Q_estimate)),
#                aes(x = qci_2.5, xend = qci_97.5, y = X2, yend = X2), size = 1) +
#   geom_segment(data = rbind(cbind(k2_q_Q_list$north$q$q_estimate_with_metadata_and_ci, k2_q_Q_list$north$Q$Q_estimate),
#                             cbind(k2_q_Q_list$mid$q$q_estimate_with_metadata_and_ci, k2_q_Q_list$mid$Q$Q_estimate)),
#                aes(y = Qci_2.5, yend = Qci_97.5, x = V2, xend = V2), size = 1)
# 

