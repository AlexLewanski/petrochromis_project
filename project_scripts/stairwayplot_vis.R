#################################################################################################
### SCRIPT NAME: stairwayplot_vis.R
### PURPOSE: Visualize estimates of historical population size Stairway Plot. These plots are 
###          included in the supplementary material.
### PRODUCT:
###     stairway_plot_results_9_20_2021.png: plot of historical population size estimated with
###                                          Stairway Plot. Population size reconstructed
###                                          separately for each species in each region (e.g.
###                                          kazumbe from the north region).
#################################################################################################


#####################
### SCRIPT SET-UP ###
#####################

### Loading libraries ###
library(tidyverse)
library(cowplot)
library(scales)
library(here)


### Custom plotting function ###
plot_stairway <- function(stairway_output, 
                          years_per_gen = 1,
                          median_col = "red",
                          CI95_col = '#ced4da',
                          CI75_col = '#adb5bd',
                          plot_title = 'Stairway results',
                          title_size = 15,
                          axis_title_size = 13) {
  
  stairway_output %>%
    mutate(generation = year/years_per_gen) %>% 
    ggplot() +
    geom_ribbon(aes(x = generation, ymin = Ne_2.5., ymax = Ne_97.5.), col = NA, fill = CI95_col, alpha = 0.5) +
    geom_ribbon(aes(x = generation, ymin = Ne_12.5., ymax = Ne_87.5.), col = NA, fill = CI75_col, alpha = 0.5) +
    geom_line(aes(x = generation, y = Ne_median), col = median_col, size = 2) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    annotation_logticks() +
    ggtitle(plot_title) +
    xlab('Generation') + ylab('Ne') +
    theme_bw() +
    theme(plot.title = element_text(size = title_size),
          axis.title = element_text(size = axis_title_size))
}


### Loading and processing data ###
#***the code for creating the path stored in results_full_path may need to be adjusted for non-mac OS
results_files <- list.files(here('demographic_modelling', 'stairway', 'results', 'stairway_9_3_2021'))[grep('final.summary$', list.files(here('demographic_modelling', 'stairway', 'results', 'stairway_9_3_2021')))]
results_full_path <- setNames(paste0(here('demographic_modelling', 'stairway', 'results', 'stairway_9_3_2021'), '/', results_files), 
                              nm = gsub("_MAF\\.final\\.summary", "", results_files))
stairway_results <- lapply(results_full_path, function(FILE) read.table(FILE, header = TRUE))



####################################
### VISUALIZING STAIRWAY RESULTS ###
####################################

stairway_results_processed <- lapply(setNames(nm = names(stairway_results)), function(NAME) {
  
  stairway_results[[NAME]] %>% 
    mutate(species = gsub("^.*_([a-z]*)", "\\1", NAME),
           region = gsub("^([a-z]*)_.*", "\\1", NAME))
  
})


stairplot_plots_list <- lapply(stairway_results_processed, function(DATA) {
  
  plot_stairway(stairway_output = DATA, 
                years_per_gen = 2,
                median_col = '#1697a6',
                CI95_col = '#ced4da',
                CI75_col = '#adb5bd',
                plot_title = paste0('Species: ', DATA$species[1], '; Region: ', DATA$region[1]),
                title_size = 20,
                axis_title_size = 17)
  
})

plot_grid(stairplot_plots_list$north_region_kazumbe, stairplot_plots_list$north_region_polyodon,
          stairplot_plots_list$mid_region_kazumbe, stairplot_plots_list$mid_region_polyodon,
          nrow = 2)

ggsave2(here('figures', 'stairway_plot_results_9_20_2021.png'), 
        width = 16, height = 10, bg = 'white')

