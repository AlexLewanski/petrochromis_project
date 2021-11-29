##############################################################################################
### SCRIPT NAME: d_stat_analysis_and_vis.R
### PURPOSE: visualize and analyze D statistics that quantified variation in interspecific
###          allele sharing between different populations of kazumbe (1a and 1b) and between
###          different populations of polyodon (2a and 2b)
### PRODUCTS:
###     dstat_matrix_plot_typea_outgroup_diagramma_9_19_2021_REARRANGE.png: visualization for
###                                    D stats with the following setups using Simochromis
###                                    diagramma as the outgroup:
###                                    1b: kazumbe_pop1, kazumbe_pop2, polyodon_pop2;
###                                    2b: polyodon_pop1, polyodon_pop2, kazumbe_pop2
###     dstat_matrix_plot_typea_outgroup_green_9_19_2021_REARRANGE.png: visualization for D
###                                    stats with the following setups using Petrochromis
###                                    green as the outgroup:
###                                    1b: kazumbe_pop1, kazumbe_pop2, polyodon_pop2;
###                                    2b: polyodon_pop1, polyodon_pop2, kazumbe_pop2
###     dstat_matrix_plot_typeb_outgroup_diagramma_9_19_2021_REARRANGE.png: visualization for
###                                    D stats with the following setups using Simochromis
###                                    diagramma as the outgroup:
###                                    1a: kazumbe_pop1, kazumbe_pop2, polyodon_all
###                                    2a: polyodon_pop1, polyodon_pop2, kazumbe_all
###     dstat_matrix_plot_typeb_outgroup_green_9_19_2021_REARRANGE.png: visualization for D 
###                                    stats with the following setups using Petrochromis
###                                    green as the outgroup:
###                                    1a: kazumbe_pop1, kazumbe_pop2, polyodon_all
###                                    2a: polyodon_pop1, polyodon_pop2, kazumbe_all
###     kazumbe_dstat_table.txt: table of all D statistics that are focused on variation
###                              in interspecific allele sharing between populations of
###                              kazumbe (1a and 1b).
###     polyodon_dstat_table.txt: table of all D statistics that are focused on variation
###                               in interspecific allele sharing between populations of
###                               polyodon (2a and 2b).
###     dstat_mantel_test_table.txt: table of results from Mantel tests quantifying the
###                                  degree of correspondence between various sets of
###                                  D statistics.
##############################################################################################



#############################
### 1. SCRIPT PREPARATION ###
#############################

### NOTES: ###
#TYPE B D statistics have the following topology: ((P1 = sp1_pop1, P2 = sp1_pop2), P3 = sp2_pop2), P4 = outgroup)
#TYPE A D statistics have the following topology: ((P1 = sp1_pop1, P2 = sp1_pop2), P3 = sp2_full), P4 = outgroup)
#this is different than in the manuscript where the a topology refers to the type b d stat structure

#the type b d stats here refers to the typa a topology in the manuscript
#confusingly, the legends are currently based on the manuscript, so the type a d statistics have the "b" legend and vice versa

#location name/letter information
# 'Gombe_South' ~ 'A'
# 'Katongwe_N' ~ 'B'
# 'Katongwe_S' ~ 'C'
# 'Ska' ~ 'D'
# 'Kalala' ~ 'E'
# 'Nondwa' ~ 'F'
# 'Hilltop' ~ 'G'
# 'Bangwe' ~ 'H'
# 'Jakob' ~ 'I'
# 'Ulombola' ~ 'J'


### Color scheme
#first column of colors is the original color scheme (red to green) --> not friendly for color blindness ###
#new color scheme (purple to green) --> more friendly for color blindness
#"group1" = "#EA1313" --> #7b3294 --> #88469e
#"group2" =  '#f8b8b8' --> #c2a5cf
#"group3" =  "#a6a6a6" --> #a6a6a6
#"group4" =  "#99e099" --> #a6dba0
#"group5" =  "#00b300" --> #008837 --> #19934b


### RESOURCES ###
#plotting
#https://stackoverflow.com/questions/20251119/increase-the-size-of-variable-size-points-in-ggplot2-scatter-plot
#https://stackoverflow.com/questions/32555531/how-to-italicize-part-one-or-two-words-of-an-axis-title
#https://rpubs.com/lgadar/matrix-visualizations
#https://stackoverflow.com/questions/20035359/diverging-size-scale-ggplot2

#latex tables
#https://stackoverflow.com/questions/45206908/kableextra-dynamic-add-header-above-labeling


### Load libraries ###
library(tidyverse)
library(vegan) #--> used for Mantel tests
library(reshape2)
library(cowplot)
library(egg)
library(kableExtra)
library(here)


### Load D statistic results ###
dstat_type_a <- read.table(here('DStats', 'Results', 'results_9_19_2021', 'dstat_results_TYPE_A_colormorphs_9_19_21.txt'), header = TRUE)
dstat_type_b_kazumbe <- read.table(here('DStats', 'Results', 'results_9_19_2021', 'dstat_results_TYPEBKazumbe_colormorphs_9_19_21.txt'), header = TRUE)
dstat_type_b_polyodon <- read.table(here('DStats', 'Results', 'results_9_19_2021', 'dstat_results_TYPEBPolyodon_colormorphs_9_19_21.txt'), header = TRUE)



#######################
### 2. PROCESS DATA ###
#######################

### (ADMIXTOOLS) 2a. create lists of the d-stats ###
type_a_list <- list(kazumbe = dstat_type_a[which(gsub("_[A-za-z]*", '', dstat_type_a$P1) == 'kazumbe' ),], 
                    polyodon = dstat_type_a[which(gsub("_[A-za-z]*", '', dstat_type_a$P1) == 'polyodon' ),])

type_b_list <- list(kazumbe = dstat_type_b_kazumbe, polyodon = dstat_type_b_polyodon)

full_data_list <- list(type_a = type_a_list,
                       type_b = type_b_list)


### (ADMIXTOOLS) 2b. initial processing of d-stat results including extracting taxa and location information for p1, p2, and p3 (ADMIXTOOLS) ###
for (type in names(full_data_list)) {
  
  for (i in names(full_data_list[[type]])) {
    full_data_list[[type]][[i]]$p1_p2_taxa <- gsub("_.*" , "", full_data_list[[type]][[i]]$P1)
    
    full_data_list[[type]][[i]]$p1_location <- str_remove(full_data_list[[type]][[i]]$P1, "[A-Za-z]*_")
    
    full_data_list[[type]][[i]]$p1_location_factor <- factor(full_data_list[[type]][[i]]$p1_location, 
                                                             levels = c("Gombe_South", "Katongwe_N", "Katongwe_S", "Ska", "Kalala", "Nondwa", "Hilltop", "Jakob"))
    
    full_data_list[[type]][[i]]$p2_p3_location <- str_remove(full_data_list[[type]][[i]]$P2, "[A-Za-z]*_")
    
    full_data_list[[type]][[i]]$p2_p3_location_factor <- factor(full_data_list[[type]][[i]]$p2_p3_location, 
                                                                levels = c("Gombe_South", "Katongwe_N", "Katongwe_S", "Ska", "Kalala", "Nondwa", "Hilltop", "Jakob"))
    
    full_data_list[[type]][[i]]$p3_taxa <- gsub("_.*" , "", full_data_list[[type]][[i]]$P3)
    
    
    full_data_list[[type]][[i]] <- full_data_list[[type]][[i]] %>% 
      mutate(p1_location_factor_plotting = factor(case_when(p1_location_factor == 'Gombe_South' ~ "Gombe S",
                                                            p1_location_factor == 'Katongwe_N' ~ "Katongwe N",
                                                            p1_location_factor == 'Katongwe_S' ~ "Katongwe S",
                                                            p1_location_factor == 'Ska' ~ "S Kagango",
                                                            p1_location_factor == 'Kalala' ~ "Kalalangabo",
                                                            p1_location_factor == 'Nondwa' ~ "Nondwa",
                                                            p1_location_factor == 'Hilltop' ~ "Hilltop",
                                                            p1_location_factor == 'Bangwe' ~ "Bangwe",
                                                            p1_location_factor == 'Jakob' ~ "Jakobsen's S",
                                                            p1_location_factor == 'Ulombola' ~ "Ulombola"),
                                                  levels= c("Gombe S", "Katongwe N", "Katongwe S", "S Kagango", "Kalalangabo", "Nondwa", "Hilltop", "Bangwe", "Jakobsen's S", "Ulombola")),
             p1_location_factor_letter = factor(case_when(p1_location_factor == 'Gombe_South' ~ "A",
                                                          p1_location_factor == 'Katongwe_N' ~ "B",
                                                          p1_location_factor == 'Katongwe_S' ~ "C",
                                                          p1_location_factor == 'Ska' ~ "D",
                                                          p1_location_factor == 'Kalala' ~ "E",
                                                          p1_location_factor == 'Nondwa' ~ "F",
                                                          p1_location_factor == 'Hilltop' ~ "G",
                                                          p1_location_factor == 'Bangwe' ~ "H",
                                                          p1_location_factor == 'Jakob' ~ "I",
                                                          p1_location_factor == 'Ulombola' ~ "J"),
                                                levels= LETTERS[1:10]),
             p2_p3_location_factor_plotting = factor(case_when(p2_p3_location_factor == 'Gombe_South' ~ "Gombe S",
                                                               p2_p3_location_factor == 'Katongwe_N' ~ "Katongwe N",
                                                               p2_p3_location_factor == 'Katongwe_S' ~ "Katongwe S",
                                                               p2_p3_location_factor == 'Ska' ~ "S Kagango",
                                                               p2_p3_location_factor == 'Kalala' ~ "Kalalangabo",
                                                               p2_p3_location_factor == 'Nondwa' ~ "Nondwa",
                                                               p2_p3_location_factor == 'Hilltop' ~ "Hilltop",
                                                               p2_p3_location_factor == 'Bangwe' ~ "Bangwe",
                                                               p2_p3_location_factor == 'Jakob' ~ "Jakobsen's S",
                                                               p2_p3_location_factor == 'Ulombola' ~ "Ulombola"),
                                                     levels= c("Gombe S", "Katongwe N", "Katongwe S", "S Kagango", "Kalalangabo", "Nondwa", "Hilltop", "Bangwe", "Jakobsen's S", "Ulombola")),
             p2_p3_location_factor_letter = factor(case_when(p2_p3_location_factor == 'Gombe_South' ~ "A",
                                                             p2_p3_location_factor == 'Katongwe_N' ~ "B",
                                                             p2_p3_location_factor == 'Katongwe_S' ~ "C",
                                                             p2_p3_location_factor == 'Ska' ~ "D",
                                                             p2_p3_location_factor == 'Kalala' ~ "E",
                                                             p2_p3_location_factor == 'Nondwa' ~ "F",
                                                             p2_p3_location_factor == 'Hilltop' ~ "G",
                                                             p2_p3_location_factor == 'Bangwe' ~ "H",
                                                             p2_p3_location_factor == 'Jakob' ~ "I",
                                                             p2_p3_location_factor == 'Ulombola' ~ "J"),
                                                   levels= LETTERS[1:10]))
    
    full_data_list[[type]][[i]]$summarized_inference <- NA
    full_data_list[[type]][[i]]$summarized_inference[full_data_list[[type]][[i]]$z_score < -3] <- "admixture"
    full_data_list[[type]][[i]]$summarized_inference[full_data_list[[type]][[i]]$z_score > 3] <- "alt. population admixture"
    full_data_list[[type]][[i]]$summarized_inference[is.na(full_data_list[[type]][[i]]$summarized_inference)] <- "none"
  }
}


### (ADMIXTOOLS) 2c. further processing of results -- creating a new hierarchical list of dataframes that include D-stat interpretations ###
outgroup_taxa <- c('diagramma', 'green') #the outgroup taxa used in the d-stats
focal_taxa <- c('kazumbe', 'polyodon') #the two species included in the d-stats

dstat_df_list <- list()

for (type in names(full_data_list)) {
  outgroup_list <- list()
  
  for (out in outgroup_taxa) {
    focal_taxa_list <- list()
    
    for (focal in focal_taxa) {
      focal_df <- full_data_list[[type]][[focal]][which(full_data_list[[type]][[focal]]$outgroup == out), c('p2_p3_location_factor_plotting', 'p2_p3_location_factor_letter', 'p1_location_factor_plotting', 'p1_location_factor_letter', 'D', 'z_score')]
      focal_df$D_absolute <- abs(focal_df$D)
      focal_df$z_score_absolute <- abs(focal_df$z_score)
      
      focal_df$inference[focal_df$z_score < - 3] <- "P2 allele sharing"
      focal_df$inference[focal_df$z_score > 3] <- "P1 allele sharing"
      focal_df$inference[ (focal_df$z_score <= 3) & (focal_df$z_score >= -3)] <- "no excess allele sharing"
      
      #c2a5cf light red
      ##a6dba0 light green
      focal_df <- focal_df %>% 
        mutate(color_z = case_when(z_score < -3 ~ "#88469e",
                                   z_score > 3 ~ "#19934b",
                                   z_score >= -3 & z_score < 0 ~ "#c2a5cf",
                                   z_score <= 3 & z_score > 0 ~ "#a6dba0",
                                   z_score == 0 ~ "#a6a6a6") )
      
      focal_taxa_list[[focal]] <- focal_df
    }
    outgroup_list[[out]] <- focal_taxa_list
  }
  dstat_df_list[[type]] <- outgroup_list
}


### (ADMIXTOOLS) 2d. create a hierarchical list of matrices of D-stat values and D-stat z-scores (used in Mantel tests) ###
dstat_matrix_list <- list()

for (type in names(dstat_df_list)) {
  outgroup_matrix_list <- list()
  
  for (out in names(dstat_df_list[[type]]) ) {
    focal_taxa_matrix_list <- list()
    
    for (focal in names(dstat_df_list[[type]][[out]])) {
      dstat_mat <- acast(dstat_df_list[[type]][[out]][[focal]], p2_p3_location_factor_plotting~p1_location_factor_plotting, value.var = "D")
      dstat_mat[is.na(dstat_mat)] <- 0
      
      zscore <- acast(dstat_df_list[[type]][[out]][[focal]], p2_p3_location_factor_plotting~p1_location_factor_plotting, value.var = "z_score")
      zscore[is.na(zscore)] <- 0
      
      focal_taxa_matrix_list[[focal]] <- list(d_stat_val = dstat_mat, 
                                              d_stat_zscore = zscore)
    }
    outgroup_matrix_list[[out]] <- focal_taxa_matrix_list
  }
  dstat_matrix_list[[type]] <- outgroup_matrix_list
}


significant_dstats <- sapply(setNames(nm = c('kazumbe', 'polyodon')), function(SPECIES, results) {
  processed1 <- results[[SPECIES]] %>% 
    filter(outgroup == 'diagramma') %>%
    rowwise() %>% 
    mutate(id_col = paste(sort(c(P1, P2)), collapse = "_") ) %>%
    ungroup() %>% 
    distinct(id_col, .keep_all = TRUE)
  
  round(
    processed1 %>% 
      filter(abs(z_score) > 3) %>% 
      nrow()/nrow(processed1),
    2
  )
  
}, results = full_data_list$type_b)



############################
### 3. VISUALIZE D-STATS ###
############################

# expression(paste(italic("P"), ". sp. 'kazumbe' (Pk)", sep = ""))
# expression(paste(italic("P"), ". cf. ",  italic("polyodon"), " (Pp)", sep = ""))
# 
# expression(paste(italic("Petrochromis"), " sp. 'kazumbe' (Pk)", sep = ""))
# expression(paste(italic("Petrochromis"), " cf. ",  italic("polyodon"), " (Pp)", sep = ""))

### 3a: CREATION OF MATRIX PLOTS FROM ADMIXTOOLS RESULTS ###
#kazumbe_corrplot_title <- expression(paste(italic("P. kazumbe"), " (Pk)"))
#polyodon_corrplot_title <- expression(paste(italic("P. polyodon"), " (Pp)"))
kazumbe_corrplot_title <- expression(paste(italic("Petrochromis"), " sp. 'kazumbe' (Pk)", sep = ""))
polyodon_corrplot_title <- expression(paste(italic("Petrochromis"), " cf. ",  italic("polyodon"), " (Pp)", sep = ""))

scaleFUN <- function(x) sprintf("%.0f", x)
title_list <- list(kazumbe = kazumbe_corrplot_title,
                   polyodon = polyodon_corrplot_title)

full_plot_list <- list()
for (type in names(dstat_df_list)) {
  outgroup_plot_list <- list()
  
  for (out in names(dstat_df_list[[type]])) {
    taxa_plot_list <- list()
    
    for (taxa in names(dstat_df_list[[type]][[out]])) {
      matrix_plot <- ggplot(data = dstat_df_list[[type]][[out]][[taxa]], 
                                                          aes(x = p2_p3_location_factor_letter, y = p1_location_factor_letter)) + 
        geom_raster(fill = 'white') +
        #geom_point(aes(size = D_absolute), 
        geom_point(aes(size = z_score_absolute), 
                   color = dstat_df_list[[type]][[out]][[taxa]]$color_z, 
                   fill = dstat_df_list[[type]][[out]][[taxa]]$color_z) +
        scale_size_continuous(range = c(3, 23)) +
        theme(plot.background = element_rect(fill = "white"), 
              panel.background = element_rect(fill = "white"),
              #panel.border = element_rect(colour = "gray", fill=NA, size=1),
              axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5, size = 19),
              axis.text.y = element_text(size = 19),
              axis.ticks = element_blank(),
              axis.title = element_text(size = 21),
              axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 0, l = 0)),
              axis.title.x = element_text(margin = margin(t = 7, r = 0, b = 0, l = 0)),
              plot.title = element_blank(),
              legend.position = "none",
              plot.margin = margin(-2, 0, 5.5, 5.5)) +
        xlab("P2 population") + ylab("P1 population")
      
      bar_plot <- ggplot(data = dstat_df_list[[type]][[out]][[taxa]] %>% 
                               mutate(bar = 1,
                                      group = factor(case_when(z_score < -3 ~ "group1",
                                                               z_score < 0 & z_score >= -3 ~ "group2",
                                                               z_score == 0 ~ "group3",
                                                               z_score > 0 & z_score <= 3 ~ "group4",
                                                               z_score > 3 ~ "group5"), levels = paste0('group', 5:1))), 
                             aes(x = p2_p3_location_factor_letter, y = bar, fill = group) ) +
        geom_bar(position="fill", stat="identity", width = 0.9)  +
        scale_fill_manual("legend", values = c("group1" = "#88469e", 
                                               "group2" =  '#c2a5cf',
                                               "group3" =  "#a6a6a6",
                                               "group4" =  "#a6dba0", 
                                               "group5" =  "#19934b")) +
        theme(plot.background = element_rect(fill = "white"), 
              panel.background = element_rect(fill = "white"),
              axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              #axis.text.y = element_text(size = 13, color = 'white'),
              axis.title.y = element_text(color = 'white'),
              axis.text.y = element_blank(),
              axis.ticks = element_blank(),
              axis.title = element_text(size = 21),
              plot.title = element_text(size = 24, hjust = 0.5),
              legend.position = "none",
              plot.margin = margin(5.5, 0, 0, 27)) + #originally right was 5.5
        scale_y_continuous(labels = scaleFUN, expand = c(0, 0)) +
        ggtitle(title_list[[taxa]])
      
      taxa_plot_list[[taxa]] <- ggarrange(bar_plot, 
                                          matrix_plot, 
                                          ncol = 1, heights = c(0.4, 1.7))
    }
    outgroup_plot_list[[out]] <- taxa_plot_list
  }
  full_plot_list[[type]] <- outgroup_plot_list
}



########################
### CREATING LEGENDS ###
########################

pos_neg_legend_vertical <- ggplot(data = data.frame(x = c(1, 1),
                                                    x_end = c(1, 1),
                                                    y = c(-0.25, 10.25),
                                                    y_end = c(10.25, -0.25)) ) +
  geom_segment(aes(x = x, y = y, xend = x_end, yend = y_end),
               arrow = arrow(length = unit(0.325, "cm"), type = "closed"), size = 1.75) +
  geom_rect(aes(xmin = 0.92, xmax = 1.08, ymin = 4.89, ymax = 5.11), color = "black", fill = "black") +
  geom_point(data = data.frame(x = 1.7,
                               y = c(0, 2.5, 5, 7.5, 10)), 
             aes(x = x, y = y), 
             size = c(13, 7, 3.5, 7, 13),
             color = c("#88469e", '#c2a5cf', "#a6a6a6", "#a6dba0", "#19934b")) +
  geom_text(data = data.frame(x = 0.5,
                              y = c(0, 2.5, 5, 7.5, 10)), 
            aes(x = x, y = y),
            label = c('-', NA, '0', NA, '+'), size = c(9, 7, 7, 7, 9), hjust = 0.5, vjust = 0.5) +
  theme_void() +
  theme(plot.margin = margin(5, 0, 5, 3.5)) +
  xlim(0.3, 2.1) +
  coord_cartesian(clip = 'off')


z_legend <- ggplot(data.frame(x = c(0, 0, 1, 1),
                              y = c(0, 1, 0, 1))) +
  geom_point(aes(x = x, y = y), size = c(7, 13, 7, 13), color = c('#c2a5cf', "#88469e", "#a6dba0", "#19934b")) +
  xlim(-0.45, 3.55) + ylim(-0.5, 2) +
  geom_text(data = data.frame(x = 2.5, y = c(0, 1)), 
            aes(x = x, y = y), 
            label = c('|z| < 3', '|z| > 3'), size = 8, vjust = 0.5, hjust = 0.5) +
  theme_void()

empty_plot <- ggplot() + theme_void()


#kazumbe_corrplot_title <- expression(paste(italic("P. kazumbe"), " (Pk)"))
#polyodon_corrplot_title <- expression(paste(italic("P. polyodon"), " (Pp)"))
kazumbe_corrplot_title <- expression(paste(italic("P"), ". sp. 'kazumbe' (Pk)", sep = ""))
polyodon_corrplot_title <- expression(paste(italic("P"), ". cf. ",  italic("polyodon"), " (Pp)", sep = ""))



title_list <- list(kazumbe = kazumbe_corrplot_title,
                   polyodon = polyodon_corrplot_title)

#to add spaces, use ~~
tip_label_list_typea <- list(kazumbe = c('P1: Pk[pop1]', 'P2: Pk[pop2]', 'P3: Pp[full]', 'Outgroup'),
                             polyodon = c('P1: Pp[pop1]', 'P2: Pp[pop2]', 'P3: Pk[full]', 'Outgroup'))

tip_label_list_typeb <- list(kazumbe = c('P1: Pk[pop1]', 'P2: Pk[pop2]', 'P3: Pp[pop2]', 'Outgroup'),
                             polyodon = c('P1: Pp[pop1]', 'P2: Pp[pop2]', 'P3: Pk[pop2]', 'Outgroup'))

tree_plot_legend_list_typea <- list()
tree_plot_legend_list_typeb <- list()

for (i in c('polyodon', 'kazumbe')) {
  
  tree_plot_legend_list_typea[[i]] <- ggplot(data.frame(x_start = c(15,   10,  5, 10,  20, 30),
                                                  y_start = c(0,   10,  20,  30,   30,  30),
                                                  x_end =   c(10, 5, 0,    5, 10, 15),
                                                  y_end =   c(10, 20,  30,   20,  10,  0))) +
    geom_segment(aes(x = x_start, y = y_start, xend = x_end, yend = y_end), size = 2.5, lineend = 'round') +
    geom_text(data = data.frame(x = c(0, 10, 20, 30),
                                y = 33),
              aes(x = x, y = y),
              label = tip_label_list_typea[[i]],
              size = 5, parse = TRUE) + #size originally 5
    geom_segment(data = data.frame(x_start = c(25/3, 55/3, 10/3, 50/3),
                                   y_start = c(80/3, 80/3, 70/3, 70/3),
                                   x_end =   c(55/3, 25/3, 50/3, 10/3),
                                   y_end =   c(80/3, 80/3, 70/3, 70/3)),
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 arrow = arrow(length = unit(0.325, "cm"), type = "closed"), 
                 size = 1.75, color = c("#88469e", "#88469e", "#19934b", "#19934b")) +
    xlim(-5, 33) +
    theme_void() +
    ggtitle(title_list[[i]]) +
    theme(plot.title = element_text(size = 22, hjust = 0.5, margin = margin(0, 0, 2, 0))) +
    coord_cartesian(clip = 'off')
  
  tree_plot_legend_list_typeb[[i]] <- ggplot(data.frame(x_start = c(15,   10,  5, 10,  20, 30),
                                                        y_start = c(0,   10,  20,  30,   30,  30),
                                                        x_end =   c(10, 5, 0,    5, 10, 15),
                                                        y_end =   c(10, 20,  30,   20,  10,  0))) +
    geom_segment(aes(x = x_start, y = y_start, xend = x_end, yend = y_end), size = 2.5, lineend = 'round') +
    geom_text(data = data.frame(x = c(0, 10, 20, 30),
                                y = 33),
              aes(x = x, y = y),
              label = tip_label_list_typeb[[i]],
              size = 5, parse = TRUE) + #size originally 5.5
    geom_segment(data = data.frame(x_start = c(25/3, 55/3, 10/3, 50/3),
                                   y_start = c(80/3, 80/3, 70/3, 70/3),
                                   x_end =   c(55/3, 25/3, 50/3, 10/3),
                                   y_end =   c(80/3, 80/3, 70/3, 70/3)),
                 aes(x = x_start, y = y_start, xend = x_end, yend = y_end),
                 arrow = arrow(length = unit(0.325, "cm"), type = "closed"), 
                 size = 1.75, color = c("#88469e", "#88469e", "#19934b", "#19934b")) +
    xlim(-5, 33) +
    theme_void() +
    ggtitle(title_list[[i]]) +
    theme(plot.title = element_text(size = 21, hjust = 0.5, margin = margin(0, 0, 2, 0))) +
    coord_cartesian(clip = 'off')

}


### creating legends ###
tree_multipanel_typea <- plot_grid(tree_plot_legend_list_typea$kazumbe, empty_plot, tree_plot_legend_list_typea$polyodon, 
                                   nrow = 3, rel_heights = c(0.45, 0.08, 0.45))

tree_multipanel_typeb <- plot_grid(tree_plot_legend_list_typeb$kazumbe, empty_plot, tree_plot_legend_list_typeb$polyodon, 
                                   nrow = 3, rel_heights = c(0.45, 0.08, 0.45))

# pos_neg_plus_tree_typea <- plot_grid(plot_grid(empty_plot, pos_neg_legend_vertical, empty_plot, 
#                                                rel_heights = c(0.2, 0.6, 0.15), ncol = 1), 
#                                      tree_multipanel_typea, ncol = 2, rel_widths = c(0.32, 0.73))


right_legend <- plot_grid(empty_plot,
          plot_grid(empty_plot, pos_neg_legend_vertical, empty_plot, rel_heights = c(0.2, 0.6, 0.05), ncol = 1) + theme(plot.margin = margin(0, 20, 0, 18)),
          plot_grid(empty_plot, z_legend, empty_plot, rel_widths = c(0.05, 0.65, 0.05), ncol = 3),
          empty_plot,
          ncol = 1, 
          rel_heights = c(0.05, 0.85, 0.2, 0.1)) + theme(plot.margin = margin(0, 0, 0, 5))


#type a with S. diagramma outgroup
plot_grid(plot_grid(empty_plot, tree_multipanel_typeb, empty_plot, nrow = 3, rel_heights = c(0.08, 0.8, 0.08)), 
          plot_grid(full_plot_list$type_a$diagramma$kazumbe, 
                    full_plot_list$type_a$diagramma$polyodon, 
                    nrow = 1, labels = c("(b)", "(c)"), label_size = 26, vjust = 1), right_legend, 
          ncol = 3, rel_widths = c(0.25, 0.8, 0.14), labels = c("(a)", "", ""), label_size = 26, vjust = 1) +
  theme(plot.margin = margin(4, 0, 0, 0))

ggsave2(here('figures', 'dstat_matrix_plot_typea_outgroup_diagramma_9_19_2021_REARRANGE.png'), 
        width = 19*1.1, height = 8*1.1, bg = 'white')


#type a with P. green outgroup
plot_grid(plot_grid(empty_plot, tree_multipanel_typeb, empty_plot, nrow = 3, rel_heights = c(0.08, 0.8, 0.08)), 
          plot_grid(full_plot_list$type_a$green$kazumbe, 
                    full_plot_list$type_a$green$polyodon, 
                    nrow = 1, labels = c("(b)", "(c)"), label_size = 26, vjust = 1), right_legend, 
          ncol = 3, rel_widths = c(0.25, 0.8, 0.14), labels = c("(a)", "", ""), label_size = 26, vjust = 1) +
  theme(plot.margin = margin(4, 0, 0, 0))

ggsave2(here('figures', 'dstat_matrix_plot_typea_outgroup_green_9_19_2021_REARRANGE.png'), 
        width = 19*1.1, height = 8*1.1, bg = 'white')


#type b with S. diagramma outgroup
plot_grid(plot_grid(empty_plot, tree_multipanel_typea, empty_plot, nrow = 3, rel_heights = c(0.08, 0.8, 0.08)), 
          plot_grid(full_plot_list$type_b$diagramma$kazumbe, 
                    full_plot_list$type_b$diagramma$polyodon, 
                    nrow = 1, labels = c("(b)", "(c)"), label_size = 26, vjust = 1), right_legend, 
          ncol = 3, rel_widths = c(0.25, 0.8, 0.14), labels = c("(a)", "", ""), label_size = 26, vjust = 1) +
  theme(plot.margin = margin(4, 0, 0, 0))

ggsave2(here('figures', 'dstat_matrix_plot_typeb_outgroup_diagramma_9_19_2021_REARRANGE.png'), 
        width = 19*1.1, height = 8*1.1, bg = 'white')


#type a with P. green outgroup
plot_grid(plot_grid(empty_plot, tree_multipanel_typea, empty_plot, nrow = 3, rel_heights = c(0.08, 0.8, 0.08)), 
          plot_grid(full_plot_list$type_b$green$kazumbe, 
                    full_plot_list$type_b$green$polyodon, 
                    nrow = 1, labels = c("(b)", "(c)"), label_size = 26, vjust = 1), right_legend, 
          ncol = 3, rel_widths = c(0.25, 0.8, 0.14), labels = c("(a)", "", ""), label_size = 26, vjust = 1) +
  theme(plot.margin = margin(4, 0, 0, 0))

ggsave2(here('figures', 'dstat_matrix_plot_typeb_outgroup_green_9_19_2021_REARRANGE.png'), 
        width = 19*1.1, height = 8*1.1, bg = 'white')



#####################################
### 4. Tables of D stats results  ###
#####################################

processed_table_info <- lapply(dstat_df_list, function(TYPE) {
  lapply(TYPE, function(OUTGROUP) {
    lapply(OUTGROUP, function(SPECIES) {
      SPECIES %>% 
        select(p1_location_factor_letter, p2_p3_location_factor_letter, D, z_score) %>%
        arrange(p1_location_factor_letter, p2_p3_location_factor_letter) %>% 
        rename("P1 location" = p1_location_factor_letter,
               "P2/P3 location" = p2_p3_location_factor_letter,
               "Z-score" = z_score)
      
    })
  })
})


dstat_table_df_list <- lapply(setNames(nm = c('kazumbe', 'polyodon')), function(SPECIES, processed_table_info) {
  
  species_table_df <- left_join(processed_table_info$type_b$diagramma[[SPECIES]],
                                processed_table_info$type_b$green[[SPECIES]], by = c("P1 location", "P2/P3 location")) %>% 
    left_join(., processed_table_info$type_a$diagramma[[SPECIES]], by = c("P1 location", "P2/P3 location")) %>% 
    left_join(., processed_table_info$type_a$green[[SPECIES]], by = c("P1 location", "P2/P3 location"))
  
  colnames(species_table_df) <- c('P1 location', 'P2/P3 location',
                                        'D', 'Z-score',
                                        'D', 'Z-score',
                                        'D', 'Z-score',
                                        'D', 'Z-score')
  return(species_table_df)
  
}, processed_table_info = processed_table_info)


table_caption <- list(kazumbe = "Results for the D statistics associated with the tests examining variation in interspecific allele sharing between populations of \textit{Petrochromis} sp. `kazumbe'. The type a tests have the following thing",
                      polyodon = "Results for the D statistics associated with the tests examining variation in interspecific allele sharing between populations of \textit{Petrochromis} sp. `kazumbe'. The type a tests have the following thing")


latex_table_list <- lapply(setNames(nm = c('kazumbe', 'polyodon')), function(SPECIES, table_list, caption_list) {
  
  species_number <- switch(SPECIES, kazumbe = '1', polyodon = '2')
  a_top_label <- paste0('topology ', species_number, 'a')
  b_top_label <- paste0('topology ', species_number, 'b')
  
  topology_label <- c(" ", " ", a_top_label = 4, b_top_label = 4)
  names(topology_label) <- c(" ", " ", a_top_label, b_top_label)
  
  kbl(table_list[[SPECIES]], 'latex', caption = caption_list[[SPECIES]], longtable = TRUE, booktabs = TRUE, align = "c") %>% 
    add_header_above(c(" ", " ", "S. diagramma" = 2, "P. green" = 2, "S. diagramma" = 2, "P. green" = 2)) %>% 
    #add_header_above(c(" ", " ", a_top_label = 4, b_top_label = 4)) %>%
    add_header_above(header = topology_label) %>%
    kable_styling(font_size = 8)
  
}, table_list = dstat_table_df_list, caption_list = table_caption)

latex_table_list_kazumbe <- kbl(dstat_table_df_list$kazumbe, 'latex', caption = "caption", longtable = TRUE, booktabs = TRUE, align = "c") %>% 
  add_header_above(c(" ", " ", "diagramma" = 2, "green" = 2, "S. diagramma" = 2, "P. green" = 2)) %>% 
  add_header_above(c(" ", " ", "topology 1a" = 4, "topology 1b" = 4)) %>%
  kable_styling(font_size = 8)

latex_table_list_polyodon <- kbl(dstat_table_df_list$polyodon, 'latex', caption = "caption", longtable = TRUE, booktabs = TRUE, align = "c") %>% 
  add_header_above(c(" ", " ", "diagramma" = 2, "green" = 2, "S. diagramma" = 2, "P. green" = 2)) %>% 
  add_header_above(c(" ", " ", "topology 2a" = 4, "topology 2b" = 4)) %>%
  kable_styling(font_size = 8)

cat(latex_table_list_kazumbe, file = here('tables', 'kazumbe_dstat_table.txt'), append = FALSE)
cat(latex_table_list_polyodon, file = here('tables', 'polyodon_dstat_table.txt'), append = FALSE)



###################################
### 5. MANTEL TESTS OF D-STATS  ###
###################################
#Note: all Mantel tests are based on the d-statistics calculated with AdmixTools

#how many permutations to use in mantel tests?
number_permutations <- 10000

### TEST 1: tests between polyodon and kazumbe matrices (within same setup type and outgroup) ###
#e.g. type a setup of polyodon with P. green outgroup and type a setup of kazumbe with P. green outgroup

test1_comparison_vec_dstat <- c()
test1_mantel_stat_dstat_vec <- c()
test1_mantel_signif_dstat_vec <- c()

test1_comparison_vec_zscore <- c()
test1_mantel_stat_zscore_vec <- c()
test1_mantel_signif_zscore_vec <- c()

for (type in names(dstat_matrix_list)) {
  for (out in names(dstat_matrix_list[[type]])) {
    
    test1_dstat <- mantel(dstat_matrix_list[[type]][[out]]$polyodon$d_stat_val, 
                          dstat_matrix_list[[type]][[out]]$kazumbe$d_stat_val, 
                          method = "spearman", permutations = number_permutations)
    
    test1_comparison_vec_dstat <- c(test1_comparison_vec_dstat, paste('test1', type, out, 'dstat_val', sep = "_"))
    test1_mantel_stat_dstat_vec <- c(test1_mantel_stat_dstat_vec, test1_dstat$statistic)
    test1_mantel_signif_dstat_vec <- c(test1_mantel_signif_dstat_vec, test1_dstat$signif)
    
    
    test1_zscore <- mantel(dstat_matrix_list[[type]][[out]]$polyodon$d_stat_zscore, 
                           dstat_matrix_list[[type]][[out]]$kazumbe$d_stat_zscore, 
                           method = "spearman", permutations = number_permutations)
    
    test1_comparison_vec_zscore <- c(test1_comparison_vec_zscore, paste('test1', type, out, 'zscore', sep = "_"))
    test1_mantel_stat_zscore_vec <- c(test1_mantel_stat_zscore_vec, test1_zscore$statistic)
    test1_mantel_signif_zscore_vec <- c(test1_mantel_signif_zscore_vec, test1_zscore$signif)
    
  }
}

test1_dstat_val_df <- data.frame(comparison = test1_comparison_vec_dstat,
                                 mantel_stat = test1_mantel_stat_dstat_vec,
                                 mantel_signif = test1_mantel_signif_dstat_vec)

test1_dstat_zscore_df <- data.frame(comparison = test1_comparison_vec_zscore,
                                    mantel_stat = test1_mantel_stat_zscore_vec,
                                    mantel_signif = test1_mantel_signif_zscore_vec)


### TEST 2: tests between same taxa but different outgroups (within same setup type) ###
#e.g. type a setup of polyodon with P. green outgroup and type a setup of polyodon with S. diagramma outgroup

test2_comparison_vec_dstat <- c()
test2_mantel_stat_dstat_vec <- c()
test2_mantel_signif_dstat_vec <- c()

test2_comparison_vec_zscore <- c()
test2_mantel_stat_zscore_vec <- c()
test2_mantel_signif_zscore_vec <- c()

for (type in names(dstat_matrix_list)) {
  for (taxa in names(dstat_matrix_list$type_a$diagramma)) {
    
    test2_dstat <- mantel(dstat_matrix_list[[type]]$diagramma[[taxa]]$d_stat_val, 
                          dstat_matrix_list[[type]]$green[[taxa]]$d_stat_val, 
                          method = "spearman", permutations = number_permutations)
    
    test2_comparison_vec_dstat <- c(test2_comparison_vec_dstat, paste('test2', type, taxa, 'dstat_val', sep = "_"))
    test2_mantel_stat_dstat_vec <- c(test2_mantel_stat_dstat_vec, test2_dstat$statistic)
    test2_mantel_signif_dstat_vec <- c(test2_mantel_signif_dstat_vec, test2_dstat$signif)
    
    test2_zscore <- mantel(dstat_matrix_list[[type]]$diagramma[[taxa]]$d_stat_zscore, 
                           dstat_matrix_list[[type]]$green[[taxa]]$d_stat_zscore, 
                           method = "spearman", permutations = number_permutations)
    
    test2_comparison_vec_zscore <- c(test2_comparison_vec_zscore, paste('test2', type, taxa, 'zscore', sep = "_"))
    test2_mantel_stat_zscore_vec <- c(test2_mantel_stat_zscore_vec, test2_zscore$statistic)
    test2_mantel_signif_zscore_vec <- c(test2_mantel_signif_zscore_vec, test2_zscore$signif)
  }
}

test2_dstat_val_df <- data.frame(comparison = test2_comparison_vec_dstat,
                                 mantel_stat = test2_mantel_stat_dstat_vec,
                                 mantel_signif = test2_mantel_signif_dstat_vec)

test2_dstat_zscore_df <- data.frame(comparison = test2_comparison_vec_zscore,
                                    mantel_stat = test2_mantel_stat_zscore_vec,
                                    mantel_signif = test2_mantel_signif_zscore_vec)


### TEST 3 tests between same taxa but different type (with same outgroup) ###
#e.g. type a setup of polyodon with P. green outgroup and type b setup of polyodon with P. green outgroup

test3_comparison_vec_dstat <- c()
test3_mantel_stat_dstat_vec <- c()
test3_mantel_signif_dstat_vec <- c()

test3_comparison_vec_zscore <- c()
test3_mantel_stat_zscore_vec <- c()
test3_mantel_signif_zscore_vec <- c()

for (out in names(dstat_matrix_list$type_a) ) {
  for (taxa in names(dstat_matrix_list$type_a$diagramma)) {
    
    test3_dstat <- mantel(dstat_matrix_list$type_a[[out]][[taxa]]$d_stat_val, 
                          dstat_matrix_list$type_b[[out]][[taxa]]$d_stat_val, 
                          method = "spearman", permutations = number_permutations)
    
    test3_comparison_vec_dstat <- c(test3_comparison_vec_dstat, paste('test3', out, taxa, 'dstat_val', sep = "_"))
    test3_mantel_stat_dstat_vec <- c(test3_mantel_stat_dstat_vec, test3_dstat$statistic)
    test3_mantel_signif_dstat_vec <- c(test3_mantel_signif_dstat_vec, test3_dstat$signif)
    
    test3_zscore <- mantel(dstat_matrix_list$type_a[[out]][[taxa]]$d_stat_zscore, 
                           dstat_matrix_list$type_b[[out]][[taxa]]$d_stat_zscore,
                           method = "spearman", permutations = number_permutations)
    
    test3_comparison_vec_zscore <- c(test3_comparison_vec_zscore, paste('test3', out, taxa, 'zscore', sep = "_"))
    test3_mantel_stat_zscore_vec <- c(test3_mantel_stat_zscore_vec, test3_zscore$statistic)
    test3_mantel_signif_zscore_vec <- c(test3_mantel_signif_zscore_vec, test3_zscore$signif)
  }
}

test3_dstat_val_df <- data.frame(comparison = test3_comparison_vec_dstat,
                                 mantel_stat = test3_mantel_stat_dstat_vec,
                                 mantel_signif = test3_mantel_signif_dstat_vec)

test3_dstat_zscore_df <- data.frame(comparison = test3_comparison_vec_zscore,
                                    mantel_stat = test3_mantel_stat_zscore_vec,
                                    mantel_signif = test3_mantel_signif_zscore_vec)


### TEST 2: tests between same taxa but different outgroups (within same setup type) ###
#e.g. type a setup of polyodon with P. green outgroup and type a setup of polyodon with S. diagramma outgroup


### TEST 3 tests between same taxa but different type (with same outgroup) ###
#e.g. type a setup of polyodon with P. green outgroup and type b setup of polyodon with P. green outgroup


#TYPE B in the output
#1a --> kazumbe with p3 as full sample polyodon
#2a --> polyodon with p3 as full sample kazumbe

#TYPE A in the output
#1b --> kazumbe with p3 as the polyodon samples from the P2 location
#2b --> polyodon with p3 as the kazumbe samples from the P2 location

test2_processed <- test2_dstat_val_df %>% 
  mutate(`Focal species` = str_extract(comparison, pattern = "kazumbe|polyodon"),
         mantel_stat = round(mantel_stat, 3),
         mantel_signif = round(mantel_signif, 5),
         Outgroup = 'diagramma and green') %>%
  mutate(`Topology comparison` = case_when(str_detect(comparison, "type_a") & `Focal species` == 'kazumbe' ~ "1b vs 1b (different outgroups)",
                                           str_detect(comparison, "type_b") & `Focal species` == 'kazumbe' ~ "1a vs 1a (different outgroups)",
                                           str_detect(comparison, "type_a") & `Focal species` == 'polyodon' ~ "2b vs 2b (different outgroups)",
                                           str_detect(comparison, "type_b") & `Focal species` == 'polyodon' ~ "2a vs 2a (different outgroups)")) %>% 
  rename(`Correlation metric` = mantel_stat,
         P = mantel_signif) %>% 
  select(`Focal species`, Outgroup, `Topology comparison`, `Correlation metric`, P) %>% 
  arrange(`Focal species`) %>% 
  mutate(test = 'test2')

test3_processed <- test3_dstat_val_df %>%
  mutate(`Focal species` = str_extract(comparison, pattern = "kazumbe|polyodon"),
         `Outgroup` = str_extract(comparison, pattern = "diagramma|green"),
         mantel_stat = round(mantel_stat, 3),
         mantel_signif = round(mantel_signif, 5)) %>%
  mutate(`Topology comparison` = case_when(`Focal species` == 'kazumbe' ~ "1a vs 1b",
                                           `Focal species` == 'polyodon' ~ "2a vs 2b")) %>% 
  rename(`Correlation metric` = mantel_stat,
         P = mantel_signif) %>% 
  select(`Focal species`, Outgroup, `Topology comparison`, `Correlation metric`, P) %>% 
  arrange(`Focal species`) %>% 
  mutate(test = 'test3')

mantel_test_table_final <- rbind(test2_processed, test3_processed) %>%
  select(!test) %>% 
  kbl(., 'latex', booktabs = TRUE, align = "c") %>%
  kable_styling(latex_options = c("hold_position", "scale_down")) %>% 
  pack_rows("Evaluating the effect of outgroup identity", 1, 4) %>% 
  pack_rows("Evaluating the effect of samples included in outgroup", 5, 8)

cat(mantel_test_table_final, file = here('tables', 'dstat_mantel_test_table.txt'), append = FALSE)



################################################################
### 6. MISCELLANEOUS (CODE THAT IS CURRENTLY NOT BEING USED) ###
################################################################

### ALTERNATIVE PLOT SET-UP ###
# ### type a legend ###
# tree_multipanel_typea <- plot_grid(tree_plot_legend_list_typea$kazumbe, empty_plot, tree_plot_legend_list_typea$polyodon, 
#                                    nrow = 3, rel_heights = c(0.45, 0.05, 0.45))
# 
# pos_neg_plus_tree_typea <- plot_grid(plot_grid(empty_plot, pos_neg_legend_vertical, empty_plot, 
#                                                rel_heights = c(0.2, 0.6, 0.15), ncol = 1), 
#                                      tree_multipanel_typea, ncol = 2, rel_widths = c(0.32, 0.73))
# 
# 
# full_legend_vertical_updated_typea  <- plot_grid(empty_plot,
#                                                  pos_neg_plus_tree_typea,
#                                                  plot_grid(empty_plot, z_legend, empty_plot, rel_widths = c(0.25, 0.6, 0.25), ncol = 3),
#                                                  empty_plot, rel_heights = c(0.05, 0.85, 0.2, 0.05), ncol = 1)
# 
# 
# ### type b legend ###
# tree_multipanel_typeb <- plot_grid(tree_plot_legend_list_typeb$kazumbe, empty_plot, tree_plot_legend_list_typeb$polyodon, 
#                                    nrow = 3, rel_heights = c(0.45, 0.05, 0.45))
# 
# pos_neg_plus_tree_typeb <- plot_grid(plot_grid(empty_plot, pos_neg_legend_vertical, empty_plot, 
#                                                rel_heights = c(0.2, 0.6, 0.15), ncol = 1), 
#                                      tree_multipanel_typeb, ncol = 2, rel_widths = c(0.32, 0.73))
# 
# 
# full_legend_vertical_updated_typeb  <- plot_grid(empty_plot,
#                                                  pos_neg_plus_tree_typeb,
#                                                  plot_grid(empty_plot, z_legend, empty_plot, rel_widths = c(0.25, 0.6, 0.25), ncol = 3),
#                                                  empty_plot, rel_heights = c(0.05, 0.85, 0.2, 0.05), ncol = 1)
# 
# 
# ### 3b: CREATE FINAL ADMIXTOOLS MULTIPANEL PLOTS AND EXPORTING FINAL PLOTS ###
# #type a with S. diagramma outgroup
# plot_grid(plot_grid(full_plot_list$type_a$diagramma$kazumbe, 
#                     full_plot_list$type_a$diagramma$polyodon, 
#                     nrow = 1, labels = c("(a)", "(b)"), label_size = 26, vjust = 1) +
#             theme(plot.margin = margin(4, 0, 0, 0)), full_legend_vertical_updated_typeb, 
#           ncol = 2, rel_widths = c(0.76, 0.24))
# 
# ggsave2(here('figures', 'dstat_matrix_plot_typea_outgroup_diagramma_9_19_2021.png'), 
#         width = 18*1.1, height = 8*1.1, bg = 'white')
# 
# #type a with P. green outgroup
# plot_grid(plot_grid(full_plot_list$type_a$green$kazumbe, 
#                     full_plot_list$type_a$green$polyodon, 
#                     nrow = 1, labels = c("(a)", "(b)"), label_size = 26, vjust = 1) +
#             theme(plot.margin = margin(4, 0, 0, 0)), 
#           full_legend_vertical_updated_typeb, 
#           ncol = 2, rel_widths = c(0.76, 0.24))
# 
# ggsave2(here('figures', 'dstat_matrix_plot_typea_outgroup_green_9_19_2021.png'), 
#         width = 18*1.1, height = 8*1.1, bg = 'white')
# 
# 
# #type b with S. diagramma outgroup
# plot_grid(plot_grid(full_plot_list$type_b$diagramma$kazumbe, 
#                     full_plot_list$type_b$diagramma$polyodon, 
#                     nrow = 1, labels = c("(a)", "(b)"), label_size = 26, vjust = 1) +
#             theme(plot.margin = margin(4, 0, 0, 0)), full_legend_vertical_updated_typea, 
#           ncol = 2, rel_widths = c(0.76, 0.24))
# 
# ggsave2(here('figures', 'dstat_matrix_plot_typeb_outgroup_diagramma_9_19_2021.png'), 
#         width = 18*1.1, height = 8*1.1, bg = 'white')
# 
# 
# #type a with P. green outgroup
# plot_grid(plot_grid(full_plot_list$type_b$green$kazumbe, 
#                     full_plot_list$type_b$green$polyodon, 
#                     nrow = 1, labels = c("(a)", "(b)"), label_size = 26, vjust = 1) +
#             theme(plot.margin = margin(4, 0, 0, 0)), full_legend_vertical_updated_typea, 
#           ncol = 2, rel_widths = c(0.76, 0.24))
# 
# ggsave2(here('figures', 'dstat_matrix_plot_typeb_outgroup_green_9_19_2021.png'), 
#         width = 18*1.1, height = 8*1.1, bg = 'white')

# full_legend_vertical_updated_typea  <- plot_grid(empty_plot,
#                                                  pos_neg_plus_tree_typea,
#                                                  plot_grid(empty_plot, z_legend, empty_plot, rel_widths = c(0.25, 0.6, 0.25), ncol = 3),
#                                                  empty_plot, rel_heights = c(0.05, 0.85, 0.2, 0.05), ncol = 1)


# ### type b legend ###
# tree_multipanel_typeb <- plot_grid(tree_plot_legend_list_typeb$kazumbe, empty_plot, tree_plot_legend_list_typeb$polyodon, 
#                                    nrow = 3, rel_heights = c(0.45, 0.05, 0.45))
# 
# pos_neg_plus_tree_typeb <- plot_grid(plot_grid(empty_plot, pos_neg_legend_vertical, empty_plot, 
#                                                rel_heights = c(0.2, 0.6, 0.15), ncol = 1), 
#                                      tree_multipanel_typeb, ncol = 2, rel_widths = c(0.32, 0.73))
# 
# 
# full_legend_vertical_updated_typeb  <- plot_grid(empty_plot,
#                                                  pos_neg_plus_tree_typeb,
#                                                  plot_grid(empty_plot, z_legend, empty_plot, rel_widths = c(0.25, 0.6, 0.25), ncol = 3),
#                                                  empty_plot, rel_heights = c(0.05, 0.85, 0.2, 0.05), ncol = 1)






### BARPLOTS FOR TYPE A ###

# barplot_type_a_kazumbe_outgroup_diagramma <- ggplot() + 
#   geom_bar(data = type_a_list[['kazumbe']] %>%
#              filter(outgroup == "Sdiagramma"), aes(x = p2_p3_location_factor_plotting, fill = summarized_inference), position = "fill") +
#   scale_fill_manual(values = c("#06d6a0", "#ef476f", "#ffd166")) +
#   theme(#legend.position = "none", 
#     plot.background = element_rect(fill = "white"), 
#     panel.background = element_rect(fill = "white"),
#     axis.text.y=element_text(size = 18, color = 'black'),
#     axis.text.x=element_text(angle = 90, margin=margin(0,0, 0,0), hjust= 1, size = 18, color = 'black', vjust = 0.5),
#     axis.ticks.x=element_blank(),
#     axis.title.x = element_text(size = 22, face = "bold"),
#     axis.title.y = element_text(size = 22, face = "bold"),
#     legend.text = element_text(size = 18),
#     legend.title = element_blank(),
#     plot.title = element_text(size = 25, hjust = 0.5),
#     legend.position = "none") +
#   xlab("Focal population location") +
#   ylab("Proportion of tests") +
#   ggtitle("((kazumbe, kazumbe), polyodon)")


# barplot_type_a_polyodon_outgroup_diagramma <- ggplot() + 
#   geom_bar(data = type_a_list[['polyodon']] %>%
#              filter(outgroup == "Sdiagramma"), aes(x = p2_p3_location_factor_plotting, fill = summarized_inference), position = "fill") +
#   scale_fill_manual(values = c("#06d6a0", "#ef476f", "#ffd166")) +
#   theme(#legend.position = "none", 
#     plot.background = element_rect(fill = "white"), 
#     panel.background = element_rect(fill = "white"),
#     axis.text.y=element_text(size = 18, color = 'black'),
#     axis.text.x=element_text(angle = 90, margin=margin(0,0, 0,0), hjust= 1, size = 18, color = 'black', vjust = 0.5),
#     axis.ticks.x=element_blank(),
#     axis.title.x = element_text(size = 22, face = "bold"),
#     axis.title.y = element_text(size = 22, face = "bold"),
#     legend.text = element_text(size = 18),
#     legend.title = element_blank(),
#     plot.title = element_text(size = 25, hjust = 0.5),
#     legend.position = "none") +
#   xlab("Focal population location") +
#   ylab("Proportion of tests") +
#   ggtitle("((polyodon, polyodon), kazumbe)") 

#barplots_type_a_outgroup_diagramma <- cowplot::plot_grid(cowplot::plot_grid(plot_type_a_kazumbe, plot_type_a_polyodon, nrow = 1, labels = c("(a)", "(b)"), label_size = 22), d_stat_barplot_legend, 
#                                                         ncol = 2, rel_widths = c(1, 0.2))


# barplot_type_b_kazumbe_outgroup_diagramma <- ggplot() + 
#   geom_bar(data = type_b_list[['kazumbe']] %>%
#              filter(outgroup == "Sdiagramma"), aes(x = p2_p3_location_factor_plotting, fill = summarized_inference), position = "fill") +
#   scale_fill_manual(values = c("#06d6a0", "#ef476f", "#ffd166")) +
#   theme(#legend.position = "none", 
#     plot.background = element_rect(fill = "white"), 
#     panel.background = element_rect(fill = "white"),
#     axis.text.y=element_text(size = 18, color = 'black'),
#     axis.text.x=element_text(angle = 90, margin=margin(0,0, 0,0), hjust= 1, size = 18, color = 'black', vjust = 0.5),
#     axis.ticks.x=element_blank(),
#     axis.title.x = element_text(size = 22, face = "bold"),
#     axis.title.y = element_text(size = 22, face = "bold"),
#     legend.text = element_text(size = 18),
#     legend.title = element_blank(),
#     plot.title = element_text(size = 25, hjust = 0.5),
#     legend.position = "none") +
#   xlab("Focal population location") +
#   ylab("Proportion of tests") +
#   ggtitle("((kazumbe, kazumbe), polyodon)")
# 
# 
# barplot_type_b_polyodon_outgroup_diagramma <- ggplot() + 
#   geom_bar(data = type_b_list[['polyodon']] %>%
#              filter(outgroup == "Sdiagramma"), aes(x = p2_p3_location_factor_plotting, fill = summarized_inference), position = "fill") +
#   scale_fill_manual(values = c("#06d6a0", "#ef476f", "#ffd166")) +
#   theme(#legend.position = "none", 
#     plot.background = element_rect(fill = "white"), 
#     panel.background = element_rect(fill = "white"),
#     axis.text.y=element_text(size = 18, color = 'black'),
#     axis.text.x=element_text(angle = 90, margin=margin(0,0, 0,0), hjust= 1, size = 18, color = 'black', vjust = 0.5),
#     axis.ticks.x=element_blank(),
#     axis.title.x = element_text(size = 22, face = "bold"),
#     axis.title.y = element_text(size = 22, face = "bold"),
#     legend.text = element_text(size = 18),
#     legend.title = element_blank(),
#     plot.title = element_text(size = 25, hjust = 0.5),
#     legend.position = "none") +
#   xlab("Focal population location") +
#   ylab("Proportion of tests") +
#   ggtitle("((polyodon, polyodon), kazumbe)") 



# barplots_type_b <- cowplot::plot_grid(cowplot::plot_grid(plot_type_b_kazumbe_outgroup_diagramma, 
#                                                          plot_type_b_polyodon_outgroup_diagramma, 
#                                                          nrow = 1, labels = c("(a)", "(b)"), label_size = 22),
#                                       d_stat_barplot_legend, ncol = 2, rel_widths = c(1, 0.2))

#barplots_type_b
#ggsave2("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/DStats/Results/results_6_19_2020/plots_6_19_2020/D_stat_diagramma_outgroup_barplots_typeb.png", width = 18, height = 10)



### (DSUITE) 2e. create a hierarchical list of matrices of D-stat values and D-stat z-scores (used in Mantel tests) ###

# dsuite_list <- list(diagramma = dsuite_out_diagramma_bin1000,
#                     green = dsuite_out_green_bin1000)
# 
# dsuite_processed_list <- list()
# 
# for (out in names(dsuite_list)) {
#   
#   #extract species identity for p1, p2, and p3
#   dsuite_list[[out]]$P1_sp <- gsub('_[A-Za-z]*', '', dsuite_list[[out]]$P1)
#   dsuite_list[[out]]$P2_sp <- gsub('_[A-Za-z]*', '', dsuite_list[[out]]$P2)
#   dsuite_list[[out]]$P3_sp <- gsub('_[A-Za-z]*', '', dsuite_list[[out]]$P3)
#   
#   #adjust p-values using Benjamini-Hochberg correction
#   dsuite_list[[out]]$p_adjust_BH <- p.adjust(dsuite_list[[out]]$p.value, method = 'BH')
#   dsuite_list[[out]]$sig_BH <- ifelse(dsuite_list[[out]]$p_adjust_BH < 0.05, 'sig', 'non_sig')
#   
#   #adjust p-values using Bonferroni correction
#   dsuite_list[[out]]$p_adjust_bonf <- p.adjust(dsuite_list[[out]]$p.value, method = 'bonferroni')
#   dsuite_list[[out]]$sig_bonf <- ifelse(dsuite_list[[out]]$p_adjust_bonf < 0.05, 'sig', 'non_sig')
#   
#   
#   #extract comparisons that contain 2 species (remove trios of only one species)
#   step1 <- dsuite_list[[out]][apply(dsuite_list[[out]][, c('P1_sp', 'P2_sp', 'P3_sp')], 1, function(x) {length(unique(x)) > 1} ),]
#   
#   
#   taxa_vec <- unique(c(step1$P1, step1$P2, step1$P3))
#   
#   d_stat_processed_list <- list()
#   for (i in taxa_vec) {
#     subset <- step1[apply(step1[, c('P1', 'P2', 'P3')], 1, function(x) {i %in% x} ) ,]
#     subset$focal_taxa <- i
#     subset$focal_taxa_location <- str_remove(i, '[A-Za-z]*_')
#     subset$focal_species <- str_extract(i, '[A-Za-z]*')
#     
#     d_stat_processed_list[[i]] <- subset
#   }
#   dsuite_processed_list[[out]] <- do.call(rbind, d_stat_processed_list)
#   
#   dsuite_processed_list[[out]]$focal_taxa_location_plotting[dsuite_processed_list[[out]]$focal_taxa_location == 'Gombe_South'] <- "Gombe S"
#   dsuite_processed_list[[out]]$focal_taxa_location_plotting[dsuite_processed_list[[out]]$focal_taxa_location == 'Katongwe_N'] <- "Katongwe N"
#   dsuite_processed_list[[out]]$focal_taxa_location_plotting[dsuite_processed_list[[out]]$focal_taxa_location == 'Katongwe_S'] <- "Katongwe S"
#   dsuite_processed_list[[out]]$focal_taxa_location_plotting[dsuite_processed_list[[out]]$focal_taxa_location == 'Ska'] <- "S Kagango"
#   dsuite_processed_list[[out]]$focal_taxa_location_plotting[dsuite_processed_list[[out]]$focal_taxa_location == 'Kalala'] <- "Kalalangabo"
#   dsuite_processed_list[[out]]$focal_taxa_location_plotting[dsuite_processed_list[[out]]$focal_taxa_location == 'Nondwa'] <- "Nondwa"
#   dsuite_processed_list[[out]]$focal_taxa_location_plotting[dsuite_processed_list[[out]]$focal_taxa_location == 'Hilltop'] <- "Hilltop"
#   dsuite_processed_list[[out]]$focal_taxa_location_plotting[dsuite_processed_list[[out]]$focal_taxa_location == 'Bangwe'] <- "Bangwe"
#   dsuite_processed_list[[out]]$focal_taxa_location_plotting[dsuite_processed_list[[out]]$focal_taxa_location == 'Jakob'] <- "Jakobsen's S"
#   dsuite_processed_list[[out]]$focal_taxa_location_plotting[dsuite_processed_list[[out]]$focal_taxa_location == 'Ulombola'] <- "Ulombola"
#   dsuite_processed_list[[out]]$focal_taxa_location_plotting[dsuite_processed_list[[out]]$focal_taxa_location == 'Harembe'] <- "Harembe"
#   dsuite_processed_list[[out]]$focal_taxa_location_plotting <- factor(dsuite_processed_list[[out]]$focal_taxa_location_plotting, 
#                                                                       levels= c("Gombe S", "Katongwe N", "Katongwe S", "S Kagango", "Kalalangabo", "Nondwa", "Hilltop", "Bangwe", "Jakobsen's S", "Ulombola"))
# }



### 3c: CREATE AND EXPORT DSUITE BARPLOTS ###

# dsuite_plot_list <- list()
# sig_col <- 'sig_BH'
# 
# 
# for (out in names(dsuite_processed_list)) {
#   
#   barplot_dsuite_kazumbe <- ggplot() + 
#     geom_bar(data = dsuite_processed_list[[out]] %>%
#                filter(focal_species == 'Pkazumbe', !(focal_taxa_location %in% c('Bangwe', 'Ulombola')) ), 
#              aes_string(x = 'focal_taxa_location_plotting', fill = sig_col), position = 'fill') +
#     scale_fill_manual(values = c("#06d6a0", "#ef476f")) +
#     theme(#legend.position = "none", 
#       plot.background = element_rect(fill = "white"), 
#       panel.background = element_rect(fill = "white"),
#       axis.text.y=element_text(size = 18, color = 'black'),
#       axis.text.x=element_text(angle = 90, margin=margin(0,0, 0,0), hjust= 1, size = 18, color = 'black', vjust = 0.5),
#       axis.ticks.x=element_blank(),
#       axis.title.x = element_text(size = 22, face = "bold"),
#       axis.title.y = element_text(size = 22, face = "bold"),
#       legend.text = element_text(size = 18),
#       legend.title = element_blank(),
#       plot.title = element_text(size = 25, hjust = 0.5),
#       legend.position = "none") +
#     xlab("Focal population location") +
#     ylab("Proportion of tests") +
#     ggtitle(expression(italic("P. kazumbe")))
#   
#   barplot_dsuite_polyodon <- ggplot() + 
#     geom_bar(data = dsuite_processed_list[[out]] %>%
#                filter(focal_species == 'Ppolyodon'), 
#              aes_string(x = 'focal_taxa_location_plotting', fill = sig_col), position = 'fill') +
#     scale_fill_manual(values = c("#06d6a0", "#ef476f")) +
#     theme(#legend.position = "none", 
#       plot.background = element_rect(fill = "white"), 
#       panel.background = element_rect(fill = "white"),
#       axis.text.y=element_text(size = 18, color = 'black'),
#       axis.text.x=element_text(angle = 90, margin=margin(0,0, 0,0), hjust= 1, size = 18, color = 'black', vjust = 0.5),
#       axis.ticks.x=element_blank(),
#       axis.title.x = element_text(size = 22, face = "bold"),
#       axis.title.y = element_text(size = 22, face = "bold"),
#       legend.text = element_text(size = 18),
#       legend.title = element_blank(),
#       plot.title = element_text(size = 25, hjust = 0.5),
#       legend.position = "none") +
#     xlab("Focal population location") +
#     ylab("Proportion of tests") +
#     ggtitle(expression(italic("P. polyodon")))
#   
#   
#   legend_dsuite_barplot <- get_legend(ggplot() + 
#                                         geom_bar(data = dsuite_processed_list[[out]] %>%
#                                                    filter(focal_species == 'Ppolyodon'), 
#                                                  aes_string(x = 'focal_taxa_location_plotting', fill = sig_col), position = 'fill') +
#                                         scale_fill_manual(values = c("#06d6a0", "#ef476f")) +
#                                         theme(#legend.position = "none", 
#                                           plot.background = element_rect(fill = "white"), 
#                                           panel.background = element_rect(fill = "white"),
#                                           axis.text.y=element_text(size = 18, color = 'black'),
#                                           axis.text.x=element_text(angle = 90, margin=margin(0,0, 0,0), hjust= 1, size = 18, color = 'black', vjust = 0.5),
#                                           axis.ticks.x=element_blank(),
#                                           axis.title.x = element_text(size = 22, face = "bold"),
#                                           axis.title.y = element_text(size = 22, face = "bold"),
#                                           legend.text = element_text(size = 18),
#                                           legend.title = element_blank(),
#                                           plot.title = element_text(size = 25, hjust = 0.5)) +
#                                         xlab("Focal population location") +
#                                         ylab("Proportion of tests"))
#   
#   dsuite_plot_list[[out]] <- cowplot::plot_grid(cowplot::plot_grid(barplot_dsuite_kazumbe, barplot_dsuite_polyodon, nrow = 1, labels = c("(a)", "(b)"), label_size = 22),
#                                                 legend_dsuite_barplot, ncol = 2, rel_widths = c(1, 0.2))
# }
# 
# #export the dsuite multipanel barplots
# dsuite_plot_list$diagramma
# ggsave2("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/DStats/Results/results_9_7_2020/plots_9_7_2020/DSuite_out_diagramma_outgroup_barplot_9_7_2020.png", width = 18, height = 8)
# 
# dsuite_plot_list$green
# ggsave2("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/DStats/Results/results_9_7_2020/plots_9_7_2020/DSuite_out_green_outgroup_barplot_9_7_2020.png", width = 18, height = 8)





#             
# kazumbe_latex_table <- kbl(kazumbe_dstat_table_df, 'latex', caption = kazumbe_caption, longtable = TRUE, booktabs = TRUE, align = "c") %>% 
#   add_header_above(c(" ", " ", "S. diagramma" = 2, "P. green" = 2, "S. diagramma" = 2, "P. green" = 2)) %>% 
#   add_header_above(c(" ", " ", "type a" = 4, "type b" = 4)) %>%
#   kable_styling(font_size = 8)
# #kable_styling(latex_options = c("repeat_header"), font_size = 8)
# 
# 
# 
# kazumbe_dstat_table_df <- left_join(kazumbe_dstat_list$type_a$diagramma,
#                                     kazumbe_dstat_list$type_a$green, by = c("P1 location", "P2/P3 location")) %>% 
#   left_join(., kazumbe_dstat_list$type_b$diagramma, by = c("P1 location", "P2/P3 location")) %>% 
#   left_join(., kazumbe_dstat_list$type_b$green, by = c("P1 location", "P2/P3 location"))
# 
# colnames(kazumbe_dstat_table_df) <- c('P1 location', 'P2/P3 location',
#                                       'D', 'Z-score',
#                                       'D', 'Z-score',
#                                       'D', 'Z-score',
#                                       'D', 'Z-score')
# 
# kazumbe_dstat_list <- lapply(dstat_df_list, function(x) {
#   lapply(x, function(y) { 
#     y[['kazumbe']] %>% 
#       select(p1_location_factor_letter, p2_p3_location_factor_letter, D, z_score) %>%
#       arrange(p1_location_factor_letter, p2_p3_location_factor_letter) %>% 
#       rename("P1 location" = p1_location_factor_letter,
#              "P2/P3 location" = p2_p3_location_factor_letter,
#              "Z-score" = z_score)
#   })
# })
# 
# 
# kazumbe_dstat_table_df <- left_join(kazumbe_dstat_list$type_a$diagramma,
#           kazumbe_dstat_list$type_a$green, by = c("P1 location", "P2/P3 location")) %>% 
#   left_join(., kazumbe_dstat_list$type_b$diagramma, by = c("P1 location", "P2/P3 location")) %>% 
#   left_join(., kazumbe_dstat_list$type_b$green, by = c("P1 location", "P2/P3 location"))
# 
# colnames(kazumbe_dstat_table_df) <- c('P1 location', 'P2/P3 location',
#                              'D', 'Z-score',
#                              'D', 'Z-score',
#                              'D', 'Z-score',
#                              'D', 'Z-score')
# 
# kazumbe_caption <- "Results for the D statistics associated with the tests examining variation in interspecific allele sharing between populations of \textit{Petrochromis} sp. `kazumbe'. The type a tests have the followign thing"
# kazumbe_latex_table <- kbl(kazumbe_dstat_table_df, 'latex', caption = kazumbe_caption, longtable = TRUE, booktabs = TRUE, align = "c") %>% 
#   add_header_above(c(" ", " ", "S. diagramma" = 2, "P. green" = 2, "S. diagramma" = 2, "P. green" = 2)) %>% 
#   add_header_above(c(" ", " ", "type a" = 4, "type b" = 4)) %>%
#   kable_styling(font_size = 8)
#   #kable_styling(latex_options = c("repeat_header"), font_size = 8)
# 
# cat(kazumbe_latex_table, file = paste0(table_dir, 'kazumbe_dstat_table.txt'), append = FALSE)
# 
# 
# 
# 
# polyodon_dstat_list <- lapply(dstat_df_list, function(x) {
#   lapply(x, function(y) { 
#     y[['polyodon']] %>% 
#       select(p1_location_factor_letter, p2_p3_location_factor_letter, D, z_score) %>%
#       arrange(p1_location_factor_letter, p2_p3_location_factor_letter) %>% 
#       rename("P1 location" = p1_location_factor_letter,
#              "P2/P3 location" = p2_p3_location_factor_letter,
#              "Z-score" = z_score)
#   })
# })
# 
# 
# polyodon_dstat_table_df <- left_join(polyodon_dstat_list$type_a$Sdiagramma,
#                                      polyodon_dstat_list$type_a$Pgreen, by = c("P1 location", "P2/P3 location")) %>% 
#   left_join(.,polyodon_dstat_list$type_b$Sdiagramma, by = c("P1 location", "P2/P3 location")) %>% 
#   left_join(., polyodon_dstat_list$type_b$Pgreen, by = c("P1 location", "P2/P3 location"))
# 
# colnames(polyodon_dstat_table_df) <- c('P1 location', 'P2/P3 location',
#                                       'D', 'Z-score',
#                                       'D', 'Z-score',
#                                       'D', 'Z-score',
#                                       'D', 'Z-score')
# 
# 
# polyodon_caption <- "Results for the D statistics associated with the tests examining variation in interspecific allele sharing between populations of \textit{Petrochromis} sp. `kazumbe'. The type a tests have the followign thing"
# polyodon_latex_table <- kbl(polyodon_dstat_table_df, 'latex', caption = kazumbe_caption, longtable = TRUE, booktabs = TRUE, align = "c") %>% 
#   add_header_above(c(" ", " ", "S. diagramma" = 2, "P. green" = 2, "S. diagramma" = 2, "P. green" = 2)) %>% 
#   add_header_above(c(" ", " ", "type a" = 4, "type b" = 4)) %>% 
#   kable_styling(latex_options = c("repeat_header"), font_size = 8)
# 
# 
# # polyodon_latex_table <- kbl(polyodon_dstat_table_df, 'latex', booktabs = TRUE, align = "c") %>%
# #   kable_styling(latex_options = "scale_down") %>% 
# #   add_header_above(c(" ", " ", "S. diagramma" = 2, "p. green" = 2, "S. diagramma" = 2, "P. green" = 2)) %>% 
# #   add_header_above(c(" ", " ", "type a" = 4, "type b" = 4))
# 
# cat(polyodon_latex_table, file = paste0(table_dir, 'kazumbe_dstat_table.txt'), append = FALSE)
# 


# ### 3b: CREATE FINAL ADMIXTOOLS MULTIPANEL PLOTS AND EXPORTING FINAL PLOTS ###
# #type a with S. diagramma outgroup
# (matrix_plot_typea_outgroup_diagramma <- cowplot::plot_grid(full_plot_list$type_a$Sdiagramma$kazumbe, 
#                                                             full_plot_list$type_a$Sdiagramma$polyodon, 
#                                                             nrow = 1, labels = c("(a)", "(b)"), label_size = 22))
# 
# ggsave2(paste0(dstat_plot_output_path, 'dstat_matrix_plot_typea_outgroup_diagramma_10_18_2020.png'), width = 18, height = 8)
# 
# 
# #type a with P. green outgroup
# (matrix_plot_typea_outgroup_green <- cowplot::plot_grid(full_plot_list$type_a$Pgreen$kazumbe, 
#                                                         full_plot_list$type_a$Pgreen$polyodon, 
#                                                         nrow = 1, labels = c("(a)", "(b)"), label_size = 22))
# 
# ggsave2(paste0(dstat_plot_output_path, 'dstat_matrix_plot_typea_outgroup_green_10_18_2020.png'), width = 18, height = 8)
# 
# 
# #type b with S. diagramma outgroup
# (matrix_plot_typeb_outgroup_diagramma <- cowplot::plot_grid(full_plot_list$type_b$Sdiagramma$kazumbe, 
#                                                             full_plot_list$type_b$Sdiagramma$polyodon, 
#                                                             nrow = 1, labels = c("(a)", "(b)"), label_size = 22))
# 
# ggsave2(paste0(dstat_plot_output_path, 'dstat_matrix_plot_typeb_outgroup_diagramma_10_18_2020.png'), width = 18, height = 8)
# 
# 
# #type b with P. green outgroup
# (matrix_plot_typeb_outgroup_green <- cowplot::plot_grid(full_plot_list$type_b$Pgreen$kazumbe, 
#                                                         full_plot_list$type_b$Pgreen$polyodon, 
#                                                         nrow = 1, labels = c("(a)", "(b)"), label_size = 22))
# 
# ggsave2(paste0(dstat_plot_output_path, 'dstat_matrix_plot_typeb_outgroup_green_10_18_2020.png'), width = 18, height = 8)




# scaleFUN <- function(x) sprintf("%.0f", x)
# title_list <- list(kazumbe = kazumbe_corrplot_title,
#                    polyodon = polyodon_corrplot_title)
# 
# bar_plot_list <- list()
# for (type in names(dstat_df_list)) {
#   outgroup_plot_list <- list()
#   
#   for (out in names(dstat_df_list[[type]])) {
#     taxa_plot_list <- list()
#     for (taxa in names(dstat_df_list[[type]][[out]])) {
#       
#       bar_plot <- ggplot(data = dstat_df_list[[type]][[out]][[taxa]] %>% 
#                            mutate(bar = 1,
#                                   group = factor(case_when(z_score < -3 ~ "group1",
#                                                            z_score < 0 & z_score >= -3 ~ "group2",
#                                                            z_score == 0 ~ "group3",
#                                                            z_score > 0 & z_score <= 3 ~ "group4",
#                                                            z_score > 3 ~ "group5"), levels = paste0('group', 5:1))), 
#                          aes(x = p2_p3_location_factor_letter, y = bar, fill = group) ) +
#         geom_bar(position="fill", stat = "identity", width = 0.9)  +
#         scale_fill_manual("legend", values = c("group1" = "#88469e", 
#                                                "group2" =  '#c2a5cf',
#                                                "group3" =  "#a6a6a6",
#                                                "group4" =  "#a6dba0", 
#                                                "group5" =  "#19934b")) +
#         theme(plot.background = element_rect(fill = "white"), 
#               panel.background = element_rect(fill = "white"),
#               axis.text.x = element_blank(),
#               axis.title.x = element_blank(),
#               #axis.text.y = element_text(size = 13, color = 'white'),
#               axis.title.y = element_text(color = 'white'),
#               axis.text.y = element_blank(),
#               axis.ticks = element_blank(),
#               axis.title = element_text(size = 21),
#               plot.title = element_text(size = 24, hjust = 0.5),
#               legend.position = "none",
#               plot.margin = margin(5.5, 0, 0, 27)) + #originally right was 5.5
#         scale_y_continuous(labels = scaleFUN, expand = c(0, 0)) +
#         ggtitle(title_list[[taxa]])
#       
#       taxa_plot_list[[taxa]] <- bar_plot
#     }
#     outgroup_plot_list[[out]] <- taxa_plot_list
#   }
#   bar_plot_list[[type]] <- outgroup_plot_list
# }
# 
# 
# bar_plot_list$type_b$diagramma$kazumbe
# ggsave2('/Users/alexlewanski/Documents/University_of_Wyoming/defense_presentation/figures_images/dstat_kazumbe_barplot.png', 
#         width = 10, height = 8, bg = 'white')
# 
# bar_plot_list$type_b$diagramma$polyodon
# ggsave2('/Users/alexlewanski/Documents/University_of_Wyoming/defense_presentation/figures_images/dstat_polyodon_barplot.png', 
#         width = 10, height = 8, bg = 'white')


# ### TEST 4: COMPARING DMIN CALCULATIONS WITH THE DIFFERENT OUTGROUPS ###
# 
# dsuite_processed_list_diagramma_reduce <- dsuite_processed_list$diagramma[, c('P1','P2', 'P3', 'Dstatistic')]
# colnames(dsuite_processed_list_diagramma_reduce)[colnames(dsuite_processed_list_diagramma_reduce) == 'Dstatistic'] <- 'Dstatistic_diagramma' 
# dsuite_processed_list_diagramma_reduce$P1_P2_P2 <- paste(dsuite_processed_list_diagramma_reduce$P1,
#                                                          dsuite_processed_list_diagramma_reduce$P2,
#                                                          dsuite_processed_list_diagramma_reduce$P3, sep = "_")
# 
# dsuite_processed_list_diagramma_reduce_noduplicates <- dsuite_processed_list_diagramma_reduce[!duplicated(dsuite_processed_list_diagramma_reduce$P1_P2_P2),]
# 
# 
# dsuite_processed_list_green_reduce <- dsuite_processed_list$green[, c('P1','P2', 'P3', 'Dstatistic')]
# colnames(dsuite_processed_list_green_reduce)[colnames(dsuite_processed_list_green_reduce) == 'Dstatistic'] <- 'Dstatistic_green' 
# dsuite_processed_list_green_reduce$P1_P2_P2 <- paste(dsuite_processed_list_green_reduce$P1,
#                                                          dsuite_processed_list_green_reduce$P2,
#                                                          dsuite_processed_list_green_reduce$P3, sep = "_")
# 
# dsuite_processed_list_green_reduce_noduplicates <- dsuite_processed_list_green_reduce[!duplicated(dsuite_processed_list_green_reduce$P1_P2_P2),]
# 
# 
# dsuite_outgroup_merge <- merge(dsuite_processed_list_diagramma_reduce_noduplicates, dsuite_processed_list_green_reduce_noduplicates,
#                                by= 'P1_P2_P2')
# 
# setdiff(dsuite_processed_list_diagramma_reduce_noduplicates$P1_P2_P2, dsuite_processed_list_green_reduce_noduplicates$P1_P2_P2)
# setdiff(dsuite_processed_list_green_reduce_noduplicates$P1_P2_P2, dsuite_processed_list_diagramma_reduce_noduplicates$P1_P2_P2)

