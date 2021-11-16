#######################################################################
#######################################################################
### FUNCTIONS FOR VISUALIZING ENTROPY q VALUES (STRUCTURE BARPLOTS) ###
#######################################################################
#######################################################################


#DESCRIPTION OF HIERARCHICHAL_ORDERING_FUNC:
#hierarchical_ordering_func takes two or three variables (the last variable is continuous and the others are categorical), 
#and then hierarchically sorts the factor levels of a another variable based on those variables. For example, when given
#three variables to sort by (e.g. Var1, Var2, Var3) and a target variable to sort (e.g Targ1), the function first sorts the
#levels of Targ1 based on Var1, and then sorts by Var2 within each level of Var1, and then sorts by either decreasing or 
#increasing order of Var3 within each Var1 level X Var2 level combination. 

#HIERARCHICHAL_ORDERING_FUNC UTILITY:
#Since ggplot orders objects based on the levels of a factor, hierarchical_ordering_func can be used to reorder a variable
#so that is visualized in a particular order. For example, for the q plots of entropy, I use hierarchical_ordering_func to
#order the individuals in a way so that they organized by location, organized by primary ancestry group within each 
#location, and order by proportion of primary ancestry group within each primary ancestry group.

#currently only handles two and three level hierarchical ordering where the lowest level is a continuous variable
#will need to make this into some sort of recursive function to generalize to any numbers of variables
hierarchical_ordering_func <- function(data, variable_order, level_order, variable_to_sort, decreasing = TRUE) {
  
  order_vec <- c()
  
  if (length(variable_order) == 3) {
    
    ifelse(!is.na(level_order[[1]]), 
           v1_ord <- level_order[[1]], 
           v1_ord <- unique(data[,variable_order[1]]) )
    
    for (i in v1_ord) {
      data_subset1 <- data[which(data[,variable_order[1]] == i),]
      
      ifelse(!is.na(level_order[[2]]), 
             v2_ord <- level_order[[2]], 
             v2_ord <- unique(data_subset1[,variable_order[2]]) )
      
      for (j in v2_ord) {
        data_subset2 <- data_subset1[which(data_subset1[,variable_order[2]] == j),]
        
        order_vec <- c(order_vec, unique(as.character(data_subset2[order(data_subset2[, variable_order[3]], decreasing = decreasing), variable_to_sort])) )
      }
    }
  }
  
  if (length(variable_order) == 2) {
    
    ifelse(!is.na(level_order[[1]]), 
           v1_ord <- level_order[[1]], 
           v1_ord <- unique(data[,variable_order[1]]) )
    
    for (i in v1_ord) {
      data_subset1 <- data[which(data[,variable_order[1]] == i),]
      
      order_vec <- c(order_vec, unique(as.character(data_subset1[order(data_subset1[, variable_order[2]], decreasing = decreasing), variable_to_sort])) )
    }
  }
  
  data[, paste0(variable_to_sort, "sorted_", paste(variable_order, collapse = "_"))] <- factor(data[, variable_to_sort], levels = order_vec )
  return(data)
  #return(data[ match(levels(data[,paste0(variable_to_sort, "sorted_", paste(variable_order, collapse = "_"))]), data[, variable_to_sort] ),])
  
}



### V1 plot creation function ###
#ANCESTRY: for each model (k = 2, 3, etc.), individuals are sorted using the following hierarchical procedure:
#1. primary ancestry group (the ancestry group with the largest proportion)
#2. (within each primary ancestry group) proportion of primary ancestry group (largest proportion --> small proportion)

#LOCATION_ANCESTRY: for each model (k = 2, 3, etc.), individuals are sorted using the following hierarchical procedure:
#1. location (N --> S)
#2. primary ancestry group (the ancestry group with the largest proportion)
#3. (within each primary ancestry group) proportion of primary ancestry group (largest proportion --> small proportion)

#arguments:
#entropy_list: list of entropy results dataframes
#q_col_pal: a vector of colors that has the same number of colors as the maximum k
#show_progress: if TRUE, will print a message after the plots for each model in entropy_list are complete

create_V1_plots <- function(entropy_list, q_col_pal, show_progress = TRUE) {
  
  V1_q_plot_list <- list()
  
  subset <- entropy_list[[1]][which(entropy_list[[1]]$cluster == "V1"),]
  #subset <- entropy_list[[1]]$q_long[which(entropy_list[[1]]$cluster == "V1"),]
  
  #https://chartio.com/resources/tutorials/how-to-sort-a-data-frame-by-multiple-columns-in-r/
  # Sort by vector name [z] then [x]
  #dataframe[ with(dataframe, order(z, x)),]
  subset$sample_code_sortedlocation <- factor(subset$sample_code, levels = as.character(subset[with(subset, order(location_plotting, decreasing = FALSE)), ]$sample_code) )
  
  #subset$sample_code_sortedlocation <- factor(metadata_subset_reorder$sample_code, levels = as.character(metadata_subset_reorder[with(metadata_subset_reorder, order(location_plotting, decreasing = FALSE)), ]$sample_code) )
  subset_sorted_by_location <- subset[match(levels(subset$sample_code_sortedlocation), subset$sample_code),]
  
  divide_between_locations_vec <- c()
  
  for (loc in levels(subset_sorted_by_location$location_plotting)[- (length(levels(subset_sorted_by_location$location_plotting))) ]) {
    
    divide_between_locations_vec <- c(divide_between_locations_vec, max(which(subset_sorted_by_location$location_plotting == loc)) )
  }
  
  
  location_vector <- c()
  lower_pos <- c()
  upper_pos <- c()
  average_pos <- c()
  
  for (loc in levels(subset_sorted_by_location$location_plotting)) {
    location_vector <- c(location_vector, loc)
    
    lower_pos_val <- min(which(subset_sorted_by_location$location_plotting == loc))
    upper_pos_val <- max(which(subset_sorted_by_location$location_plotting == loc))
    average_pos_val <- mean(c(lower_pos_val, upper_pos_val))
    
    lower_pos <- c(lower_pos, lower_pos_val )
    upper_pos <- c(upper_pos, upper_pos_val )
    average_pos <- c(average_pos, average_pos_val)
  }
  
  
  location_position_dataframe <- data.frame(location = location_vector,
                                            lower_pos = lower_pos,
                                            upper_pos = upper_pos, 
                                            average_pos = average_pos)
  location_position_dataframe$label = as.character(1:nrow(location_position_dataframe))
  
  sample_order_list <- list()
  
  for (dat in 1:length(entropy_list) ) {
    
    k_mod_num <- dat + 1
    
    entropy_list[[dat]]$major_genetic_group <- NA
    for (indiv in unique(entropy_list[[dat]]$sample_code)) {
      
      data_subset <- entropy_list[[dat]][which(entropy_list[[dat]]$sample_code == indiv),]
      
      entropy_list[[dat]]$major_genetic_group[which(entropy_list[[dat]]$sample_code == indiv)] <- data_subset[which(data_subset$proportion == max(data_subset$proportion) ),]$cluster
    }
    
    sample_order_majorgeneticgroup <- c()
    for (i in unique(entropy_list[[dat]]$major_genetic_group)) {
      
      data_subset <- entropy_list[[dat]][which( (entropy_list[[dat]]$major_genetic_group == i) & (entropy_list[[dat]]$cluster == i) ),]
      
      sample_order_majorgeneticgroup <- c(sample_order_majorgeneticgroup, as.character(data_subset[order(data_subset$proportion, decreasing = TRUE),]$sample_code) )
    }
    
    
    ### sorted by: 1. major ancestry group, 2. major ancestry group proportion (descending) ###
    genetic_group <- hierarchical_ordering_func(data = entropy_list[[dat]], 
                                                variable_order = c("major_genetic_group", "proportion"), 
                                                level_order = list(paste0("V", 1:k_mod_num)), 
                                                variable_to_sort = 'sample_code', 
                                                decreasing = TRUE)
    
    dataframe_sorted_by_ancestry <- genetic_group[match(levels(genetic_group$sample_codesorted_major_genetic_group_proportion), genetic_group$sample_code),]
    
    sample_order_list[[paste0("K", k_mod_num, "_ancestry")]] <- levels(genetic_group$sample_codesorted_major_genetic_group_proportion)
    
    
    if (k_mod_num == 2) {
      dividing_line_vec <- min(which(dataframe_sorted_by_ancestry$proportion < 0.5)) - 0.5
      
    } else {
      
      dividing_line_vec <- c()
      for (i in unique(dataframe_sorted_by_ancestry$major_genetic_group)[-1] ) {
        dividing_line_vec <- c(dividing_line_vec, min(which(dataframe_sorted_by_ancestry$major_genetic_group == i)) )
      }
    }
    
    
    V1_q_plot_list[[paste0("K_", k_mod_num, "genetic")]] <- ggplot() + 
      geom_bar(data = genetic_group, aes(fill= cluster, y=proportion, x=sample_codesorted_major_genetic_group_proportion,  color = cluster), stat="identity", width=1) +
      scale_fill_manual(name = "Species", values = q_col_pal[1:k_mod_num]) +
      scale_color_manual(name = "Species", values = q_col_pal[1:k_mod_num]) +
      theme(#legend.position = "none", 
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"), 
        axis.text.x=element_blank(),
        axis.text.y=element_text(size = 15), 
        axis.ticks.x=element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(5, 8, 2, 2, unit = "pt"),
        legend.position = "none") +
      geom_segment(data = data.frame(x1 = dividing_line_vec,
                                     x2 = dividing_line_vec, 
                                     y1 = 0, 
                                     y2 = 1 ), 
                   aes(x = x1, y = y1, xend = x2, yend = y2) ) +
      scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))
    
    
    
    ### sorted by: 1. Location, 2. major ancestry group, 3. major ancestry group proportion (descending) ###
    location_genetic_group <- hierarchical_ordering_func(data = entropy_list[[dat]], 
                                                         variable_order = c("location_plotting", "major_genetic_group", "proportion"), 
                                                         level_order = list(levels(entropy_list[[dat]]$location_plotting), paste0("V", 1:k_mod_num)), 
                                                         variable_to_sort = 'sample_code', 
                                                         decreasing = TRUE)
    
    sample_order_list[[paste0("K", k_mod_num, "_location_ancestry")]] <- levels(location_genetic_group$sample_codesorted_location_plotting_major_genetic_group_proportion)
    
    #info on adding text outside of plot: https://github.com/tidyverse/ggplot2/issues/3400 (used for indicating location identity)
    V1_q_plot_list[[paste0("K_", k_mod_num, "location_genetic")]] <- ggplot() + 
      geom_bar(data = location_genetic_group, aes(fill= cluster, y=proportion, x=sample_codesorted_location_plotting_major_genetic_group_proportion, color = cluster), 
               stat="identity", 
               width=1) +
      scale_fill_manual(name = "Species", values = q_col_pal[1:k_mod_num]) +
      scale_color_manual(name = "Species", values = q_col_pal[1:k_mod_num]) +
      theme(#legend.position = "none", 
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"), 
        axis.text.x=element_blank(),
        axis.text.y=element_text(size = 15), 
        axis.ticks.x=element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(5, 8, 2, 2, unit = "pt"),
        legend.position = "none") +
      geom_segment(data = data.frame(x1 = divide_between_locations_vec + 0.5, 
                                     x2 = divide_between_locations_vec + 0.5, 
                                     y1 = rep(0, length(divide_between_locations_vec)), 
                                     y2 = rep(1, length(divide_between_locations_vec))), 
                   aes(x = x1, y = y1, xend = x2, yend = y2) ) +
      scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
      geom_text(aes(x = average_pos, y = Inf, label = label), size = 7, data = location_position_dataframe, vjust = 0.1) +
      coord_cartesian(clip = "off")
    
    
    
    if (show_progress == TRUE)  print(paste0("finished K = ", k_mod_num, " model plots"))
    
  }
  
  
  sample_order_list_collapsed <- as.data.frame(do.call(cbind, sample_order_list))
  
  
  
  return(list(plot_list = V1_q_plot_list,
              sample_order_list = sample_order_list_collapsed))
  
  
}

#checking to make sure that order of individuals makes sense
#metadata_subset_reorder[match( as.character(V1_plots$sample_order_list$K3_location_ancestry ), metadata_subset_reorder$sample_code),]$location_plotting
#metadata_subset_reorder[match( as.character(V1_plots$sample_order_list$K6_ancestry), metadata_subset_reorder$sample_code),]$species2


### V2 plot creation function ###: all models with k > 2 are sorted based on the sorting of the k = 2 model
#ANCESTRY:
#first, the k = 2 model is sorted using the following procedure:
#1. primary ancestry group (the ancestry group with the largest proportion)
#2. (within each primary ancestry group) proportion of primary ancestry group (largest proportion --> small proportion)

#second, each model with k > 2 is sorted based on the sorting of the k = 2 model


#LOCATION_ANCESTRY: all models with k > 2 are sorted based on the sorting of the k = 2 model
#first, the k = 2 model is sorted using the following procedure:
#1. location (N --> S)
#2. primary ancestry group (the ancestry group with the largest proportion)
#3. (within each primary ancestry group) proportion of primary ancestry group (largest proportion --> small proportion)

#second, each model with k > 2 is sorted based on the sorting of the k = 2 model

#arguments:
#entropy_list: list of entropy results dataframes
#q_col_pal: a vector of colors that has the same number of colors as the maximum k
#show_progress: if TRUE, will print a message after the plots for each model in entropy_list are complete

create_V2_plots <- function(entropy_list, q_col_pal, show_progress = TRUE) {
  
  V2_q_plot_list <- list()
  
  
  subset <- entropy_list[[1]][which(entropy_list[[1]]$cluster == "V1"),]
  
  #https://chartio.com/resources/tutorials/how-to-sort-a-data-frame-by-multiple-columns-in-r/
  # Sort by vector name [z] then [x]
  #dataframe[ with(dataframe, order(z, x)),]
  subset$sample_code_sortedlocation <- factor(subset$sample_code, levels = as.character(subset[with(subset, order(location_plotting, decreasing = FALSE)), ]$sample_code) )
  
  #subset$sample_code_sortedlocation <- factor(metadata_subset_reorder$sample_code, levels = as.character(metadata_subset_reorder[with(metadata_subset_reorder, order(location_plotting, decreasing = FALSE)), ]$sample_code) )
  subset_sorted_by_location <- subset[match(levels(subset$sample_code_sortedlocation), subset$sample_code),]
  
  divide_between_locations_vec <- c()
  
  for (loc in levels(subset_sorted_by_location$location_plotting)[- (length(levels(subset_sorted_by_location$location_plotting))) ]) {
    
    divide_between_locations_vec <- c(divide_between_locations_vec, max(which(subset_sorted_by_location$location_plotting == loc)) )
  }
  
  location_vector <- c()
  lower_pos <- c()
  upper_pos <- c()
  average_pos <- c()
  
  for (loc in levels(subset_sorted_by_location$location_plotting)) {
    location_vector <- c(location_vector, loc)
    
    lower_pos_val <- min(which(subset_sorted_by_location$location_plotting == loc))
    upper_pos_val <- max(which(subset_sorted_by_location$location_plotting == loc))
    average_pos_val <- mean(c(lower_pos_val, upper_pos_val))
    
    lower_pos <- c(lower_pos, lower_pos_val )
    upper_pos <- c(upper_pos, upper_pos_val )
    average_pos <- c(average_pos, average_pos_val)
  }
  
  
  location_position_dataframe <- data.frame(location = location_vector,
                                            lower_pos = lower_pos,
                                            upper_pos = upper_pos, 
                                            average_pos = average_pos)
  #location_position_dataframe$label = as.character(1:nrow(location_position_dataframe))
  location_position_dataframe$label <- levels(subset_sorted_by_location$location_plotting)
  
  
  
  
  entropy_list[[1]]$major_genetic_group <- NA
  for (indiv in unique(entropy_list[[1]]$sample_code)) {
    
    data_subset <- entropy_list[[1]][which(entropy_list[[1]]$sample_code == indiv),]
    
    entropy_list[[1]]$major_genetic_group[which(entropy_list[[1]]$sample_code == indiv)] <- data_subset[which(data_subset$proportion == max(data_subset$proportion) ),]$cluster
  }
  
  ### sorted by: 1. major ancestry group, 2. major ancestry group proportion (descending) ###
  K2_genetic_group <- hierarchical_ordering_func(data = entropy_list[[1]], 
                                                 variable_order = c("major_genetic_group", "proportion"), 
                                                 level_order = list(paste0("V", 1:2)), 
                                                 variable_to_sort = 'sample_code', 
                                                 decreasing = TRUE)
  
  k2_dataframe_sorted_by_ancestry <- K2_genetic_group[match(levels(K2_genetic_group$sample_codesorted_major_genetic_group_proportion), K2_genetic_group$sample_code),]
  
  
  ### sorted by: 1. Location, 2. major ancestry group, 3. major ancestry group proportion (descending) ###
  K2_location_genetic_group <- hierarchical_ordering_func(data = entropy_list[[1]], 
                                                          variable_order = c("location_plotting", "major_genetic_group", "proportion"), 
                                                          level_order = list(levels(entropy_list[[1]]$location_plotting), paste0("V", 1:2)), 
                                                          variable_to_sort = 'sample_code', 
                                                          decreasing = TRUE)
  
  
  
  for (dat in 1:length(entropy_list) ) {
    
    k_mod_num <- dat + 1
    
    
    #genetic_group_alt <- hierarchical_ordering_func(data = entropy_list[[dat]], 
    #                                                variable_order = c("major_genetic_group", "proportion"), 
    #                                                level_order = list(paste0("V", 1:k_mod_num)), 
    #                                                variable_to_sort = 'sample_code', 
    #                                                decreasing = TRUE)
    
    genetic_group_alt <- entropy_list[[dat]]
    
    genetic_group_alt$sample_codesorted_genetic_K2 <- factor(genetic_group_alt$sample_code, levels = levels(K2_genetic_group$sample_codesorted_major_genetic_group_proportion))
    genetic_group_alt$sample_codesorted_location_genetic_K2 <- factor(genetic_group_alt$sample_code, levels = levels(K2_location_genetic_group$sample_codesorted_location_plotting_major_genetic_group_proportion))
    
    
    V2_q_plot_list[[paste0("K_", k_mod_num, "genetic")]] <- ggplot() + 
      geom_bar(data = genetic_group_alt, aes(fill= cluster, y=proportion, x=sample_codesorted_genetic_K2,  color = cluster), stat="identity", width=1) +
      scale_fill_manual(name = "Species", values = q_col_pal[1:k_mod_num]) +
      scale_color_manual(name = "Species", values = q_col_pal[1:k_mod_num]) +
      theme(#legend.position = "none", 
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"), 
        axis.text.x=element_blank(),
        axis.text.y=element_text(size = 15),
        axis.ticks.x=element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(5, 8, 2, 2, unit = "pt"),
        legend.position = "none") +
      geom_abline(slope = 0, intercept = c(0.25, 0.5, 0.75),  col = "#e5e5e5", linetype = 'dashed') +
      geom_segment(data = data.frame(x1 = min(which(k2_dataframe_sorted_by_ancestry$proportion < 0.5)) - 0.5, 
                                     x2 = min(which(k2_dataframe_sorted_by_ancestry$proportion < 0.5)) - 0.5, 
                                     y1 = 0, 
                                     y2 = 1 ), 
                   aes(x = x1, y = y1, xend = x2, yend = y2), color = 'black', linetype = 'dashed') +
      scale_y_continuous(breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1"))
    
    
    
    #info on adding text outside of plot: https://github.com/tidyverse/ggplot2/issues/3400 (used for indicating location identity)
    V2_q_plot_list[[paste0("K_", k_mod_num, "location_genetic")]] <- ggplot() + 
      geom_bar(data = genetic_group_alt, aes(fill= cluster, y=proportion, x=sample_codesorted_location_genetic_K2, color = cluster), 
               stat="identity", 
               width=1) +
      scale_fill_manual(name = "Species", values = q_col_pal[1:k_mod_num]) +
      scale_color_manual(name = "Species", values = q_col_pal[1:k_mod_num]) +
      theme(#legend.position = "none", 
        plot.background = element_rect(fill = "white"), 
        panel.background = element_rect(fill = "white"), 
        axis.text.x=element_blank(),
        axis.text.y=element_text(size = 15),
        axis.ticks.x=element_blank(),
        axis.title = element_blank(),
        plot.margin = margin(5, 8, 2, 2, unit = "pt"),
        legend.position = "none") +
      geom_abline(slope = 0, intercept = c(0.25, 0.5, 0.75),  col = "#e5e5e5", linetype = 'dashed') +
      geom_segment(data = data.frame(x1 = divide_between_locations_vec + 0.5, 
                                     x2 = divide_between_locations_vec + 0.5, 
                                     y1 = rep(0, length(divide_between_locations_vec)),
                                     y2 = rep(1, length(divide_between_locations_vec))), 
                   aes(x = x1, y = y1, xend = x2, yend = y2), color = 'black', linetype = 'dashed') +
      scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = c("0", "0.25", "0.5", "0.75", "1")) +
      #geom_text(aes(x = average_pos, y = Inf, label = label), size = 7, data = location_position_dataframe, vjust = 0.1) +
      geom_text(aes(x = average_pos, y = -Inf, label = label), size = 7, data = location_position_dataframe, vjust = 0.95) +
      coord_cartesian(clip = "off")
    
    if (show_progress == TRUE) print(paste0("finished K = ", k_mod_num, " model plots"))
  }
  
  return(list(plot_list = V2_q_plot_list,
              sample_order_list = data.frame(K2_major_gen_group = levels(K2_genetic_group$sample_codesorted_major_genetic_group_proportion),
                                             K2_location_gen = levels(K2_location_genetic_group$sample_codesorted_location_plotting_major_genetic_group_proportion)) ) )
}


#checking to make sure that order of individuals makes sense
#metadata_subset_reorder[match( as.character(V2_plots$sample_order_list$K2_location_gen), metadata_subset_reorder$sample_code),]$location_plotting
#metadata_subset_reorder[match( as.character(V2_plots$sample_order_list$K2_major_gen_group), metadata_subset_reorder$sample_code),]$species2

