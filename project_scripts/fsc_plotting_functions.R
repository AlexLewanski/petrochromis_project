#######################################################################
#######################################################################
### FUNCTIONS FOR EVALUATING FIT OF FASTSIMCOAL2 DEMOGRAPHIC MODELS ###
#######################################################################
#######################################################################


##############################################
### CUSTOM FUNCTIONS USED IN parse_fsc_par ###
##############################################

#eval_numeric --> determine whether a value can be coerced to numeric  (returns logical)
eval_numeric <- function(x, factor2character = TRUE) {
  if (factor2character) x <- as.character(x)
  return(!is.na(suppressWarnings(as.numeric(x))))
}


#str_to_mat --> convert vector of matrix row strings (e.g. c("0 1 1", "3 1 7", "9 7 2")) to a matrix
#convert_numeric = TRUE results in matrix being converted to numeric if possible
str_to_mat <- function(str_vec, convert_numeric = TRUE, separator = NULL) {
  if (is.null(separator)) separator <- " "
  mat_step1 <- do.call(rbind, lapply(str_vec, function(x, separator) {
    #unlist1 <- unlist(strsplit(x, separator))
    unlist1 <- unlist(strsplit(x, "\t| ")) #add hoc change on 9/27
    unlist1[!grepl("^$", unlist1)] #added 9/13/21 to deal with >1 spaces between values in mig mat
  }, separator = separator) )
  
  if (convert_numeric && all(eval_numeric(mat_step1))) {
    return(apply(mat_step1, 2, as.numeric))
  } else {
    return(mat_step1)
  }
}


#extract_number_samples --> from par/tpl file, extract number of simulated samples optionally the separator between elements (" " or \t)
extract_number_samples <- function(par, eval_separator = TRUE) {
  num_samples <- par[2]
  sp.l <- unlist(strsplit(num_samples, split = ' '))
  tab.l <- unlist(strsplit(num_samples, split = '\t'))
  
  output_list <- list()
  
  if (length(sp.l) >= length(tab.l)) {
    output_list[['number_samples']] <- as.numeric(sp.l[1])
    if (eval_separator) output_list[['separator']] <- " "
  } else {
    output_list[['number_samples']] <- as.numeric(tab.l[1])
    if (eval_separator) output_list[['separator']] <- "\t"
  }
  
  return(output_list)
}


#extract_mig_mats --> extract migration matrices and optionally return position of final matrix row
extract_mig_mats <- function(par,
                             num_mig_mat,
                             header_line,
                             num_samples,
                             include_final_pos = FALSE,
                             separator) {
  
  #position of migration section + 1 (# matrices) + number of headers + number of lines taken up by previous matrices
  mig_pos_vec_head <- header_line + 1 + 1:num_mig_mat + num_samples*0:(num_mig_mat - 1)
  
  mig_mat_extract <- mapply(function(x, y, par_file) {par_file[x:y]},
                            mig_pos_vec_head + 1, #beginning positions of the matrices
                            mig_pos_vec_head + num_samples, #ending positions of the matrices
                            SIMPLIFY = FALSE,
                            MoreArgs = list(par_file = par))
  
  if (include_final_pos) {
    mig_mat_processed <- lapply(mig_mat_extract, function(x, separator) {
      str_to_mat(str_vec = x, convert_numeric = TRUE, separator = separator)
    }, separator = separator)
    
    return(list(mig_mat_list = mig_mat_processed,
                end_position = max(mig_pos_vec_head + num_samples) ))
  } else {
    return(lapply(mig_mat_extract, function(x, separator) {
      str_to_mat(str_vec = x, convert_numeric = TRUE, separator = separator)
    }, separator = separator))
  }
}


#extract_hist_events --> extract historical events as a dataframe (each row represents an event)
extract_hist_events <- function(par, hist_header_pos, separator) {
  num_hist_events <- as.numeric(unlist(strsplit(par[hist_header_pos + 1], split = separator))[1])
  
  if (num_hist_events > 0) {
    hist_events <- par[hist_header_pos + 2:(num_hist_events + 1)]
    hist_event_step1 <- as.data.frame(do.call(rbind, strsplit(hist_events, split = separator)), stringsAsFactors = FALSE) #stringsAsFactors = FALSE is necessary for R versions (<4) where default is factor instead of character
    hist_event_step1[, 7] <- gsub('nomig', '-1', hist_event_step1[, 7])
    
    hist_event_step2 <- cbind.data.frame(lapply(hist_event_step1, function(x) {
      if (all(eval_numeric(x))) {
        return(as.numeric(x))
      } else {
        return(x)
      }
    }))
    
    colnames(hist_event_step2) <- c('time', 'source', 'sink', 'migrants', 'new_deme_size', 'new_growth_rate', 'mig_mat_index')
    
  } else {
    hist_event_step2 <- 'NONE'
  }
  
  return(list(num_hist_events = num_hist_events,
              processed_hist_events = hist_event_step2) )
}



###############################################################
### parse_fsc_par: PROCESSS TPL/PAR FILE FROM FASTSIMCOAL2  ###
###############################################################

### DIRECTIONS FOR LOADING FILE: ###
#parse_fsc_par was built to work with files read in using one of the following approaches:
#tpl_file <- readLines('/file/path/model.tpl', skip=0,)
#tpl_file <- scan('/file/path/model.tpl', character(0), sep = "\n", strip.white = TRUE)

#Other functions that resulted in the format for the data once it is read in should also work.

### parse_fsc_par arguments ###
#par: par or tpl file uploaded with scan or readLines (or equivalent method)
#start (default = 2): the position right before population sizes header (equivalently, the line
#                     specifying the number of populations)
#   - unless you have a good reason to switch this (e.g. there are additional lines that you added
#     to the beginning of the file), KEEP THIS VALUE AT 2.

parse_fsc_par <- function(par, start = 2) {
  
  #convert file to character if it is input as a factor
  if (is.factor(par)) par <- as.character(par)
  
  #extract information on number of samples and the separator (" " or \t) that separates elements on a line
  starting_info <- extract_number_samples(par = par, eval_separator = TRUE)
  
  ## EXTRACT THE FIRST SET OF ELEMENTS (THROUGH THE NUMBER OF MIGRATION MATRICES) ###
  #calculate header positions based on number of samples (used to locate elements)
  #starting position + number of headers + number of samples * number of times samples are involved (counting the number of previous rows of params)
  header_location <- start + 1:4 + starting_info$number_samples*0:3
  names(header_location) <- c('pop_size', 'samp_info', 'growth_rate', 'num_mig_mat')
  
  #create vectors specifying starting and ending locations for each element
  start_location <- header_location + 1 #starting locations of elements
  end_location <- header_location + ifelse(names(header_location) == 'num_mig_mat', 1, starting_info$number_samples) #ending locations of elements
  
  #extract elements
  initial_extract <- mapply(function(x, y, PAR) {return(PAR[x:y])},
                            start_location,
                            end_location,
                            MoreArgs = list(PAR = par),
                            SIMPLIFY = FALSE)
  
  
  ### PROCESS THE INITIAL LIST OF EXTRACTED ELEMENTS ###
  #convert sample info (sample size and optionally sample age and inbreeding) into a matrix
  initial_extract[['samp_info']] <- str_to_mat(initial_extract[['samp_info']], convert_numeric = TRUE, separator = starting_info$separator)
  colnames(initial_extract[['samp_info']]) <- c('sample_size', 'sample_age', 'inbreeding')[seq_len(ncol(initial_extract[['samp_info']]))]
  
  #convert elements to numeric if they meet to criteria:
  #(1) coercible to numeric
  #(2) not already numeric (this prevents as.numeric from collapsing the samp_info matrix into a vector)
  extract_set1 <- lapply(initial_extract, function(x) {
    if (all(eval_numeric(x)) && !is.numeric(x)) {
      return(as.numeric(x))
    } else {
      return(x)
    }
  })
  
  
  ### EXTRACT THE REMAINING ELEMENTS (MIGRATION MATRICES AND HISTORICAL EVENTS) ###
  #extract migration matrices
  if (extract_set1[['num_mig_mat']] > 0) {
    mig_mat_extract <- extract_mig_mats(par = par,
                                        num_mig_mat = extract_set1[['num_mig_mat']],
                                        header_line = header_location['num_mig_mat'],
                                        num_samples = starting_info$number_samples,
                                        include_final_pos = TRUE,
                                        separator = starting_info$separator)
    mig_mat_output <- mig_mat_extract$mig_mat_list
    hist_event_pos <- mig_mat_extract$end_position + 1 #position of historical event header when matrices are present
  } else {
    mig_mat_output <- 'NONE'
    hist_event_pos <- header_location['num_mig_mat'] + 2 #position of historical event header when matrices are absent
  }
  
  #extract historical events
  hist_event_extract <- extract_hist_events(par = par, 
                                            hist_header_pos = hist_event_pos, 
                                            separator = starting_info$separator)
  
  
  ### OUTPUT THE PROCESSED INFORMATION ###
  return(c(number_samples = starting_info$number_samples,
           extract_set1, 
           list(mig_mat_list = mig_mat_output), 
           hist_event_info = list(hist_event_extract) ))
}

generic_par <- function(par, generic_element = c('pop_size', 'mig_mat_list', 'hist_event_info')) {
  
  for (element in names(par)) {
    
    if (element == 'pop_size' && element %in% generic_element) {
      par[[element]] <- rep(10, length(par[[element]]))
    }
    
    if (element == 'mig_mat_list' && element %in% generic_element) {
      par[[element]] <- lapply(par[[element]], function(y) {
        y[y != 0] <- 10
        return(y)
      })
    }
    
    if (element == 'hist_event_info' && element %in% generic_element) {
      par[[element]]$processed_hist_events$time <- seq(10, length(par[[element]]$processed_hist_events$time)*10, by = 10)
    }
    
  }
  
  return(par)
}



time_interval <- function(mig_mat_list, historical_event, diverge_cushion) {
  time_df_list <- list(data.frame(mig_mat = 0, start = 0, end = NA))
  
  new_migmat_pos <- which(!duplicated(historical_event$mig_mat_index) & historical_event$mig_mat_index != 0)
  
  if (length(new_migmat_pos) == 0) {
    time_df_list[[1]][1, 'end']<- historical_event[historical_event$source != historical_event$sink,]$time - diverge_cushion
  } else {
    subset_historical_events <- historical_event[new_migmat_pos,]
    for (i in seq_len(nrow(subset_historical_events))) {
      if (subset_historical_events$source[i] == subset_historical_events$sink[i]) {
        time_df_list[[i]]$end <- subset_historical_events$time[i]
        time_df_list <- c(time_df_list, list(data.frame(mig_mat = i, start = subset_historical_events$time[i], end = NA)))
      } else {
        time_df_list[[i]]$end <- subset_historical_events$time[i] - diverge_cushion
        time_df_list <- c(time_df_list, list(data.frame(mig_mat = i, start = subset_historical_events$time[i] - diverge_cushion, end = NA)))
      }
    }
    
  }
  
  return(time_df_list)
}


migmat_arrow_info <- function(mig_mat,
                              start_y, 
                              end_y,
                              offset = 0,
                              arrow_replicates = 2,
                              max_val = NULL,
                              max_arrow_size = 2,
                              min_arrow_size = 0.2,
                              max_arrow_head_size = 0.15,
                              min_arrow_head_size = 0.04,
                              pop0_x,
                              pop1_x,
                              direction = c('forward', 'backward')) {
  
  if (all(mig_mat == 0) | any(mig_mat == 'NONE') ) {
    return(data.frame(from = 0, to = 0, mig_est = 0, arrow_size = 0, arrow_head_size = 0, arrow_pos = 0, arrow_start = 0, arrow_end = 0))
  }
  
  ### 1. argument checks ###
  if (start_y >= end_y)
    stop('start_y is more recent than end_y (migmat_arrow_info has a present --> past perspective) so start_y needs to be smaller than end_y')
  direction <- match.arg(direction, several.ok = FALSE)
  
  
  ### 2. process migration matrix ###
  dimnames(mig_mat) <- replicate(2, 0:(ncol(mig_mat) - 1), simplify = FALSE) #add row and column names 
  
  #transform matrix to a "long" dataframe where each row represents an entry in the matrix (and remove zero entries)
  melted_mig_mat <- melt(mig_mat) %>% 
    filter(value > 0)
  #Var1 --> from; Var2 --> to; value --> migration param
  colnames(melted_mig_mat) <- c('from', 'to', 'mig_est')
  
  mig_params <- unique(melted_mig_mat$mig_est[melted_mig_mat$mig_est > 0]) #vector of unique non-zero migration parameters
  number_unique_arrows <- length(mig_params) #how many unique migration parameters?
  
  
  ### 3. arrow size specifications ###
  if (is.null(max_val)) max_val <- max(mig_params) #set max value if one isn't provided
  
  #arrow size
  melted_mig_mat$arrow_size <- (melted_mig_mat$mig_est/max_val)*max_arrow_size #standardize 
  melted_mig_mat$arrow_size[melted_mig_mat$arrow_size < min_arrow_size] <- min_arrow_size
  
  #arrow head size
  arrow_head_size_dif <- max_arrow_head_size - min_arrow_head_size
  mig_mat_prop <- arrow_head_size_dif*(melted_mig_mat$mig_est/max_val)
  melted_mig_mat$arrow_head_size <- min_arrow_head_size + mig_mat_prop
  
  
  ### 4. arrow positions ###
  start_end_dif <- abs(end_y - start_y)
  updated_start <- start_y + (start_end_dif*offset)
  updated_end <- end_y - (start_end_dif*offset)
  updated_dif <- abs(updated_end - updated_start)
  interval_vec <- seq(updated_start, updated_end, by = abs(updated_start - updated_end)/((arrow_replicates*nrow(melted_mig_mat)) + 1))
  melted_mig_mat_updated <- do.call(rbind, replicate(arrow_replicates, melted_mig_mat, simplify = FALSE)) %>%
    #mutate(arrow_pos = seq(updated_start, updated_end, length.out = nrow(melted_mig_mat)*arrow_replicates))
    mutate(arrow_pos = interval_vec[-c(1, length(interval_vec))])
  
  
  ### 5. arrow directionality -- forward or backward in time representation ###
  if (direction == 'forward') {
    melted_mig_mat_updated <- melted_mig_mat_updated %>% 
      mutate(arrow_end = ifelse(from == 0, pop0_x, pop1_x),
             arrow_start = ifelse(from == 0, pop1_x, pop0_x))
  } else {
    melted_mig_mat_updated <- melted_mig_mat_updated %>% 
      mutate(arrow_start = ifelse(from == 0, pop0_x, pop1_x),
             arrow_end = ifelse(from == 0, pop1_x, pop0_x))
  }
  
  
  ### 6. symmetric migration -- add arrow head on other side of arrow with symmetric migration ###
  if ( (length(unique(melted_mig_mat_updated$mig_est)) == 1) && (nrow(melted_mig_mat) > 1) )  {
    melted_mig_mat_opposite <- melted_mig_mat_updated %>% 
      mutate(arrow_end1 = arrow_start,
             arrow_start1 = arrow_end) %>% 
      select(from, to, mig_est, arrow_size, arrow_head_size, arrow_pos, arrow_end = arrow_end1, arrow_start = arrow_start1)
    
    melted_mig_mat_updated <- rbind(melted_mig_mat_updated, melted_mig_mat_opposite)
  }
  
  return(melted_mig_mat_updated)
}

rescale <- function(x,min,max) {
  (((x - min(data$pop_size)) * (max - min)) / diff(range(data$pop_size))) + min
}

rescale_alt <- function(x, min, max, pop_size_vec) {
  (((x - min(c(pop_size_vec, x) )) * (max - min)) / diff(range(c(pop_size_vec, x)))) + min
}

theme_demovis <- function() {
  theme(panel.background = element_rect(color = 'white', fill = 'white'),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line.y = element_line(color = 'black', size = 1.25, lineend = 'square'),
        axis.ticks.x = element_blank(),
        axis.ticks.length.y = unit(2.5, "mm"),
        axis.ticks.y = element_line(size = 1.25),
        axis.text.y = element_text(color = 'black', size = 11),
        legend.position = 'none')
}



####################################################################################
### plot_demographic_mod: plotting demographic model (without population growth) ###
####################################################################################

plot_demographic_mod <- function(par,
                                 colors = TRUE,
                                 pop_colors = NULL,
                                 scale_pops = TRUE,
                                 pop_color = 'gray',
                                 mig_band_color = "#bfbfbf",
                                 generic_plotting = FALSE,
                                 max_arrow_head_size = 0.15,
                                 min_arrow_head_size = 0.05,
                                 generic_arrow_head_size = 0.15,
                                 label_size = 7) {
  
  scale_pops <- scale_pops & !generic_plotting
  
  if (generic_plotting) {
    data <- generic_par(par, generic_element = c('pop_size', 'mig_mat_list', 'hist_event_info')[c(1, 3)])
  } else {
    data <- par
  }
  
  events <- data$hist_event_info$processed_hist_events
  
  if (!is.null(pop_colors) & colors & length(pop_colors) == data$number_samples) {
    pop.cols <- pop_colors
  } else {
    pop.cols <- viridis::viridis(data$number_samples)
  }
  
  ############################
  # subsetting data into migration vs div events
  ############################
  mig_events <- events %>%
    filter(migrants < 1) %>%
    filter(source != sink)
  div_events <- events %>%
    filter(migrants == 1) %>%
    filter(source != sink)
  
  #UNCOMMENTED THIS
  # if (scale_pops) {
  #   pops_rel_ne <- rescale(data$pop_size, 0.5, 1.2)
  #   #root_rel_ne <- rescale(div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[1], 0.5, 1.2)
  #   root_rel_ne <- rescale(div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[div_events$sink + 1], 0.5, 1.2)
  #   div_events$source_rel_ne <- pops_rel_ne[div_events$source + 1]
  #   div_events$sink_rel_ne <- pops_rel_ne[div_events$sink + 1]
  # 
  #   mig_events$source_rel_ne <- pops_rel_ne[mig_events$source + 1]
  #   mig_events$sink_rel_ne <- pops_rel_ne[mig_events$sink + 1]
  # }
  
  #COMMENTED THIS
  if (isTRUE(scale_pops)) {

    pops_rel_ne <- rescale_alt(data$pop_size, 0.5, 1.2, pop_size_vec = c(data$pop_size, div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[div_events$sink + 1]) )

    #root_rel_ne <- rescale(div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[1], 0.5, 1.2)
    root_rel_ne <- rescale_alt( div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[div_events$sink + 1], 0.5, 1.2, pop_size_vec = c(data$pop_size, div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[div_events$sink + 1]) )
    div_events$source_rel_ne <- pops_rel_ne[div_events$source + 1]
    div_events$sink_rel_ne <- pops_rel_ne[div_events$sink + 1]

    mig_events$source_rel_ne <- pops_rel_ne[mig_events$source + 1]
    mig_events$sink_rel_ne <- pops_rel_ne[mig_events$sink + 1]
  }

  
  
  ##########################################
  # first, set up the skeleton using divergences
  # it'd be good if we can make this more flexible (i.e. works with any number of div events)
  ##########################################
  tgap <- max(events$time)*0.10
  troot <- max(events$time)*1.33
  
  p <- ggplot() +
    #xlim(0,length(data$pop_size)*1.5) +
    ylim(0,troot)
  
  #width <- 0.75
  width <- 0.6
  if (scale_pops) {
    if (nrow(div_events) == 1){
      #UNCOMMENTED THIS
      # p1 <- p +
      #   geom_rect(mapping=aes(xmin=min(div_events$source[1],div_events$sink[1]) + 1 - width*pops_rel_ne[1]/2,
      #                         xmax=max(div_events$source[1],div_events$sink[1]) + 1 + width*pops_rel_ne[2]/2,
      #                         ymin=div_events$time[1],ymax=div_events$time[1]-tgap), fill=pop_color) +
      #   # root
      #   geom_rect(mapping=aes(xmin=(div_events$source[1]+div_events$sink[1]+2)/2 - width*root_rel_ne/2,
      #                         xmax=(div_events$source[1]+div_events$sink[1]+2)/2 + width*root_rel_ne/2,
      #                         ymin=div_events$time[1],
      #                         ymax=troot), fill=pop_color) +
      #   # pop0
      #   geom_rect(mapping=aes(xmin=div_events$source[1] + 1 - width*pops_rel_ne[1]/2,
      #                         xmax=div_events$source[1] + 1 + width*pops_rel_ne[1]/2,
      #                         ymin=0,
      #                         ymax=div_events$time[1]), fill=pop_color) +
      #   # pop1
      #   geom_rect(mapping=aes(xmin=div_events$sink[1] + 1 - width*pops_rel_ne[2]/2,
      #                         xmax=div_events$sink[1] + 1 + width*pops_rel_ne[2]/2,
      #                         ymin=0,
      #                         ymax=div_events$time[1]), fill=pop_color)
      
      #COMMENTED THIS
      # p1 <- p +
      #   geom_rect(mapping=aes(xmin=min(div_events$source[1],div_events$sink[1]) + 1 - width*pops_rel_ne[1]/2,
      #                         xmax=max(div_events$source[1],div_events$sink[1]) + 1 + width*pops_rel_ne[2]/2,
      #                         ymin=div_events$time[1],ymax=div_events$time[1]-tgap), fill=pop_color) +
      #   # root
      #   geom_rect(mapping=aes(xmin=(div_events$source[1]+div_events$sink[1]+2)/2 - width*root_rel_ne/2,
      #                         xmax=(div_events$source[1]+div_events$sink[1]+2)/2 + width*root_rel_ne/2,
      #                         ymin=div_events$time[1],
      #                         ymax=troot), fill=pop_color) +
      #   # pop0
      #   geom_rect(mapping=aes(xmin=div_events$source[1] + 0.9 - width*pops_rel_ne[1]/2,
      #                         xmax=div_events$source[1] + 0.9 + width*pops_rel_ne[1]/2,
      #                         ymin=0,
      #                         ymax=div_events$time[1]), fill=pop_color) +
      #   # pop1
      #   geom_rect(mapping=aes(xmin=div_events$sink[1] + 1.1 - width*pops_rel_ne[2]/2,
      #                         xmax=div_events$sink[1] + 1.1 + width*pops_rel_ne[2]/2,
      #                         ymin=0,
      #                         ymax=div_events$time[1]), fill=pop_color)
      
      #4/1/2022
      p1 <- p +
        geom_rect(mapping=aes(xmin=min(div_events$source[1],div_events$sink[1]) + 1 - width*pops_rel_ne[1]/2,
                              xmax=max(div_events$source[1],div_events$sink[1]) + 1 + width*pops_rel_ne[2]/2,
                              ymin=div_events$time[1],ymax=div_events$time[1]-tgap), fill=pop_color) +
        # root
        geom_rect(mapping=aes(xmin=(div_events$source[1]+div_events$sink[1]+2)/2 - width*root_rel_ne/2,
                              xmax=(div_events$source[1]+div_events$sink[1]+2)/2 + width*root_rel_ne/2,
                              ymin=div_events$time[1],
                              ymax=troot), fill=pop_color) +
        # pop0
        geom_rect(mapping=aes(xmin= 0 + 0.9 - width*pops_rel_ne[1]/2,
                              xmax= 0 + 0.9 + width*pops_rel_ne[1]/2,
                              ymin=0,
                              ymax=div_events$time[1]), fill=pop_color) +
        # pop1
        geom_rect(mapping=aes(xmin= 1 + 1.1 - width*pops_rel_ne[2]/2,
                              xmax= 1 + 1.1 + width*pops_rel_ne[2]/2,
                              ymin=0,
                              ymax=div_events$time[1]), fill=pop_color)
      
      #print(p1)
      
    } else if (nrow(div_events) == 2) { 
      p1 <- p + 
        geom_rect(mapping=aes(xmin=min(div_events$source[1],div_events$sink[1]) + 1 - width*pops_rel_ne[1]/2,
                              xmax=max(div_events$source[1],div_events$sink[1]) + 1 + width*pops_rel_ne[2]/2,
                              ymin=div_events$time[1],ymax=div_events$time[1]-tgap), fill=pop_color) +
        # root
        geom_rect(mapping=aes(xmin=(div_events$source[1]+div_events$sink[1]+2)/2 - width*root_rel_ne/2,
                              xmax=(div_events$source[1]+div_events$sink[1]+2)/2 + width*root_rel_ne/2,
                              ymin=div_events$time[1],
                              ymax=troot), fill=pop_color) + 
        # pop0
        geom_rect(mapping=aes(xmin=div_events$source[1] + 1 - width*pops_rel_ne[1]/2,
                              xmax=div_events$source[1] + 1 + width*pops_rel_ne[1]/2,
                              ymin=0,
                              ymax=div_events$time[1]), fill=pop_color) +
        # pop1
        geom_rect(mapping=aes(xmin=div_events$sink[1] + 1 - width*pops_rel_ne[2]/2,
                              xmax=div_events$sink[1] + 1 + width*pops_rel_ne[2]/2,
                              ymin=0,
                              ymax=div_events$time[1]), fill=pop_color) +
        
        geom_rect(mapping=aes(xmin=max(div_events$source[2],div_events$sink[2]) + 1 - width*pops_rel_ne[3]/2,
                              xmax=max(div_events$source[2],div_events$sink[2]) + 1 + width*pops_rel_ne[3]/2,
                              ymin=0,
                              ymax=div_events$time[2]), fill=pop_color) +
        
        geom_rect(mapping=aes(xmin=min(div_events$source[2],div_events$sink[2]) + 1 - width*pops_rel_ne[2]/2,
                              xmax=max(div_events$source[2],div_events$sink[2]) + 1 + width*pops_rel_ne[3]/2,
                              ymin=div_events$time[2]-tgap,
                              ymax=div_events$time[2]), fill=pop_color)
      #print(p1)
    }
  } else {
    if (nrow(div_events) == 1){
      p1 <- p + 
        geom_rect(mapping=aes(xmin=min(div_events$source[1],div_events$sink[1]) + 1 - width/2,
                              xmax=max(div_events$source[1],div_events$sink[1]) + 1 + width/2,
                              ymin=div_events$time[1],ymax=div_events$time[1]-tgap), fill=pop_color) + 
        # root
        geom_rect(mapping=aes(xmin=(div_events$source[1]+div_events$sink[1]+2)/2 - width/2,
                              xmax=(div_events$source[1]+div_events$sink[1]+2)/2 + width/2,
                              ymin=div_events$time[1],
                              ymax=troot), fill=pop_color) + 
        # pop0
        geom_rect(mapping=aes(xmin=div_events$source[1] + 1 - width/2,
                              xmax=div_events$source[1] + 1 + width/2,
                              ymin=0,
                              ymax=div_events$time[1]), fill=pop_color) +
        # pop1
        geom_rect(mapping=aes(xmin=div_events$sink[1] + 1 - width/2,
                              xmax=div_events$sink[1] + 1 + width/2,
                              ymin=0,
                              ymax=div_events$time[1]), fill=pop_color)
      #print(p1)
      
    } else if (nrow(div_events) == 2) { 
      p1 <- p + 
        geom_rect(mapping=aes(xmin=min(div_events$source[1],div_events$sink[1]) + 1 - width/2,
                              xmax=max(div_events$source[1],div_events$sink[1]) + 1 + width/2,
                              ymin=div_events$time[1],ymax=div_events$time[1]-tgap), fill=pop_color) +
        # root
        geom_rect(mapping=aes(xmin=(div_events$source[1]+div_events$sink[1]+2)/2 - width/2,
                              xmax=(div_events$source[1]+div_events$sink[1]+2)/2 + width/2,
                              ymin=div_events$time[1],
                              ymax=troot), fill=pop_color) + 
        # pop0
        geom_rect(mapping=aes(xmin=div_events$source[1] + 1 - width/2,
                              xmax=div_events$source[1] + 1 + width/2,
                              ymin=0,
                              ymax=div_events$time[1]), fill=pop_color) +
        # pop1
        geom_rect(mapping=aes(xmin=div_events$sink[1] + 1 - width/2,
                              xmax=div_events$sink[1] + 1 + width/2,
                              ymin=0,
                              ymax=div_events$time[1]), fill=pop_color) +
        
        geom_rect(mapping=aes(xmin=max(div_events$source[2],div_events$sink[2]) + 1 - width/2,
                              xmax=max(div_events$source[2],div_events$sink[2]) + 1 + width/2,
                              ymin=0,
                              ymax=div_events$time[2]), fill=pop_color) +
        
        geom_rect(mapping=aes(xmin=min(div_events$source[2],div_events$sink[2]) + 1 - width/2,
                              xmax=max(div_events$source[2],div_events$sink[2]) + 1 + width/2,
                              ymin=div_events$time[2]-tgap,
                              ymax=div_events$time[2]), fill=pop_color)
      #print(p1)
    }
  }
  
  
  if (scale_pops) {
    #UNCOMMENTED THIS
    # pop0_x_left <- div_events$source[1] + 1 - width*pops_rel_ne[1]/2
    # pop0_x_right <- div_events$source[1] + 1 + width*pops_rel_ne[1]/2
    # pop1_x_left <- div_events$sink[1] + 1 - width*pops_rel_ne[2]/2
    # pop1_x_right <-div_events$sink[1] + 1 + width*pops_rel_ne[2]/2
    
    #COMMENTED THIS
    # pop0_x_left <- div_events$source[1] + 0.9 - width*pops_rel_ne[1]/2
    # pop0_x_right <- div_events$source[1] + 0.9 + width*pops_rel_ne[1]/2
    # pop1_x_left <- div_events$sink[1] + 1.1 - width*pops_rel_ne[2]/2
    # pop1_x_right <-div_events$sink[1] + 1.1 + width*pops_rel_ne[2]/2
    
    #4/1/2022
    pop0_x_left <- 0 + 0.9 - width*pops_rel_ne[1]/2
    pop0_x_right <- 0 + 0.9 + width*pops_rel_ne[1]/2
    pop1_x_left <- 1 + 1.1 - width*pops_rel_ne[2]/2
    pop1_x_right <-1 + 1.1 + width*pops_rel_ne[2]/2
    
  } else {
    # pop0_x_left <- div_events$source[1] + 1 - width/2
    # pop0_x_right <- div_events$source[1] + 1 + width/2
    # pop1_x_left <- div_events$sink[1] + 1 - width/2
    # pop1_x_right <- div_events$sink[1] + 1 + width/2
    
    pop0_x_left <- 0 + 1 - width/2
    pop0_x_right <- 0 + 1 + width/2
    pop1_x_left <- 1 + 1 - width/2
    pop1_x_right <- 1 + 1 + width/2
  }
  
  
  ############################
  # add colors for populations, if desired
  ############################
  if (colors){
    if (scale_pops) {
      p1.2 <- p1 +
        geom_rect(data=div_events, mapping=aes(xmin=source + 1 - width*source_rel_ne/2,
                                               xmax=source + 1 + width*source_rel_ne/2,
                                               ymin=0,
                                               ymax=tgap/2, 
                                               fill=factor(source))) +
        geom_rect(data=div_events, mapping=aes(xmin=sink + 1 - width*sink_rel_ne/2,
                                               xmax=sink + 1 + width*sink_rel_ne/2,
                                               ymin=0,
                                               ymax=tgap/2,
                                               fill=factor(sink))) +
        scale_fill_manual(values=pop.cols)
    } else {
      p1.2 <- p1 +
        geom_rect(data=div_events, mapping=aes(xmin=source + 1 - width/2,
                                               xmax=source + 1 + width/2,
                                               ymin=0,
                                               ymax=tgap/2, 
                                               fill=factor(source))) +
        geom_rect(data=div_events, mapping=aes(xmin=sink + 1 - width/2,
                                               xmax=sink + 1 + width/2,
                                               ymin=0,
                                               ymax=tgap/2,
                                               fill=factor(sink))) +
        scale_fill_manual(values=pop.cols)
    }
    
  } else {
    p1.2 <- p1
  }
  
  
  ############################  
  # then, add migration events
  # needs to be extended for > 2 pops
  ############################
  if (nrow(mig_events) > 0) {
    if (length(pops) == 2) {
      p2 <- p1.2 +
        geom_segment(aes(x = if_else(mig_events$source < mig_events$sink,
                                     (mig_events$source*2+2)/2+0.25,
                                     (mig_events$source*2+2)/2-0.25), 
                         xend = if_else(mig_events$source < mig_events$sink,
                                        (mig_events$sink*2+2)/2-0.25,
                                        (mig_events$sink*2+2)/2+0.25),
                         y = mig_events$time,
                         yend = mig_events$time),
                     lineend = "round", # See available arrow types in example above
                     linejoin = "round",
                     size = 1.5, 
                     arrow = arrow(length = unit(0.1, "inches")),
                     colour = "black" # Also accepts "red", "blue' etc
        ) 
      print(p2)
    } else if (length(pops) == 3) {
      p2 <- p1.2 +
        geom_segment(data=mig_events[mig_events$time > div_events$time[2],], aes(x = if_else(source < sink,
                                                                                             (source*2+2)/2+0.25,
                                                                                             (source*2+2)/2-0.25), 
                                                                                 xend = if_else(source < sink,
                                                                                                (sink*2+2)/2-0.25,
                                                                                                (sink*2+2)/2+0.25),
                                                                                 y = time,
                                                                                 yend = time),
                     lineend = "round", # See available arrow types in example above
                     linejoin = "round",
                     size = 1.5, 
                     arrow = arrow(length = unit(0.1, "inches")),
                     colour = "black" # Also accepts "red", "blue' etc
        ) +
        geom_segment(data=mig_events[mig_events$time < div_events$time[2],], 
                     aes(x = if_else(source < sink,
                                     (source*2+2)/2+0.24,
                                     (source*2+2)/2-0.51), 
                         xend = if_else(source < sink,
                                        (sink*2+2)/2-0.51,
                                        (sink*2+2)/2+0.24),
                         y = time,
                         yend = time),
                     lineend = "round", # See available arrow types in example above
                     linejoin = "round",
                     size = 1.5, 
                     arrow = arrow(length = unit(0.1, "inches")),
                     colour = "black" # Also accepts "red", "blue' etc
        ) 
      #print(p2)
    }
  } else {
    p2 <- p1.2
  }
  
  events_sorted <- events %>%
    arrange((time))
  
  #MIGRATION MATRICES
  time_interval_list <- time_interval(mig_mat_list = data$mig_mat_list,
                                      historical_event = events_sorted,
                                      diverge_cushion = tgap)
  
  
  migmat_arrow_info_list <- lapply(split(seq_len(length(data$mig_mat_list)), seq_len(length(data$mig_mat_list))), function(i, mig_mat, time_int_list) {
    
    if (all(!is.na(time_int_list[[i]]$end))) {
      
      migmat_arrow_info(mig_mat = mig_mat[[i]],
                        start_y = time_int_list[[i]]$start,
                        #end_y = div_events$time[1],
                        end_y = time_int_list[[i]]$end,
                        offset = 0.05,
                        arrow_replicates = 1,
                        max_val = NULL,
                        max_arrow_size = 1,
                        min_arrow_size = 0.2,
                        max_arrow_head_size = max_arrow_head_size,
                        min_arrow_head_size = min_arrow_head_size,
                        pop0_x = pop0_x_right,
                        pop1_x = pop1_x_left,
                        direction = c('forward', 'backward')[2])
    }
    
  }, time_int_list = time_interval_list, mig_mat = data$mig_mat_list)
  
  migmat_arrow_info_list_filtered <- migmat_arrow_info_list[sapply(migmat_arrow_info_list, function(x) !is.null(x) & !all(x == 0))]
  
  
  if (generic_plotting) {
    migmat_arrow_info_list_filtered <- lapply(migmat_arrow_info_list_filtered, function(x, arrow_head_size) {
      x$arrow_size <- 1
      x$arrow_head_size <- arrow_head_size
      return(x)
    }, arrow_head_size = generic_arrow_head_size)
  }
  
  
  mig_mat_eval <- sapply(data$mig_mat_list, function(x) !all(x == 0) & x != 'NONE')
  if (any(mig_mat_eval)) {
    #migration_period <- time_interval_list[sapply(data$mig_mat_list, function(x) !all(x == 0) & x != 'NONE')][[1]]
    migration_period <- time_interval_list[sapply(data$mig_mat_list, function(x) !all(x == 0) & any(x != 'NONE'))][[1]]
    
    p3 <- p2 + 
      geom_rect(data = migration_period,
                aes(xmin = pop0_x_right, 
                    xmax = pop1_x_left, 
                    ymin = start, 
                    ymax = end), 
                color = NA, fill = mig_band_color, alpha = 1)
    
  } else {
    p3 <- p2
  }
  
  #migration_period <- time_interval_list[sapply(data$mig_mat_list, function(x) !all(x == 0))][[1]]
  
  
  # p3 <- p2 + 
  #   geom_rect(data = migration_period,
  #             aes(xmin = pop0_x_right, 
  #                 xmax = pop1_x_left, 
  #                 ymin = start, ymax = end), 
  #             color = NA, fill = "#bfbfbf", alpha=0.5)
  
  if (length(migmat_arrow_info_list_filtered) != 0) {
    for (i in seq_len(length(migmat_arrow_info_list_filtered))) {
      
      p3 <- p3 +
        geom_segment(aes(x = migmat_arrow_info_list_filtered[[i]]$arrow_start, 
                         xend = migmat_arrow_info_list_filtered[[i]]$arrow_end, 
                         y = migmat_arrow_info_list_filtered[[i]]$arrow_pos, # might not want to hard code this? or maybe it's fine?
                         yend = migmat_arrow_info_list_filtered[[i]]$arrow_pos), # might not want to hard code this? or maybe it's fine?
                     lineend = "round", # See available arrow types in example above
                     linejoin = "round",
                     size = migmat_arrow_info_list_filtered[[i]]$arrow_size, 
                     arrow = arrow(length = unit(migmat_arrow_info_list_filtered[[i]]$arrow_head_size, "inches"), type = "closed"),
                     lty=1,
                     colour = "black" # Also accepts "red", "blue' etc
        )
    }
  }
  
  
  if (generic_plotting) {
    return(p3 + theme_void() + theme(legend.position = 'none'))
  } else {
    
    return(
      p3 +
        theme_demovis() +
        ylab('Time (generations)') +
        theme(axis.title.y = element_text(size = 16)) +
        scale_x_continuous(expand = c(0.03, 0.03)) +
        scale_y_continuous(labels = function(x) format(x, scientific = TRUE), expand = c(0, 0)) +
        geom_segment(data=data$hist_event_info$processed_hist_events,
                     aes(x = pop0_x_left, xend = pop1_x_right,
                         y=time, yend=time),
                     lty=2, size = 1.25, lineend = "round") +
        geom_text(aes(x = pop1_x_right*1.05,
                      y = data$hist_event_info$processed_hist_events$time, label = formatC(data$hist_event_info$processed_hist_events$time, format = "d", big.mark = ",")), size = label_size, hjust = 0) +
        coord_cartesian(clip = 'off')
      #xlim(pop0_x_left, pop1_x_right*1.1)
    )
    
  }
  
}



########################################################################################
### ADDITIONAL FUNCTIONS FOR PLOTTING 2 POP MODELS WITH CONTINUOUS POPULATION GROWTH ###
########################################################################################

calc_historical_popsize <- function(current_pop, growth_rate, time) {
  return(current_pop * exp(growth_rate * time))
}

id_xval <- function(y_val, pointa, pointb) {
  slope <- (pointb[2] - pointa[2])/(pointb[1] - pointa[1])
  
  return(((y_val - pointa[2])/slope) + pointa[1])
}


migmat_arrow_info_alt <- function(mig_mat,
                                  start_y, 
                                  end_y,
                                  offset = 0,
                                  arrow_replicates = 2,
                                  max_val = NULL,
                                  max_arrow_size = 2,
                                  min_arrow_size = 0.2,
                                  max_arrow_head_size = 0.15,
                                  min_arrow_head_size = 0.04,
                                  pop0_x,
                                  pop1_x,
                                  direction = c('forward', 'backward')) {
  
  if (all(mig_mat == 0) | any(mig_mat == 'NONE') ) {
    return(data.frame(from = 0, to = 0, mig_est = 0, arrow_size = 0, arrow_head_size = 0, arrow_pos = 0, arrow_start = 0, arrow_end = 0))
  }
  
  ### 1. argument checks ###
  if (start_y >= end_y)
    stop('start_y is more recent than end_y (migmat_arrow_info has a present --> past perspective) so start_y needs to be smaller than end_y')
  direction <- match.arg(direction, several.ok = FALSE)
  
  
  ### 2. process migration matrix ###
  dimnames(mig_mat) <- replicate(2, 0:(ncol(mig_mat) - 1), simplify = FALSE) #add row and column names 
  
  #transform matrix to a "long" dataframe where each row represents an entry in the matrix (and remove zero entries)
  melted_mig_mat <- melt(mig_mat) %>% 
    filter(value > 0)
  #Var1 --> from; Var2 --> to; value --> migration param
  colnames(melted_mig_mat) <- c('from', 'to', 'mig_est')
  
  mig_params <- unique(melted_mig_mat$mig_est[melted_mig_mat$mig_est > 0]) #vector of unique non-zero migration parameters
  number_unique_arrows <- length(mig_params) #how many unique migration parameters?
  
  
  ### 3. arrow size specifications ###
  if (is.null(max_val)) max_val <- max(mig_params) #set max value if one isn't provided
  
  #arrow size
  melted_mig_mat$arrow_size <- (melted_mig_mat$mig_est/max_val)*max_arrow_size #standardize 
  melted_mig_mat$arrow_size[melted_mig_mat$arrow_size < min_arrow_size] <- min_arrow_size
  
  #arrow head size
  arrow_head_size_dif <- max_arrow_head_size - min_arrow_head_size
  mig_mat_prop <- arrow_head_size_dif*(melted_mig_mat$mig_est/max_val)
  melted_mig_mat$arrow_head_size <- min_arrow_head_size + mig_mat_prop
  
  
  ### 4. arrow positions ###
  start_end_dif <- abs(end_y - start_y)
  updated_start <- start_y + (start_end_dif*offset)
  updated_end <- end_y - (start_end_dif*offset)
  updated_dif <- abs(updated_end - updated_start)
  interval_vec <- seq(updated_start, updated_end, by = abs(updated_start - updated_end)/((arrow_replicates*nrow(melted_mig_mat)) + 1))
  melted_mig_mat_updated <- do.call(rbind, replicate(arrow_replicates, melted_mig_mat, simplify = FALSE)) %>%
    #mutate(arrow_pos = seq(updated_start, updated_end, length.out = nrow(melted_mig_mat)*arrow_replicates))
    mutate(arrow_pos = interval_vec[-c(1, length(interval_vec))])
  
  
  ####################################################################
  ### THIS PART NEEDS TO BE CHANGED FOR ADDRESSING POP SIZE CHANGE ###
  ####################################################################
  ### 5. arrow directionality -- forward or backward in time representation ###
  if (direction == 'forward') {
    # melted_mig_mat_updated <- melted_mig_mat_updated %>%
    #   mutate(arrow_end = ifelse(from == 0, pop0_x, pop1_x),
    #          arrow_start = ifelse(from == 0, pop1_x, pop0_x))
    melted_mig_mat_updated <- melted_mig_mat_updated %>%
      mutate(arrow_end = ifelse(from == 0,
                                id_xval(y_val = melted_mig_mat_updated$arrow_pos,
                                        pointa = pop0_x[['pointa']],
                                        pointb = pop0_x[['pointb']]),
                                id_xval(y_val = melted_mig_mat_updated$arrow_pos,
                                        pointa = pop1_x[['pointa']],
                                        pointb = pop1_x[['pointb']])),
             arrow_start = ifelse(from == 0,
                                  id_xval(y_val = melted_mig_mat_updated$arrow_pos,
                                          pointa = pop1_x[['pointa']],
                                          pointb = pop1_x[['pointb']]),
                                  id_xval(y_val = melted_mig_mat_updated$arrow_pos,
                                          pointa = pop0_x[['pointa']],
                                          pointb = pop0_x[['pointb']])))
  } else {
    # melted_mig_mat_updated <- melted_mig_mat_updated %>% 
    #   mutate(arrow_start = ifelse(from == 0, pop0_x, pop1_x),
    #          arrow_end = ifelse(from == 0, pop1_x, pop0_x))
    melted_mig_mat_updated <- melted_mig_mat_updated %>%
      mutate(arrow_start = ifelse(from == 0, 
                                  id_xval(y_val = melted_mig_mat_updated$arrow_pos, 
                                          pointa = pop0_x[['pointa']], 
                                          pointb = pop0_x[['pointb']]), 
                                  id_xval(y_val = melted_mig_mat_updated$arrow_pos, 
                                          pointa = pop1_x[['pointa']], 
                                          pointb = pop1_x[['pointb']])),
             arrow_end = ifelse(from == 0, 
                                id_xval(y_val = melted_mig_mat_updated$arrow_pos, 
                                        pointa = pop1_x[['pointa']], 
                                        pointb = pop1_x[['pointb']]), 
                                id_xval(y_val = melted_mig_mat_updated$arrow_pos, 
                                        pointa = pop0_x[['pointa']], 
                                        pointb = pop0_x[['pointb']]))) 
  }
  
  ####################################################################
  ####################################################################
  ####################################################################
  
  
  ### 6. symmetric migration -- add arrow head on other side of arrow with symmetric migration ###
  if ( (length(unique(melted_mig_mat_updated$mig_est)) == 1) && (nrow(melted_mig_mat) > 1) )  {
    melted_mig_mat_opposite <- melted_mig_mat_updated %>% 
      mutate(arrow_end1 = arrow_start,
             arrow_start1 = arrow_end) %>% 
      select(from, to, mig_est, arrow_size, arrow_head_size, arrow_pos, arrow_end = arrow_end1, arrow_start = arrow_start1)
    
    melted_mig_mat_updated <- rbind(melted_mig_mat_updated, melted_mig_mat_opposite)
  }
  
  return(melted_mig_mat_updated)
}








plot_demographic_mod_popgrow <- function(par,
                                         colors = TRUE,
                                         pop_colors = NULL,
                                         scale_pops = TRUE,
                                         pop_color = 'gray',
                                         generic_plotting = FALSE,
                                         max_arrow_head_size = 0.15,
                                         min_arrow_head_size = 0.05,
                                         generic_arrow_head_size = 0.15,
                                         label_size = 7,
                                         label_nudge = 0,
                                         perspective = c('forward', 'backward')) {
  
  scale_pops <- scale_pops & !generic_plotting
  
  if (generic_plotting) {
    data <- generic_par(par, generic_element = c('pop_size', 'mig_mat_list', 'hist_event_info')[c(1, 3)])
  } else {
    data <- par
  }
  
  events <- data$hist_event_info$processed_hist_events
  
  if (!is.null(pop_colors) & colors & length(pop_colors) == data$number_samples) {
    pop.cols <- pop_colors
  } else {
    pop.cols <- viridis::viridis(data$number_samples)
  }
  
  ############################
  # subsetting data into migration vs div events
  ############################
  mig_events <- events %>%
    filter(migrants < 1) %>%
    filter(source != sink)
  div_events <- events %>%
    filter(migrants == 1) %>%
    filter(source != sink)
  
  pop_size_list <- list(current = data$pop_size,
                        historical = calc_historical_popsize(current_pop = data$pop_size, 
                                                             growth_rate = data$growth_rate, 
                                                             time = div_events$time[1]))
  
  # if (scale_pops) {
  #   pops_rel_ne <- rescale(data$pop_size, 0.5, 1.2)
  #   #root_rel_ne <- rescale(div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[1], 0.5, 1.2)
  #   root_rel_ne <- rescale(div_events$new_deme_size[div_events$time == max(div_events$time)]*data$pop_size[div_events$sink + 1], 0.5, 1.2)
  #   div_events$source_rel_ne <- pops_rel_ne[div_events$source + 1]
  #   div_events$sink_rel_ne <- pops_rel_ne[div_events$sink + 1]
  #   
  #   mig_events$source_rel_ne <- pops_rel_ne[mig_events$source + 1]
  #   mig_events$sink_rel_ne <- pops_rel_ne[mig_events$sink + 1]
  # }
  
  if (scale_pops == TRUE) {
    root_pop_size <- div_events$new_deme_size[div_events$time == max(div_events$time)]*pop_size_list$historical[div_events$sink + 1]
    root_rel_ne <- rescale_alt(root_pop_size, 0.1, 1.2, c(unlist(pop_size_list), root_pop_size))
    
    pops_rel_ne <- rescale_alt(pop_size_list$current, 0.1, 1.2, pop_size_vec = c(unlist(pop_size_list), root_pop_size) )
    hist_pops_rel_ne <- rescale_alt(pop_size_list$historical, 0.1, 1.2, pop_size_vec = c(unlist(pop_size_list), root_pop_size) )
    
    mig_events$source_rel_ne <- pops_rel_ne[mig_events$source + 1]
    mig_events$sink_rel_ne <- pops_rel_ne[mig_events$sink + 1]
  }
  
  width <- 0.75
  pop_poly_list <- lapply(1:2, function(X, pops_rel_ne, hist_pops_rel_ne, width, div_events) {
    
    data.frame(y = rep(c(0, div_events$time[1]), each  = 2),
              # x = c(div_events[[switch(X, `1` = 'source', `2` = 'sink')]][1] + 1 - width*pops_rel_ne[X]/2,
              #       div_events[[switch(X, `1` = 'source', `2` = 'sink')]][1] + 1 + width*pops_rel_ne[X]/2,
              #       div_events[[switch(X, `1` = 'source', `2` = 'sink')]][1] + 1 + width*hist_pops_rel_ne[X]/2,
              #       div_events[[switch(X, `1` = 'source', `2` = 'sink')]][1] + 1 - width*hist_pops_rel_ne[X]/2)
              x = c((X - 1) + 1 - width*pops_rel_ne[X]/2,
                    (X - 1) + 1 + width*pops_rel_ne[X]/2,
                    (X - 1) + 1 + width*hist_pops_rel_ne[X]/2,
                    (X - 1) + 1 - width*hist_pops_rel_ne[X]/2)
               )
    
  }, pops_rel_ne = pops_rel_ne, hist_pops_rel_ne = hist_pops_rel_ne, width = width, div_events = div_events)
  
  #return(pop_poly_list)
  #3/30/2022
  #pop_poly_list <- pop_poly_list[match(div_events[c('source', 'sink')], 0:1)]
  # pop1_pos <- which(div_events[c('source', 'sink')] == 0)
  # pop2_pos <- which(div_events[c('source', 'sink')] == 1)
  
  
  ############################
  # first, set up the skeleton using divergences
  # it'd be good if we can make this more flexible (i.e. works with any number of div events)
  ############################
  tgap <- max(events$time)*0.10
  troot <- max(events$time)*1.33
  
  p <- ggplot() +
    xlim(0,length(data$pop_size)*1.5)+
    ylim(0,troot)
  
  #width <- 0.75
  p1 <- p + 
    geom_rect(mapping=aes(xmin=min(pop_poly_list[[1]][pop_poly_list[[1]]$y != 0,]$x),
                         xmax=max(pop_poly_list[[2]][pop_poly_list[[1]]$y != 0,]$x),
                         ymin=div_events$time[1],ymax=div_events$time[1]+tgap), fill=pop_color) +
    # geom_rect(mapping=aes(xmin=min(pop_poly_list[[pop1_pos]][pop_poly_list[[pop1_pos]]$y != 0,]$x),
    #                       xmax=max(pop_poly_list[[pop2_pos]][pop_poly_list[[pop2_pos]]$y != 0,]$x),
    #                       ymin=div_events$time[1],ymax=div_events$time[1]+tgap), fill=pop_color) + 
    # root
    geom_rect(mapping=aes(xmin=(div_events$source[1]+div_events$sink[1]+2)/2 - width*root_rel_ne/2,
                          xmax=(div_events$source[1]+div_events$sink[1]+2)/2 + width*root_rel_ne/2,
                          ymin=div_events$time[1],
                          ymax=troot), fill=pop_color) + 
    geom_polygon(data = pop_poly_list[[1]],
                 aes(x = x, y = y), fill=pop_color)  +
    geom_polygon(data = pop_poly_list[[2]],
                 aes(x = x, y = y), fill=pop_color)
   # geom_polygon(data = pop_poly_list[[pop1_pos]],
   #              aes(x = x, y = y), fill=pop_color)  + 
   # geom_polygon(data = pop_poly_list[[pop2_pos]],
   #              aes(x = x, y = y), fill=pop_color) 
  
  ########################################## 
  # add colors for populations, if desired #
  ##########################################
  if (colors) {
    
    p2 <- p1 +
      geom_polygon(data = rbind(pop_poly_list[[1]][1:2,],
                                data.frame(y = rep(tgap/2, 2),
                                           x = c(id_xval(y_val = tgap/2,
                                                         pointa = rev(unlist(pop_poly_list[[1]][2,])),
                                                         pointb = rev(unlist(pop_poly_list[[1]][3,])) ),
                                                 id_xval(y_val = tgap/2,
                                                         pointa = rev(unlist(pop_poly_list[[1]][1,])),
                                                         pointb = rev(unlist(pop_poly_list[[1]][4,]))) )  )),
                   aes(x = x, y = y), fill = pop_colors[1]) +
      geom_polygon(data = rbind(pop_poly_list[[2]][1:2,],
                            data.frame(y = rep(tgap/2, 2),
                                       x = c(id_xval(y_val = tgap/2,
                                                     pointa = rev(unlist(pop_poly_list[[2]][2,])),
                                                     pointb = rev(unlist(pop_poly_list[[2]][3,])) ),
                                             id_xval(y_val = tgap/2,
                                                     pointa = rev(unlist(pop_poly_list[[2]][1,])),
                                                     pointb = rev(unlist(pop_poly_list[[2]][4,]))) )  )),
               aes(x = x, y = y ), fill = pop_colors[2] ) #+
  #scale_fill_manual(values=pop.cols)
    
    # p2 <- p1 +
    #   geom_polygon(data = rbind(pop_poly_list[[pop1_pos]][1:2,],
    #                             data.frame(y = rep(tgap/2, 2),
    #                                        x = c(id_xval(y_val = tgap/2, 
    #                                                      pointa = rev(unlist(pop_poly_list[[pop1_pos]][2,])), 
    #                                                      pointb = rev(unlist(pop_poly_list[[pop1_pos]][3,])) ),
    #                                              id_xval(y_val = tgap/2, 
    #                                                      pointa = rev(unlist(pop_poly_list[[pop1_pos]][1,])), 
    #                                                      pointb = rev(unlist(pop_poly_list[[pop1_pos]][4,]))) )  )),
    #                aes(x = x, y = y), fill = pop_colors[1]) +
    #   geom_polygon(data = rbind(pop_poly_list[[pop2_pos]][1:2,],
    #                             data.frame(y = rep(tgap/2, 2),
    #                                        x = c(id_xval(y_val = tgap/2, 
    #                                                      pointa = rev(unlist(pop_poly_list[[pop2_pos]][2,])), 
    #                                                      pointb = rev(unlist(pop_poly_list[[pop2_pos]][3,])) ),
    #                                              id_xval(y_val = tgap/2, 
    #                                                      pointa = rev(unlist(pop_poly_list[[pop2_pos]][1,])), 
    #                                                      pointb = rev(unlist(pop_poly_list[[pop2_pos]][4,]))) )  )),
    #                aes(x = x, y = y ), fill = pop_colors[2] ) 

  } else {
    p2 <- p1
  }

  
  events_sorted <- events %>%
    arrange((time))
  
  #MIGRATION MATRICES
  time_interval_list <- time_interval(mig_mat_list = data$mig_mat_list,
                                      historical_event = events_sorted,
                                      diverge_cushion = tgap)
  
  
  migmat_arrow_info_list <- lapply(split(seq_len(length(data$mig_mat_list)), seq_len(length(data$mig_mat_list))), function(i, mig_mat, time_int_list) {
    
    if (all(!is.na(time_int_list[[i]]$end))) {
      
      migmat_arrow_info_alt(mig_mat = mig_mat[[i]],
                            start_y = time_int_list[[i]]$start,
                            #end_y = div_events$time[1],
                            end_y = time_int_list[[i]]$end,
                            offset = 0.05,
                            arrow_replicates = 1,
                            max_val = NULL,
                            max_arrow_size = 1,
                            min_arrow_size = 0.2,
                            max_arrow_head_size = max_arrow_head_size,
                            min_arrow_head_size = min_arrow_head_size,
                            pop0_x = list(pointa = rev(unlist(pop_poly_list[[1]][2,])), pointb = rev(unlist(pop_poly_list[[1]][3,]))),
                            pop1_x = list(pointa = rev(unlist(pop_poly_list[[2]][1,])), pointb = rev(unlist(pop_poly_list[[2]][4,]))),
                            direction = perspective)
      
      # migmat_arrow_info_alt(mig_mat = mig_mat[[i]],
      #                       start_y = time_int_list[[i]]$start,
      #                       #end_y = div_events$time[1],
      #                       end_y = time_int_list[[i]]$end,
      #                       offset = 0.05,
      #                       arrow_replicates = 1,
      #                       max_val = NULL,
      #                       max_arrow_size = 1,
      #                       min_arrow_size = 0.2,
      #                       max_arrow_head_size = max_arrow_head_size,
      #                       min_arrow_head_size = min_arrow_head_size,
      #                       pop0_x = list(pointa = rev(unlist(pop_poly_list[[pop1_pos]][2,])), pointb = rev(unlist(pop_poly_list[[pop1_pos]][3,]))),
      #                       pop1_x = list(pointa = rev(unlist(pop_poly_list[[pop2_pos]][1,])), pointb = rev(unlist(pop_poly_list[[pop2_pos]][4,]))),
      #                       direction = perspective)
    }
    
  }, time_int_list = time_interval_list, mig_mat = data$mig_mat_list)
  
  migmat_arrow_info_list_filtered <- migmat_arrow_info_list[sapply(migmat_arrow_info_list, function(x) !is.null(x) & !all(x == 0))]
  
  
  if (generic_plotting) {
    migmat_arrow_info_list_filtered <- lapply(migmat_arrow_info_list_filtered, function(x, arrow_head_size) {
      x$arrow_size <- 1
      x$arrow_head_size <- arrow_head_size
      return(x)
    }, arrow_head_size = generic_arrow_head_size)
  }
  
  
  mig_mat_eval <- sapply(data$mig_mat_list, function(x) !all(x == 0) & x != 'NONE')
  if (any(mig_mat_eval)) {
    #migration_period <- time_interval_list[sapply(data$mig_mat_list, function(x) !all(x == 0) & x != 'NONE')][[1]]
    migration_period <- time_interval_list[sapply(data$mig_mat_list, function(x) !all(x == 0) & any(x != 'NONE'))][[1]]
    
    p3 <- p2 +
      geom_polygon(data = data.frame(y = rep(unlist(migration_period[1, c('start', 'end')]), 2),
                                     x = c(id_xval(y_val = unlist(migration_period[1, c('start', 'end')]),
                                                   pointa = rev(unlist(pop_poly_list[[1]][2,])),
                                                   pointb = rev(unlist(pop_poly_list[[1]][3,])) ),
                                           id_xval(y_val = unlist(migration_period[1, c('start', 'end')]),
                                                   pointa = rev(unlist(pop_poly_list[[2]][1,])),
                                                   pointb = rev(unlist(pop_poly_list[[2]][4,])) )   ))[c(1, 2, 4, 3),],
                   aes(x = x, y = y), color = NA, fill = "#bfbfbf", alpha=0.5)
    
    # p3 <- p2 + 
    #   geom_polygon(data = data.frame(y = rep(unlist(migration_period[1, c('start', 'end')]), 2),
    #                                  x = c(id_xval(y_val = unlist(migration_period[1, c('start', 'end')]), 
    #                                                pointa = rev(unlist(pop_poly_list[[pop1_pos]][2,])), 
    #                                                pointb = rev(unlist(pop_poly_list[[pop1_pos]][3,])) ),
    #                                        id_xval(y_val = unlist(migration_period[1, c('start', 'end')]), 
    #                                                pointa = rev(unlist(pop_poly_list[[pop2_pos]][1,])), 
    #                                                pointb = rev(unlist(pop_poly_list[[pop2_pos]][4,])) )   ))[c(1, 2, 4, 3),],
    #                aes(x = x, y = y), color = NA, fill = "#bfbfbf", alpha=0.5)
    
  } else {
    p3 <- p2
  }
  
  #migration_period <- time_interval_list[sapply(data$mig_mat_list, function(x) !all(x == 0))][[1]]
  
  
  # p3 <- p2 + 
  #   geom_rect(data = migration_period,
  #             aes(xmin = pop0_x_right, 
  #                 xmax = pop1_x_left, 
  #                 ymin = start, ymax = end), 
  #             color = NA, fill = "#bfbfbf", alpha=0.5)
  
  
  if (length(migmat_arrow_info_list_filtered) != 0) {
    for (i in seq_len(length(migmat_arrow_info_list_filtered))) {
      
      p3 <- p3 +
        geom_segment(aes(x = migmat_arrow_info_list_filtered[[i]]$arrow_start, 
                         xend = migmat_arrow_info_list_filtered[[i]]$arrow_end, 
                         y = migmat_arrow_info_list_filtered[[i]]$arrow_pos, # might not want to hard code this? or maybe it's fine?
                         yend = migmat_arrow_info_list_filtered[[i]]$arrow_pos), # might not want to hard code this? or maybe it's fine?
                     lineend = "round", # See available arrow types in example above
                     linejoin = "round",
                     size = migmat_arrow_info_list_filtered[[i]]$arrow_size, 
                     arrow = arrow(length = unit(migmat_arrow_info_list_filtered[[i]]$arrow_head_size, "inches"), type = "closed"),
                     lty=1,
                     colour = "black" # Also accepts "red", "blue' etc
        )
    }
  }
  
  
  if (generic_plotting) {
    return(p3 + theme_void() + theme(legend.position = 'none'))
  } else {
    
    return(
      list(plot = p3 +
             theme_demovis() +
             ylab('Time (generations)') +
             theme(axis.title.y = element_text(size = 16)) +
             scale_x_continuous(expand = c(0.03, 0.03)) +
             scale_y_continuous(labels = function(x) format(x, scientific = TRUE), expand = c(0, 0)) +
             geom_segment(data=data$hist_event_info$processed_hist_events,
                          aes(x = min(do.call(rbind, pop_poly_list)$x), xend = max(do.call(rbind, pop_poly_list)$x),
                              y=time, yend=time),
                          lty=2, size = 1.25, lineend = "round") +
             geom_text(aes(x = max(do.call(rbind, pop_poly_list)$x)*(1 + label_nudge),
                           y = data$hist_event_info$processed_hist_events$time, label = formatC(round(data$hist_event_info$processed_hist_events$time), format = "d", big.mark = ",")), size = label_size, hjust = 0) +
             coord_cartesian(clip = 'off'),
           right_side = max(do.call(rbind, pop_poly_list)$x))
      
      #xlim(pop0_x_left, pop1_x_right*1.1)
    )
    
  }
  
}

