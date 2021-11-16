#######################################################################
#######################################################################
### FUNCTIONS FOR EVALUATING FIT OF FASTSIMCOAL2 DEMOGRAPHIC MODELS ###
#######################################################################
#######################################################################


#############################
### OVERVIEW OF FUNCTIONS ###
#############################

#name: relative_lik
#purpose: calculates the relative likelihood (Akaike's weight of evidence)
#input: vector of AIC values
#required package(s): none

#name: rounding_func
#purpose: rounds values up or down. Used internally in plot_2d_sfs.
#input: vector of numbers
#required package(s): none

#name: plot_2d_sfs
#purpose: visualize the expected and observed 2D SFS including individual
#         plots of each and a comparison of the two
#input: list of 2D SFS.
#required package(s): ggplot2, reshape2

#name: plot_1d_sfs
#purpose: visualize the expected and observed marginal 1D SFSs from 2D SFS
#input: list of 2D SFS.
#required package(s): tidyverse


### INFORMATION ON USING THESE FUNCTIONS IN ANOTHER SCRIPT ###
#To use this function, run the following function at the top of your script:
#source('/path/to/script/sfs_compare_vis_funcs.R')



#################
### FUNCTIONS ###
#################

relative_lik <- function(aic) {
  delta_aic <- aic - min(aic)
  exp_delta_aic <- exp(-0.5 * delta_aic)
  
  return(exp_delta_aic/sum(exp_delta_aic))
}


rounding_func <- function(val, resolution, direction = c('up', 'down')) {
  direction <- match.arg(direction, several.ok = FALSE)
  
  if (direction == 'up') return(ceiling(val/(resolution))*(resolution))
  if (direction == 'down') return(floor(val/(resolution))*(resolution))
}


plot_2d_sfs <- function(sfs_list,
                        plot_obs = TRUE,
                        plot_exp = TRUE,
                        plot_comparison = TRUE,
                        log10 = TRUE,
                        drop_monomorph = TRUE,
                        minimum_count = 1,
                        divergent_col_scheme = c('blue', 'white', 'red'),
                        sequential_col_scheme = c('white', 'red'),
                        scale_resolution = c(obs_exp = 1, comparison = 0.5)) {
  
  ### FUNCTIONS ###
  #taken from: https://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation
  check_color <- function(x) {
    sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)), 
               error = function(e) FALSE)
    })
  }
  
  ### REQUISITE PACKAGES ###
  if (any(!require(ggplot2, reshape2))) 
    stop('Please load the two required packages: ggplot2 and reshape2.', call. = FALSE)
  
  
  ### CHECKING INPUTS ###
  #plot_obs, plot_exp, plot_comparison, log10, drop_monomorph
  if(any(!(c(plot_obs, plot_exp, plot_comparison, log10, drop_monomorph) %in% c(TRUE, FALSE))))
    stop('plot_obs, plot_exp, plot_comparison, and log10 all need to be booleans', call. = FALSE)
  
  #sfs_list
  if ( (plot_obs | plot_comparison) & !( 'observed' %in% names(sfs_list)) )
    stop('If you want to plot the observed SFS, you need to include an entry in sfs_list named "obseved".', call. = FALSE)
  if ((plot_exp | plot_comparison) & !( 'expected' %in% names(sfs_list)) )
    stop('If you want to plot the expected SFS, you need to include an entry in sfs_list named "expected".', call. = FALSE)
  if(all(sapply(sfs_list, function(x) any(x < 0) )) | !all(sapply(sfs_list, is.numeric )) | !all(sapply(sfs_list, is.matrix )) )
    stop('sfs_list needs to include only numeric matrices with exclusively non-negative entries.', call. = FALSE)
  
  #minimum_count
  if (!is.numeric(minimum_count) | length(minimum_count) != 1)
    stop('minimum_count needs to be a single, numeric entry.', call. = FALSE)
  
  #divergent_col_scheme
  if (length(divergent_col_scheme[check_color(x = divergent_col_scheme)]) != 3 )
    stop('divergent_col_scheme needs to include exactly three valid colors.', call. = FALSE)
  
  #sequential_col_scheme
  if (length(sequential_col_scheme[check_color(x = sequential_col_scheme)]) != 2 )
    stop('sequential_col_scheme needs to include exactly two valid colors.', call. = FALSE)
  
  
  ### PROCESSING INPUTS ###
  if (all(c('observed', 'expected') %in% names(sfs_list))) {
    #sfs count in observed sfs (excluding sites monomorphic across both species)
    sfscount_observed <- sum(sfs_list[['observed']]) - sfs_list[['observed']][1, 1]
    sfs_list[['expected']] <- round(sfs_list[['expected']]*sfscount_observed)
  }
  
  sfs_list <- lapply(sfs_list, function(x) { 
    colnames(x) <- 0:(ncol(x) - 1)
    #rownames(x) <- 0:(ncol(x) - 1)
    rownames(x) <- 0:(nrow(x) - 1)
    return(x)
  })
  
  sfs_list_processed <- lapply(sfs_list, function(x, minimum_count, drop_monomorph, log10) {
    x[x < minimum_count] <- 0
    if (drop_monomorph) x[1, 1] <- 0
    if (log10) x <- log10(x)
    x[is.infinite(x)] <- 0
    return(x)
  }, minimum_count = minimum_count, drop_monomorph = drop_monomorph, log10 = log10)
  
  #setting up the list of plots to create
  plot_vec <- c(observed = plot_obs, expected = plot_exp, comparison = plot_comparison)
  plot_vec_processed <- names(plot_vec[plot_vec]) #reduce the plotting list down to TRUE elements
  
  
  ### PLOTTING ###
  if (length(plot_vec_processed) > 0) {
    
    if ('comparison' %in% plot_vec_processed) { 
      sfs_list_processed[['exp_obs_prop_sfs']] <- sfs_list[['expected']]/sfs_list[['observed']]
      sfs_list_processed[['exp_obs_prop_sfs']][is.infinite(sfs_list_processed[['exp_obs_prop_sfs']]) | is.nan(sfs_list_processed[['exp_obs_prop_sfs']])] <- NA
      sfs_list_processed[['exp_obs_prop_sfs']][1, 1] <- NA
      sfs_list_processed[['exp_obs_prop_sfs']][sfs_list[['observed']] < minimum_count] <- NA
    }
    
    sfs_plot_list <- lapply(split(plot_vec_processed, factor(plot_vec_processed, levels = plot_vec_processed ) ),
                            function(sfs, sfs_list_processed, sfs_list, scale_resolution) {
                              
                              if (sfs %in% c('observed', 'expected')) {
                                
                                return(ggplot(melt(sfs_list_processed[[sfs]]), aes(x = Var1, y = Var2)) +
                                         geom_raster(aes(fill = value)) +
                                         scale_fill_gradient(low = sequential_col_scheme[1], 
                                                             high = sequential_col_scheme[2],
                                                             na.value = 'white',
                                                             #limits = c(0, round_by_half(val = max(sfs_list_processed[[sfs]], na.rm = TRUE), direction = 'up'))
                                                             limits = c(0, rounding_func(val = max(sfs_list_processed[[sfs]], na.rm = TRUE), resolution = scale_resolution[['obs_exp']], direction = 'up'))
                                         ) +
                                         theme(axis.text.x=element_text(size = 9, angle=0, vjust = 0.3),
                                               axis.text.y=element_text(size = 9),
                                               plot.title=element_text(size = 11)) +
                                         scale_x_continuous(expand = c(0, 0)) +
                                         scale_y_continuous(expand = c(0, 0)) + 
                                         theme_bw()
                                )
                              }
                              
                              if (sfs == 'comparison') {
                                melt_df <- melt(sfs_list_processed[['exp_obs_prop_sfs']])
                                return(ggplot(melt_df, aes(x = Var1, y = Var2)) +
                                         geom_raster(aes(fill = value)) +
                                         scale_fill_gradient2(low = divergent_col_scheme[1],
                                                              mid = divergent_col_scheme[2],
                                                              high = divergent_col_scheme[3],
                                                              na.value = 'white',
                                                              midpoint = 1,
                                                              space = "Lab",
                                                              limits = c(0, rounding_func(val = max(melt_df$value, na.rm = TRUE), resolution = scale_resolution[['comparison']], direction = 'up'))
                                         ) +
                                         theme(axis.text.x=element_text(size = 9, angle = 0, vjust = 0.3),
                                               axis.text.y=element_text(size = 9),
                                               plot.title=element_text(size = 11)) +
                                         scale_x_continuous(expand = c(0, 0)) +
                                         scale_y_continuous(expand = c(0, 0)) + 
                                         theme_bw()
                                )
                              }
                              
                            }, sfs_list_processed = sfs_list_processed, sfs_list = sfs_list, scale_resolution = scale_resolution)
  } else {
    sfs_plot_list <- 'No plots were created.'
  }
  
  ### OUTPUTTING RESULTS ###
  return(list(processed_sfs = sfs_list_processed,
              plots = sfs_plot_list) )
}



plot_1d_sfs <- function(sfs_list,
                        max_freq,
                        log10 = FALSE,
                        plot = TRUE) {
  
  ### REQUISITE PACKAGES ###
  if (!require(tidyverse))
    stop('Please load the tidyverse.', call. = FALSE)
  
  
  ### CHECKING INPUTS ###
  #sfs_list
  if ( !all(c('observed', 'expected') %in%  names(sfs_list)) )
    stop('The sfs list input needs to contain a matrix named "observed" and a list or matrix named "expected."')
  if ( !is.matrix(sfs_list[['observed']]) | !is.numeric(sfs_list[['observed']]) | any(sfs_list[['observed']] < 0) )
    stop('The observed entry of sfs_list needs to be a numeric matrix with exclusively non-negative entries.', call. = FALSE)
  
  check_exp_list <- ifelse(is.list(sfs_list[['expected']]), TRUE, FALSE)
  if (check_exp_list) {
    if(all(sapply(sfs_list[['expected']], function(x) any(x < 0) )) | !all(sapply(sfs_list[['expected']], is.numeric )) | !all(sapply(sfs_list[['expected']], is.matrix )) )
      stop('All of the elements of sfs_list[["expected"]] need to be numeric matrices with exclusively non-negative entries.', call. = FALSE)
  } else {
    if ( !is.matrix(sfs_list[['expected']]) | !is.numeric(sfs_list[['expected']]) | any(sfs_list[['expected']] < 0) )
      stop('The expected entry of sfs_list needs to be a numeric matrix with exclusively non-negative entries.', call. = FALSE)
  }
  
  #max_freq
  if (!is.numeric(max_freq) | length(max_freq) != 1 | max_freq < 1)
    stop('max_freq needs to be a single, numeric entry that is >= 1.', call. = FALSE)
  
  #drop_monomorph, log10, plot
  #if(any(!(c(drop_monomorph, log10, plot) %in% c(TRUE, FALSE))))
  #  stop('plot_obs, plot_exp, plot_comparison, and log10 all need to be booleans', call. = FALSE)
  if(any(!(c(log10, plot) %in% c(TRUE, FALSE))))
    stop('plot_obs, plot_exp, plot_comparison, and log10 all need to be booleans', call. = FALSE)
  
  
  ### PROCESSING SFS MATRICES ###
  #sfscount_observed <- sum(sfs_list[['observed']]) - sfs_list[['observed']][1, 1] #count of sites in the observed sfs (excluding monomorphic sites)
  sfs_list[['observed']][c(1, length(sfs_list[['observed']]))] <- 0
  sfscount_observed <- sum(sfs_list[['observed']])
  
  #obs.SFS[c(1,length(obs.SFS))] <- 0
  #exp.SFS <- round(exp.SFS*sum(obs.SFS))
  #exp.SFS[1,1] <- 0; exp.SFS[nrow(exp.SFS),ncol(exp.SFS)] <- 0
  
  
  marg_counts <- list()
  
  if (check_exp_list) {
    
    #scaling the expected SFS based on the observed SFS
    sfs_list[['expected']] <- lapply(sfs_list[['expected']], function(sfs, count) {
      sfs_processed <- round(sfs*count)
      sfs_processed[1,1] <- 0
      sfs_processed[nrow(sfs_processed), ncol(sfs_processed)] <- 0
      return(sfs_processed)
    }, count = sfscount_observed )
    
    #ensure that the specified max_freq is not greater than the minimum frequency in the expected SFSs (and adjust accordingly if it is invalid)
    min_obs_freq <- min(sapply(sfs_list[['expected']], function(x) min(dim(x) - 1)))
    max_freq <- ifelse(min_obs_freq < max_freq, min_obs_freq, max_freq)
    
    #marginalize the 2D SFS based on each taxa by summing SFS counts across rows and columns
    exp_sfs_s1 <- lapply(list(row = 1, column = 2), function(y, max_freq) {
      do.call(rbind, lapply(sfs_list[['expected']], function(x, max_freq) {
        apply(x, MARGIN = y, sum)[1:(max_freq + 1)]
      }, max_freq = max_freq))
    }, max_freq = max_freq)
    
    
    exp_sfs_s2 <- lapply(split(names(exp_sfs_s1), names(exp_sfs_s1)), function(x, sfs_info, log10) {
      
      count_col <- apply(sfs_info[[x]], 2, mean)
      max_col <- apply(sfs_info[[x]], 2, max)
      min_col <- apply(sfs_info[[x]], 2, min)
      
      if (log10) {
        count_col <- log10(count_col); count_col[is.infinite(count_col)] <- 0
        max_col <- log10(max_col); max_col[is.infinite(max_col)] <- 0
        min_col <- log10(min_col); min_col[is.infinite(min_col)] <- 0
      }
      
      return(data.frame(count = count_col,
                        freq = 0:(ncol(sfs_info[[x]]) - 1),
                        direction = x,
                        sfs_type = 'expected',
                        max = max_col,
                        min = min_col))
    }, sfs_info =  exp_sfs_s1, log10 = log10)
    
    marg_counts[['expected']] <- do.call(rbind, exp_sfs_s2)
    
  } else {
    sfs_list[['expected']] <- round(sfs_list[['expected']]*sfscount_observed)
    sfs_list[['expected']][1,1] <- 0
    sfs_list[['expected']][nrow(sfs_list[['expected']]), ncol(sfs_list[['expected']])] <- 0
    
    max_freq <- ifelse(min(dim(sfs_list[['expected']]) - 1) < max_freq, min(dim(sfs_list[['expected']]) - 1), max_freq)
    
    if (log10) {
      count_col <- log10(c(apply(sfs_list[['expected']], 1, sum)[1:(max_freq + 1)], apply(sfs_list[['expected']], 2, sum)[1:(max_freq + 1)]))
      count_col[is.infinite(count_col)] <- 0
    } else {
      count_col <- c(apply(sfs_list[['expected']], 1, sum)[1:(max_freq + 1)], apply(sfs_list[['expected']], 2, sum)[1:(max_freq + 1)])
    }
    
    marg_counts[['expected']] <- data.frame(count = count_col,
                                            freq = rep(0:max_freq, 2),
                                            direction = rep( c('row', 'column'), each = max_freq + 1),
                                            sfs_type = 'expected')
  }
  
  
  #if (drop_monomorph) sfs_list[['observed']][1, 1] <- 0
  max_freq <- ifelse(min(dim(sfs_list[['observed']]) - 1) < max_freq, min(dim(sfs_list[['observed']]) - 1), max_freq)
  
  if (log10) {
    count_col <- log10(c(apply(sfs_list[['observed']], 1, sum)[1:(max_freq + 1)], apply(sfs_list[['observed']], 2, sum)[1:(max_freq + 1)]))
    count_col[is.infinite(count_col)] <- 0
  } else {
    count_col <- c(apply(sfs_list[['observed']], 1, sum)[1:(max_freq + 1)], apply(sfs_list[['observed']], 2, sum)[1:(max_freq + 1)])
  }
  
  marg_counts[['observed']] <- data.frame(count = count_col,
                                          freq = rep(0:max_freq, 2),
                                          direction = rep( c('row', 'column'), each = max_freq + 1),
                                          sfs_type = 'observed')
  
  if (check_exp_list) {
    marg_counts[['observed']]$max <- NA
    marg_counts[['observed']]$min <- NA
  }
  
  
  marg_counts <- do.call(rbind, marg_counts) #collapse the processed data into a single dataframe
  rownames(marg_counts) <- NULL #remove row names
  
  
  ### PLOTTING PROCESSED SFS DATA ###
  if (plot) {
    plot_list <- lapply( split(unique(marg_counts$direction), unique(marg_counts$direction)), function(direc, sfs_data, check_exp_list) {
      
      if (check_exp_list) {
        
        return(ggplot(data = sfs_data[sfs_data$direction == direc,] %>% 
                        mutate(sfs_type = factor(sfs_type, levels = c('observed', 'expected'))),
                      aes(fill = sfs_type, y = count, x = freq) ) +
                 geom_bar(position = "dodge", stat = "identity") +
                 geom_errorbar(position = position_dodge(width = 0.9), aes(ymax = min, ymin = max), width = 0.3) +
                 scale_x_continuous(breaks = 0:max(sfs_data[sfs_data$direction == direc,]$freq)))
        
      } else {
        
        return(ggplot(data = sfs_data[sfs_data$direction == direc,] %>% 
                        mutate(sfs_type = factor(sfs_type, levels = c('observed', 'expected'))),
                      aes(fill = sfs_type, y = count, x = freq) ) +
                 geom_bar(position = "dodge", stat = "identity") +
                 scale_x_continuous(breaks = 0:max(sfs_data[sfs_data$direction == direc,]$freq)))
        
      }
    }, sfs_data = marg_counts, check_exp_list = check_exp_list)
    
  } else {
    plot_list <- 'No plots were created.'
  }
  
  #output the result of plotting and the processed data
  return(list(plots = plot_list,
              marg_counts = marg_counts))
}
