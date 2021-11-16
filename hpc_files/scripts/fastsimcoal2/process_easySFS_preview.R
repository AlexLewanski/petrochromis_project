#!/usr/bin/env Rscript

# Usage: process_easySFS_preview.R -f input_file_name.txt -o output_file_name.txt -n NUMBER
# Objective: based on the output of the easySFS --preview , identifies the projection(s) that lead to (1/2)
# the higher number of segregating sites for each population (2/2)

#load required libraries
library(optparse) #--> functionality to pass command line arguments to script
library(stringr) #--> some regular expression functionality

#create command line arguments to pass to script
parse_list <- list(make_option(c("-f", "--infile"), type = "character", default = NULL, 
                               help = "file of projection information produced from easySFS", metavar = "character"),
                   make_option(c("-o", "--output"), type = "character", default = NULL, 
                               help = "name of output file", metavar = "character"),
                   make_option(c("-n", "--number"), type = "double", default = 1, 
                               help = "number of projections to be retained per population starting with the projection resulting in the most segregating sites and decreasing sequentially based on number segregating sites)", metavar = "character"))



######################################
### CHECKING AND PROCESSING INPUTS ###
######################################

#handle command line arguments
option1 <- OptionParser(option_list = parse_list)
option2 <- parse_args(option1)

if (any(sapply(list(option2$infile, option2$number), is.null))) {
  stop("Make sure that you provide all necssary arguments (input file is the only obligatory argument)\nUsage: process_easySFS_preview.R -f 'input_file_name.txt' -o 'output_file_name.txt' -n NUMBER\nRequired arguments: input file (-f)\nOptional arguments: output file name (-o; default: easySFS_processed_projection_info.txt), number (-n; default: 1)")
}


#load file
infile <- readLines(option2$infile)

if (is.null(option2$output)) option2$output <- 'easySFS_processed_projection_info.txt'



####################################################
### PROCESSSING THE easy_SFS PREVIEW INFORMATION ###
####################################################

#extract population names --> updated to handle verbose output (which doesn't put population names on the first line) (2/21/2021)
#pop_extract <- str_extract_all(infile[1], "\\'([:alnum:]*)\\'")[[1]]
pop_extract <- str_extract_all(infile[which(str_detect(infile, 'Processing [0-9]+ population[s]* - odict_keys'))], "\\'([:alnum:]*)\\'")[[1]]
pop_extract <- gsub("'", "", pop_extract)


#process input file into a dataframe with the following columns:
#1. population name (pop)
#2. number of individuals in the projection (num_indiv)
#3. number of segregating at a particular projection value (num_seg_sites)

#order dataframe from most number segregating sites to least and subset down to the number of
#specified rows
project_list <- list()

for (NAME in pop_extract) {
  
  pop_projection_info <- infile[which(infile == NAME) + 1]
  
  individual_number <- str_extract_all(pop_projection_info, "\\([0-9]+\\,")[[1]]
  individual_number <- as.numeric(gsub("[\\(\\,]*", "", individual_number))
  
  segregating_sites <- str_extract_all(pop_projection_info, "\\, [0-9]+\\.*[0-9]+\\)")[[1]]
  segregating_sites <- as.numeric(gsub("[[:blank:]]*[\\)\\,]*", "", segregating_sites))
  
  project_df <- data.frame(pop = NAME,
                           num_indiv = individual_number,
                           num_seg_sites = segregating_sites)
  
  final_number <- ifelse(nrow(project_df) >= option2$number, option2$number, nrow(project_df))
  
  project_list[[NAME]] <- project_df[order(project_df$num_seg_sites, decreasing = TRUE),][1:final_number,]
  project_list[[NAME]]$num_seg_rank <- 1:final_number
}

#export file
write.table(do.call(rbind, project_list), option2$output, quote = FALSE, row.names = FALSE, col.names = TRUE)
