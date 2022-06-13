#!/usr/bin/env Rscript

### parse command line arguments ###
args <- commandArgs(trailingOnly = TRUE)
ind <- args[1]
output <- args[2]
keyword <- args[3]


#load file
ind_file <- read.delim(ind, header = FALSE)

#randomize sample codes of taxa identified by the input keyword
ind_file[grepl(pattern = keyword, ind_file$V3),]$V3 <- sample(ind_file[grepl(pattern = keyword, ind_file$V3),]$V3)

#export randomized ind file
write.table(ind_file, file = output,
            row.names = FALSE, col.names = FALSE, quote = FALSE)
