#!/usr/bin/env Rscript

### PARSING INPUTS ###
args <- commandArgs(trailingOnly = TRUE)
sfs_file <- args[1] #the name of the sfs file that is read into the script
output_name <- args[2] #the name of the output file
fsc_sfs <- readLines(sfs_file) #reading in sfs file



### PROCESSING INPUTS ###
fsc_sfs_strip <- as.numeric(strsplit(fsc_sfs[3], split = " ")[[1]]) #split the sfs line into a vector

#create the snp sfs by removing the monomorphic sitres entry (1st) and the everything past the halfway point entry (these are all empty)
#NOTE: this is currently built for folded sfs --> future step could be to make script so it handles folded and unfolded sfs
snp_sfs <- fsc_sfs_strip[2:(((length(fsc_sfs_strip) - 1)/2) + 1)] #create the snp sfs by removing the monomorphic sites entry (1st) and all entries after the (nseq/2 + 1)th entry (the entries corresponding to major alleles)

total_sites <- round(sum(fsc_sfs_strip)) #total sites the sum of all entries in the sfs including monomorphic sites
nseq <- length(fsc_sfs_strip) - 1 #number of individuals (nseq for stairway plot blueprint)
nrand_vec <- round(c((nseq - 2)/4, (nseq - 2)/2, (nseq - 2)*(3/4), nseq - 2)) #creating the vector of random breakpoint values based on suggestion in stairway plot manual



### OUTPUTTING RESULTS ###
write(paste0("processed_sfs: ", paste(snp_sfs, collapse = " ")), file = output_name, append = FALSE)
write(paste0("number_sites: ", total_sites), file = output_name, append = TRUE)
write(paste0("number_individuals: ",  nseq), file = output_name, append = TRUE)
write(paste0("nrand: ",  paste(nrand_vec, collapse = " ")), file = output_name, append = TRUE)
