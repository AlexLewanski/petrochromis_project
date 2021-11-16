#########################################
### SUMMARY OF THE SEQUENCING EFFORTS ###
#########################################

#OUTSTANDING QUESTIONS/THINGS TO DO:
#1. INCLUDE THE OUTROUPS IN THE SUMMARIES?
#2. COMPLETE THE READ DEPTH SUMMARY FOR THE FINAL FILTERED VCF FILE

### Reading in data ###
#metadata
cichlid_metadata <- read.csv("/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/working_metadata_file/monster_23jul20.csv")
metadata_kazumbe_polyodon_green_diagramma <- cichlid_metadata[ (cichlid_metadata$sciname2 %in% c("Petrochromis kazumbe", "Petrochromis polyodon", "Petrochromis green", "Simochromis diagramma")) & cichlid_metadata$region == 'Kigoma',]

file_path <- "/Users/alexlewanski/Documents/University_of_Wyoming/Research/Projects/Petrochromis_Project/Cichlid_FULLDATASET/sequencing_summary/"

#read information (all cichlids)
read_info <- read.table(paste0(file_path, 'assembled_per_ind_PRACTICE.txt'), header = TRUE)

#sequencing information (NOTE: THIS IS A PLACEHOLDER FOR THE FINAL DATASET)
depth_info <- read.delim(paste0(file_path, 'depth_info.txt'))



### 1. Reads summary ###
read_info_kazumbe_polyodon_green_diagramma <- read_info[read_info$ind %in% metadata_kazumbe_polyodon_green_diagramma$sample_code,]

#1a. number of reads that were matched to individuals
raw_reads <- sum(read_info_kazumbe_polyodon_green_diagramma$raw)

#at the individual level
mean(read_info_kazumbe_polyodon_green_diagramma$raw)
sd(read_info_kazumbe_polyodon_green_diagramma$raw)

#1a. number of reads that were successfully assembled to the reference genome
assembled_reads <- sum(read_info_kazumbe_polyodon_green_diagramma$assembled)

#proportion of the raw reads that were assembled to the reference genome
prop_assembled <- sum(read_info_kazumbe_polyodon_green_diagramma$assembled)/sum(read_info_kazumbe_polyodon_green_diagramma$raw)

#put info in dataframe
reads_summary_df <- data.frame(total_raw_reads = raw_reads,
                               total_assembled_reads = assembled_reads,
                               prop_assembled_reads = prop_assembled)



### 2. Average sequencing depth ###
#extracted from the final vcf file using vcftools: vcftools --vcf vcf_file --extract-FORMAT-info DP
#from the VCF manual (https://samtools.github.io/hts-specs/VCFv4.2.pdf): DP : read depth at this position for this sample (Integer)
#each row is a site and each column is the information for an individual sample
#if you take the mean of each column, this is equivalent to vcftools function: vcftools --vcf vcf_file --depth --out output_file_path



#the first two columns are chromosome and site position

#mean read depth per individual
mean_depth <- apply(depth_info[, -c(1,2)], MARGIN = 2, mean)

#mean read depth averaged across individuals
average_mean_rd <- mean(mean_depth)

#standard deviation of mean read depth
sd_mean_rd <- sd(mean_depth)


rd_summary_df <- data.frame(average_mean_rd = average_mean_rd,
                            sd_mean_rd = sd_mean_rd)



#SUMMARY DATAFRAMES
reads_summary_df
rd_summary_df

