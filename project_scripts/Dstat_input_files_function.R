######################################################################
######################################################################
### CREATING INPUT FILES FOR D STATISTIC COMPUATIONS IN ADMIXTOOLS ###
######################################################################
######################################################################


### INFORMATION ON dstat_inputfiles FUNCTION ###
#The dstat_inputfiles function takes a dataframe containing information on the population and species membership of samples and returns a list of files needed to compute D-stats using AdmixTools to compare differential sharing of alleles between sympatric populations of two species and a disjunct population of one of the species. 
#Specifically, following the form (((P1, P2), P3), Outgroup), the function will create files to compute D-stats with the following structures: (((SP1_popY, SP1_popX), SP2_popX), Outgroup) and (((SP2_popY, SP2_popX), SP1_popX), Outgroup). 
#The function can take a dataframe with other species that will not be included in the D-stats and the user specifies which taxa that they want to include. Currently, dstat_inputfiles only handles two species.
#OUTPUT: a list containing 3 elements: 
    #(1) information on the populations included in D-stat input files including which populations and the sample sizes for each species in each population (not used as an AdmixTools input)
    #(2) popfile that contains the different combinations of species/populations used to in D-stat calculations (AdmixTools input)
    #(3) indfile that contains the species/population identifier for each sample (AdmixTools input)


### INFORMATION ON USING THIS FUNCTION IN ANOTHER SCRIPT ###
#To use this function, run the following function at the top of your script:
#source('/path/to/script/Dstat_input_files_function.R')


### RESOURCES ###
#https://speciationgenomics.github.io/ADMIXTOOLS_admixr/
#https://github.com/DReichLab/AdmixTools/blob/master/README.Dstatistics


dstat_inputfiles <- function(pop_info,
                             pop_column, 
                             species_column, 
                             sampleID_column, 
                             species1, 
                             species2, 
                             outgroup_samples, 
                             minimum_sample_size = 2,
                             type = c('a', 'b')) {
  
  #========================================
  #== explanation of function arguments === 
  #========================================
  
  #pop_info: dataframe that contains information on the populations and species membership of individuals (needs to have a column or each)
  #pop_column: name of the column with population membership information
  #species_column: name of the column with species identity information
  #sampleID_column: name of column that contains the sample IDs for each sample
  #species1: first species to include in D-stat calculations
  #species2: second species to include in D-stat calculations
  #outgroup_samples: the samples to be used as the outgroup for D-stat calculations
  #minimum_sample_size (default = 2): the minimum number of individuals for each species that need to be found at a site for the site to be included in D-stats
  #type (either 'a' or 'b')"
  #a creates inputs for D-stats with the following topologies: ((sp1_disjunct, sp1_sympatric), sp2_sympatric) and ((sp2_disjunct, sp2_sympatric), sp1_sympatric)
  #b creates inputs for D-stats with the following topology: ((sp1_disjunct, sp1_sympatric), sp2_fullpop)
  
  
  
  #==============
  #== updates === 
  #==============
  
  #UPDATE 10/17/2020:
  #Previously, I was ordering the pop file in the following way: outgroup, p3, p2, p1. However, I reversed the order to p1, p2, p3, outgroup, which follows
  #the approach described in https://speciationgenomics.github.io/ADMIXTOOLS_admixr/. My previous implementation method calculates the D-stat numerator as 
  #(out - 3)(p2 - p1) and this more recent change calculates the D-stat numerator as (p1 - p2)(p3 - out). The changing of the order results in the same 
  #numerator value (and it doesn't affect the denominator), and thus this altered pop file order results in identical D-stat values.
  
  #AdmixTools calculates D-stats in the following way:
  #With a tree ((w, x), y), z, the D-statistic is calculated as 
  #D = sum(numerator)/sum(denominator),
  #where numerator = (w - x)(y - z) and denominator = (w + x - 2wx)(y + z - 2yz)
  
  #With D > 0 and |z| > 3, w and y show evidence of excess of allele sharing
  #With D < 0 and |z| > 3, x and y show evidence of excess of allele sharing
  
  #TO RECOVER THE PREVIOUS ORDERING, UNCOMMENT THE LAST LINE IN THE FOLLOWING SECTION: create the popfile for D-Stat computation in AdmixTools
  
  
  
  #===============================
  #== check function inputs of === 
  #===============================

  type <- match.arg(type, several.ok = FALSE)
  #if ( !(type %in% c("a", "b")) )
  #  stop("type only accepts 'a' or 'b'")
  
  if( FALSE %in% (outgroup_samples %in% pop_info[, sampleID_column] == TRUE) )
    stop("all outgroup samples need to be included in the pop_info dataframe")
  
  if( FALSE %in% ( c(pop_column, species_column, sampleID_column) %in% colnames(pop_info) ) )
    stop("pop_info dataframe does not contain all of the specific columns")
  
  
  
  #=================================================================================
  #== create the dataframe of species and populations to be included in pop file ===
  #=================================================================================
  
  #subset pop_info dataframe down to focal species and remove levels of species that are empty (if the column is a factor)
  if (is.factor(pop_info[, species_column])) {
    focal_taxa_info <- droplevels(pop_info[pop_info[, species_column] %in% c(species1, species2),])
  } else {
    focal_taxa_info <- pop_info[pop_info[, species_column] %in% c(species1, species2),]
  }
  
  #create contingency table of the frequency of the focal species across the different populations (convert to data frame)
  species_freq_table <- as.data.frame.matrix(table(focal_taxa_info[,pop_column], focal_taxa_info[,species_column]))
  species_freq_table$pop <- rownames(species_freq_table) #add a column of for population in the dataframe
  
  #subset the frequency species x population frequency dataframe to the populations where both species have the minimum sample size
  pops_to_include <- species_freq_table[which(species_freq_table[, species1] >= minimum_sample_size & species_freq_table[,species2] >= minimum_sample_size),]
  
  #outgroup info
  outgroup_info <- pop_info[pop_info[,sampleID_column] %in% outgroup_samples,]
  
  
  
  #==============================================================
  #== create the popfile for D-Stat computation in AdmixTools ===
  #==============================================================
  
  #empty dataframe to create the popfile
  pop_dataframe <- data.frame(p1 = character(), p2 = character(), p3 = character())
  
  if (type == "a") {
    #type a
    
    for (i in 1:nrow(pops_to_include)) {
      focal_pop <- pops_to_include[,'pop'][i]
      pop_dataframe <- rbind(pop_dataframe, data.frame(p1 = paste0(species1, '_', pops_to_include[,'pop'][-i]), 
                                                       p2 = paste0(species1, '_', focal_pop), 
                                                       p3 = paste0(species2, '_', focal_pop)))
      
      pop_dataframe <- rbind(pop_dataframe, data.frame(p1 = paste0(species2, '_', pops_to_include[,'pop'][-i]), 
                                                       p2 = paste0(species2, '_', focal_pop), 
                                                       p3 = paste0(species1, '_', focal_pop)))
    }
    
  } else {
    #type b
    
    for (i in 1:nrow(pops_to_include)) {
      focal_pop <- pops_to_include[,'pop'][i]
      pop_dataframe <- rbind(pop_dataframe, data.frame(p1 = paste0(species1, '_', pops_to_include[,'pop'][-i]), 
                                                       p2 = paste0(species1, '_', focal_pop), 
                                                       p3 = species2))
    }
  }
  
  #new pop_popdataframe that will contain info on the outgroup(s)
  pop_dataframe_final <- data.frame(p1 = character(), p2 = character(), p3 = character())
  
  #loop through each species included as an outgroup
  for (sp in unique(outgroup_info[, species_column])) {
    pop_dataframe$outgroup <- sp
    pop_dataframe_final <- rbind(pop_dataframe_final, pop_dataframe)
  }
  
  
  #PREVIOUS POP FILE ORDER:
  #create final popfile dataframe by reordering columns in the following order: outgroup, p3, p2, p1
  #pop_dataframe_final <- pop_dataframe_final[,c(4, 3, 2, 1)] #UNCOMMENT THIS LINE TO RECOVER THE PREVIOUS (i.e. PRE-10/17/2020 UPDATE) ORDERING OF TAXA
  
  
  
  #==============================================================
  #== create the indfile for D-Stat computation in AdmixTools ===
  #==============================================================
  
  if (type == 'a') {
    #type a
    indfile <- data.frame(individual = focal_taxa_info[, sampleID_column], 
                          u = "U", 
                          taxa = paste(focal_taxa_info[, species_column], focal_taxa_info[, pop_column], sep = "_"))
    
  } else {
    #type b
    indfile <- data.frame(individual = focal_taxa_info[, sampleID_column], 
                          u = "U", 
                          taxa = ifelse(focal_taxa_info[, species_column] == species1, 
                                        paste(focal_taxa_info[, species_column], focal_taxa_info[, pop_column], sep = "_"),
                                        as.character(focal_taxa_info[, species_column])) )
  }
  
  indfile <- rbind(indfile, data.frame(individual = outgroup_info[, sampleID_column], 
                                       u = "U", 
                                       taxa = outgroup_info[, species_column]) )
  
  
  
  #=====================================================================================================
  #== output a list that contains: (1) information on included populations, (2) popfile, (3) indfile ===
  #=====================================================================================================
  
  #output is a list that contains the info on each species' sample sizes per location, the popfile, and the indfile
  return(list(pop_info = pops_to_include, 
              popfile = pop_dataframe_final, 
              indfile = indfile))
}



######################
### EXAMPLE OF USE ###
######################

#on a mac to comment/uncomment a set of lines, highlight the lines and use the following command: command + shift + c

### Creating example dataframe ###
# example_dataframe1 <- data.frame(pop = rep('pop1', 5), species = rep('sp1', 5), sampleID_column = paste0('samp_', 1:5))
# example_dataframe2 <- data.frame(pop = rep('pop1', 5), species = rep('sp2', 5), sampleID_column = paste0('samp_', 6:10))
# example_dataframe3 <- data.frame(pop = rep('pop2', 5), species = rep('sp1', 5), sampleID_column = paste0('samp_', 11:15))
# example_dataframe4 <- data.frame(pop = rep('pop2', 5), species = rep('sp2', 5), sampleID_column = paste0('samp_', 16:20))
# example_dataframe5 <- data.frame(pop = rep('pop3', 5), species = rep('sp1', 5), sampleID_column = paste0('samp_', 21:25))
# example_dataframe6 <- data.frame(pop = rep('pop3', 2), species = rep('sp2', 2), sampleID_column = paste0('samp_', 26:27))
# example_dataframe7 <- data.frame(pop = rep('pop4', 10), species = rep('sp2', 10), sampleID_column = paste0('samp_', 30:39))
# example_dataframe8 <- data.frame(pop = rep('pop4', 10), species = rep('sp4', 10), sampleID_column = paste0('samp_', 130:139))
# example_dataframe9 <- data.frame(pop = rep('pop6', 25), species = rep('sp1', 25), sampleID_column = paste0('samp_', 221:225))
# example_dataframe10 <- data.frame(pop = rep('pop6', 5), species = rep('sp2', 5), sampleID_column = paste0('samp_', 226:230))
# 
# example_dataframe_full <- rbind(example_dataframe1, example_dataframe2, example_dataframe3, example_dataframe4, example_dataframe5, example_dataframe6, example_dataframe7, example_dataframe8, example_dataframe9, example_dataframe10)
# 
# ###
# dfoil_input_function_example <- dstat_inputfiles(pop_info = example_dataframe_full,
#                                                 pop_column = 'pop',
#                                                 species_column = 'species',
#                                                 sampleID_column = 'sampleID_column',
#                                                 species1 = 'sp1',
#                                                 species2 = 'sp2',
#                                                 outgroup_samples = c('samp_132', 'samp_133'),
#                                                 minimum_sample_size = 2,
#                                                 type = "a")

