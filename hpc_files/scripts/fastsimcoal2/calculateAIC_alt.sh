#!/usr/bin/env Rscript


# Usage: calculateAIC.sh modelprefix


# This script calculates AIC from fsc modeling results
# Run in the folder with the highest likelihood
#YOU NEED TO LOAD R (i.e. module load r) BEFORE RUNNING THIS SCRIPT

# Read model name
args=commandArgs(TRUE)

# Checks if model name was given
if(length(args)<1){
  stop("ERROR: No input / model name given\nUsage: fsc-calculateAIC.R modelname")
}


relative_lik <- function(aic_info) {
  delta_AIC <- aic_info - min(aic_info)
  exp_delta_AIC <- exp(-0.5 * delta_AIC)
  
  return(exp_delta_AIC/sum(exp_delta_AIC))
}


# Check if model.bestlhoods file exists
if(file.exists(paste(args[1],".bestlhoods",sep=""))){
  bestlhoods<-read.delim(paste(args[1],".bestlhoods",sep=""))
}else{
  stop(paste("ERROR: Aborted. No file ",args[1],".bestlhoods file exists",sep=""))
}

# Check if model.est file exists
if(file.exists(paste(args[1],".est",sep=""))){
  est<-readLines(paste(args[1],".est",sep=""))
}else{
  stop(paste("ERROR: Aborted. No file ",args[1],".est file exists in this directory!\nUsage: fsc-calculateAIC.R modelname",sep=""))
}

# Count number of parameters
#3/18/2022: new hopefully more robust and flexible way to count the number of parameters
#k<-(grep("RULES",est))-(grep("//all Ns* are",est)+1) #previous way (from Joana's original version)
est_reduced <- est[!(seq_along(est) %in% grep('^$|^//', est))] #remove empty lines and comment lines (starting with //)
start_pos <- grep('^\\[PARAMETERS\\]$', est_reduced) #this should be 1
end_pos <- min(grep('RULES|COMPLEX PARAMETERS', est_reduced)) - 1 #allows for optional presence of RULES section
k <- end_pos - start_pos #number of estimated parameters

# Calculate AIC
AIC<-2*k-2*(bestlhoods$MaxEstLhood/log10(exp(1)))

# Calculate delta-likelihood
deltaL<-bestlhoods$MaxObsLhood-bestlhoods$MaxEstLhood

# Output model.AIC file in simulation folder
write.table(cbind(k,deltaL,AIC,bestlhoods$MaxObsLhood,bestlhoods$MaxEstLhood),paste(args[1],".AIC",sep=""),row.names = F,col.names = T,sep = "\t",quote = F)
