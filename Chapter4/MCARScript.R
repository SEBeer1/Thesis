# base script
# take arguments from command line
arg = commandArgs(trailingOnly = T)
starting.i = as.numeric(arg[1])
print(starting.i)

# upload test script
library(nlme)
library(tidyr)
library(dplyr)
library(ggplot2)
library(MASS)

# load in method functions
source("/depot/socgen/data/solomon_beer/uploadFolder/methodFunctions.R")

# function to run one iteration
MCAR_one_iteration = function(i, path.to.simFolder, path.to.outfolder, save.name = "onlyCSIteration"){
  
  phenoData = read.csv(paste0(path.to.simFolder, "baselineSNPeffects/data/data.snp.", i, ".csv"))
  SNP.matrix = as.matrix(read.csv(paste0(path.to.simFolder, "SNPMatrices/SNP.matrix.", i, ".csv"), row.names = 1))
  
  # include semifactor: increase me
  phenoData = MCAR.drop(phenoData, 0.10, 697)
  
  phenoData = phenoData[,c(1:2,4:6)]
  
  # centre time
  colnames(phenoData) = c("ID", "Sex", "Time", "Height", "Weight")
  phenoData$Time = phenoData$Time - mean(phenoData$Time, na.rm=T)
  
  # fit LMM without SNPs
  mod1 = baseline_mod1(phenoData)
  
  lmm.out = lmm_function1(drop_na(phenoData), SNP.matrix, colnames(SNP.matrix))
  gallop.out = gallop_function1_missing_data(phenoData, SNP.matrix, colnames(SNP.matrix), mod1 = mod1)
  scebe.out = scebe_function1(drop_na(phenoData), SNP.matrix, colnames(SNP.matrix), mod1 = mod1)
  
  out.df = rbind(lmm.out$output, gallop.out$output, scebe.out$output)
  out.df$iteration = rep(i, nrow(out.df))
  write.csv(out.df, paste0(path.to.outfolder, save.name, i, "out.csv"), row.names = F)
  
  time.df = time_taken_table(lmm.out, gallop.out, scebe.out)
  write.csv(time.df, paste0(path.to.outfolder, save.name, i, "time.csv"), row.names = F)
  
  return(paste0("Iteration ", i, " completed"))
}

timediff <- function(time){
  return(as.difftime(time, units = 'secs'))
}

time_taken_table = function(lmm.out, gallop.out, scebe.out){
  
  tb = data.frame(
    method = c("LMM", "GALLOP", "SCEBE"),
    total.time = as.numeric(c(timediff(lmm.out$time.taken.total), timediff(gallop.out$time.taken.total), timediff(scebe.out$time.taken.total))),
    one.off.time = c(NA, timediff(gallop.out$time.taken.oneoff), timediff(scebe.out$time.taken.oneoff)),
    SNPs.time = c(NA, timediff(gallop.out$time.taken.SNPs), timediff(scebe.out$time.taken.SNPs)),
    per.SNP.time = c(timediff(lmm.out$time.taken.total)/nrow(lmm.out$output), timediff(gallop.out$time.taken.SNPs)/nrow(gallop.out$output), timediff(scebe.out$time.taken.SNPs)/nrow(scebe.out$output))
  )
  
  return(tb)
}

# MCAR function
MCAR.drop = function(data, drop.p, n){
  
  drop.row.array = c(2 + 6*sample(1:n,n*(1*drop.p)), 3 + 6*sample(1:n,n*(2*drop.p)), 4 + 6*sample(1:n,n*(3*drop.p)),
                     5 + 6*sample(1:n,n*(4*drop.p)), 6 + 6*sample(1:n,n*(5*drop.p)))
  
  data[row.names(data) %in% drop.row.array, 2] = rep(NA, length(which(row.names(data) %in% drop.row.array))) # <--- may need to add this back in
  data[row.names(data) %in% drop.row.array, 4] = rep(NA, length(which(row.names(data) %in% drop.row.array)))
  data[row.names(data) %in% drop.row.array, 5] = rep(NA, length(which(row.names(data) %in% drop.row.array)))
  data[row.names(data) %in% drop.row.array, 6] = rep(NA, length(which(row.names(data) %in% drop.row.array)))
  data[row.names(data) %in% drop.row.array, 7] = rep(NA, length(which(row.names(data) %in% drop.row.array)))
  
  return(data)
}

pt.simFolder = "/depot/socgen/data/solomon_beer/uploadFolder/"
pt.outFolder = "/depot/socgen/data/solomon_beer/uploadFolder/sf6MCAR/"
sn = "MCARIteration"

st = Sys.time()

for(i in (10*starting.i-9):(10*starting.i)){
  print(i)
  MCAR_one_iteration(i, pt.simFolder, pt.outFolder, save.name = sn)
}

et = Sys.time()
print(et - st)
