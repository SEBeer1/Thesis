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
half_sample_one_iteration = function(i, path.to.simFolder, path.to.outfolder, save.name = "onlyCSIteration"){
  
  phenoData = read.csv(paste0(path.to.simFolder, "baselineSNPeffects/data/data.snp.", i, ".csv"))
  SNP.matrix = as.matrix(read.csv(paste0(path.to.simFolder, "SNPMatrices/SNP.matrix.", i, ".csv"), row.names = 1))
  
  # include semifactor: half sample size
  phenoData = half_sample_size(phenoData)
  # remove individuals from SNP matrix
  sv.col.names = colnames(SNP.matrix)
  SNP.matrix = matrix(SNP.matrix[row.names(SNP.matrix) %in% unique(phenoData$ID)], nrow = length(unique(phenoData$ID)))
  row.names(SNP.matrix) = unique(phenoData$ID)
  colnames(SNP.matrix) = sv.col.names
  phenoData = phenoData[,c(1:2,4:6)]
  
  # centre time
  colnames(phenoData) = c("ID", "Sex", "Time", "Height", "Weight")
  phenoData$Time = phenoData$Time - mean(phenoData$Time, na.rm=T)
  
  # fit LMM without SNPs
  mod1 = baseline_mod1(phenoData)
  
  lmm.out = lmm_function1(phenoData, SNP.matrix, colnames(SNP.matrix))
  gallop.out = gallop_function1(phenoData, SNP.matrix, colnames(SNP.matrix), mod1 = mod1)
  scebe.out = scebe_function1(phenoData, SNP.matrix, colnames(SNP.matrix), mod1 = mod1)
  
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

half_sample_size = function(data){
  
  data.length = data %>% dplyr::select(ID) %>% unique() %>% pull() %>% length()
  unique.names = data %>% dplyr::select(ID) %>% unique() %>% pull()
  keep.list = sample(unique.names, round(data.length/2))
  data = data %>% filter(ID %in% keep.list)
  
  return(data)
}

pt.simFolder = "/depot/socgen/data/solomon_beer/uploadFolder/"
pt.outFolder = "/depot/socgen/data/solomon_beer/uploadFolder/sf3halfSampleSNPeffects/"
sn = "halfSampleIteration"

st = Sys.time()

for(i in (10*starting.i-9):(10*starting.i)){
  print(i)
  half_sample_one_iteration(i, pt.simFolder, pt.outFolder, save.name = sn)
}

et = Sys.time()
print(et - st)
