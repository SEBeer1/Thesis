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
MAR_one_iteration = function(i, path.to.simFolder, path.to.outfolder, save.name = "onlyCSIteration"){
  
  phenoData = read.csv(paste0(path.to.simFolder, "baselineSNPeffects/data/data.snp.", i, ".csv"))
  SNP.matrix = as.matrix(read.csv(paste0(path.to.simFolder, "SNPMatrices/SNP.matrix.", i, ".csv"), row.names = 1))
  
  # include semifactor: MAR
  LR.params = get_LR_param(phenoData)
  
  test.df = phenoData %>% mutate(Wave = rep(c(4:9), length(unique(ID)))) %>%
    mutate(Wave = ifelse(Wave == 9, 3, Wave)) %>%
    mutate(prev.weight = ifelse(Wave == 3, NA, Weight)) %>%
    dplyr::select(-c(Age.years, Weight, Height))
  colnames(test.df) = c("ID", "Sex", "Wave", "prev.weight")
  
  phenoData = left_join(phenoData %>% mutate(Wave = rep(c(3:8), length(unique(ID)))), test.df) %>%
    mutate(prob_miss = ifelse(is.na(prev.weight), 0,
                              prob_missing(prev.weight, LR.params$alpha, LR.params$beta))) %>% 
    mutate(Weight = miss_1_0(prob_miss)*Weight) %>%
    mutate(Sex = ifelse(is.na(Weight), NA, Sex)) %>%
    mutate(Age.years = ifelse(is.na(Weight), NA, Age.years)) %>%
    mutate(Height = ifelse(is.na(Weight), NA, Height)) %>%
    dplyr::select(-c(Wave, prev.weight, prob_miss)) %>%
    dplyr::select(c(ID,Sex,Age.years,Height,Weight))
  
  #phenoData = phenoData[,c(1:2,4:6)]
  
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

get_LR_param = function(data){
  mean.w1 = data %>% filter(Wave == 3) %>% dplyr::select(Weight) %>% pull() %>% mean()
  mean.w5 = data %>% filter(Wave == 7) %>% dplyr::select(Weight) %>% pull() %>% mean()
  beta = (log(0.9/0.1)/(mean.w5 - mean.w1))
  alpha = beta*mean.w5
  return(list(alpha = alpha, beta = beta))
}

# MAR functions
prob_missing = function(y_m1, alpha, beta){
  
  prob_miss = ifelse(is.na(y_m1), yes = 1, no = exp(alpha-beta*y_m1)/(1+exp(alpha-beta*y_m1)))
  
  return(1-prob_miss)
}

miss_1_0 = function(prob_miss){
  
  i=1
  t = rep(NA, length(prob_miss))
  for(i in 1:length(prob_miss)){
    t[i] = ifelse(runif(1)>prob_miss[i],1,NA)
  }
  #t = ifelse(runif(1)>prob_miss,1,NA)
  #print(t)
  return(t)
}

pt.simFolder = "/depot/socgen/data/solomon_beer/uploadFolder/"
pt.outFolder = "/depot/socgen/data/solomon_beer/uploadFolder/sf7MAR/"
sn = "MARIteration"

st = Sys.time()

for(i in (10*starting.i-9):(10*starting.i)){
  print(i)
  MAR_one_iteration(i, pt.simFolder, pt.outFolder, save.name = sn)
}

et = Sys.time()
print(et - st)
