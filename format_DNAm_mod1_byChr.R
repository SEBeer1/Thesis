# ------------------------------------------------------------------------------
# Decompose Chr mod1 output into input for SCEBE
# ------------------------------------------------------------------------------

# packages
suppressMessages(library(argparse))
suppressMessages(library(miceadds))

# get chr from argument
parser <- ArgumentParser()
parser$add_argument("--chr", type="character")
Args <- parser$parse_args()
chr <- Args$chr

# load in .Rdata
load.Rdata(paste0("/data2/sbeer/lmQTL/mod1andresid/Chr", chr, "Output.Rdata"), "list.data")

# find first non singular cpg - hope there is one in first 1000!
first.nonsingular.i = 0 # first non singular
for(i in 1:1000){
  if(list.data[[i]]$singular == FALSE){
    first.nonsingular.i = i
    break
  }
}

# find num samples & sample names
rnames = list.data[[first.nonsingular.i]][[2]]$data$Sample_Name
num.betas = length(rnames)

# iteration through list of lists
# i  num cpgs
# j  num non singular cpgs
num.cpg = length(list.data)
j = 1
cpg.nonsing = rep(NA, num.cpg)
cpg.resid = matrix(rep(NA, num.cpg*num.betas), ncol = num.cpg)
cpg.mod1 = vector("list", num.cpg)

# loop through each cpg in list of lists
for(i in 1:num.cpg){
  
  # if singular is false:
  if(list.data[[i]]$singular == FALSE){
    
    # record cpg name
    cpg.nonsing[j] = list.data[[i]][[1]][[1]]
    
    # add column to resid matrix
    cpg.resid[,j] = as.numeric(list.data[[i]][[3]])
    
    # add model fit to list
    mod.temp = list.data[[i]][[2]]
    cpg.mod1[j] = list(mod.temp)
    
    # increase nonsing cpg count
    j = j + 1
    
    # progress check
    print(i)
    
  }
}

# remove singular cpgs
cpg.nonsing = cpg.nonsing[!is.na(cpg.nonsing)]
cpg.resid = cpg.resid[,1:length(cpg.nonsing)]
colnames(cpg.resid) = cpg.nonsing
row.names(cpg.resid) = rnames
cpg.mod1 = Filter(Negate(is.null), cpg.mod1)
names(cpg.mod1) = cpg.nonsing

print(paste0("Number of non-singular cpg sites: ", length(cpg.nonsing)))

# save output
write.csv(cpg.nonsing, file = paste0("/data2/sbeer/lmQTL/DNAm/nsDNAm/candidateLists/candidateListChr", chr, ".csv"), row.names = F)
write.csv(cpg.resid, file = paste0("/data2/sbeer/lmQTL/DNAm/nsDNAm/cpgResid/residMatrixChr", chr, ".csv"), row.names = T)
save(cpg.mod1, file = paste0("/data2/sbeer/lmQTL/mod1andresid/mod1/mod1Chr", chr, ".Rdata"))

