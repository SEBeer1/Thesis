# need my version of this
library(data.table)

load("/data2/sbeer/lmQTL/DNAm/DCHS_betas_longitudinal.Rdata")

for (chr in 1:22){
  
  print(paste0("Chromosome ", chr))
  
  candidates <- fread(paste0("/data2/sbeer/lmQTL/DNAm/DNAmSubsets/CpGbyChr/cpgChr", chr, ".csv"))
  print(length(unique(candidates$CpG)))
  candidates <- candidates[candidates$CpG %in% row.names(dnam)]
  print(length(unique(candidates$CpG)))

  dnam_subset <- dnam[unique(candidates$CpG),] 
  save(dnam_subset,
       file = paste0("/data2/sbeer/lmQTL/DNAm/DCHS_betas_longitudinal_chr", chr, ".Rdata"))
  
}
