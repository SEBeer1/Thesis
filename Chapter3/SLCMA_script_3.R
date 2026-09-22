## script to run on purdue cluster

# packages
library(dplyr)
library(tidyr)
library(slcma)
library(gdsfmt)

# functions

# tryCatch for post selection inference
sel_inf = function(myslcma){
  tryCatch(
    #this is the chunk of code we want to run
    {slcmaInfer(myslcma, 1)
      #when it throws an error, the following block catches the error
    }, error = function(msg){
      return(NA)
    })
}

#tryCatch for efficient storage
si.tc = function(si){
  tryCatch(
    {temp = c(names(si$fli$vars), si$fli$coef0, si$fli$ci[1,1:2], si$fli$pv)
    return(temp)
    }, error = function(msg){
      return(c(NA, NA, NA, NA, NA))
    })
}

# tryCatch for max-t
max_t = function(myslcma){
  tryCatch(
    #this is the chunk of code we want to run
    {slcmaInfer(myslcma, 1, method = "maxt", do.maxtCI=FALSE)
      #when it throws an error, the following block catches the error
    }, error = function(msg){
      return(NA)
    })
}

#tryCatch for efficient storage
maxt.tc = function(maxt){
  tryCatch(
    {temp = c(row.names(maxt$maxt$coefficients), maxt$maxt$coefficients[1,1], NA, NA, maxt$maxt$coefficients[1,2])
    return(temp)
    }, error = function(msg){
      return(c(NA, NA, NA, NA, NA))
    })
}

# DME SLCMA & post selection inference
DME_SLCMA = function(data, cpg.name){
  
  data$r01 = lm(mat.smok.preg ~ 1, data=data)$resid
  data$r02 = lm(mat.age.birth ~ 1, data=data)$resid
  data$r03 = lm(parity ~ 1, data=data)$resid
  data$r04 = lm(birthweight.kg ~ 1, data=data)$resid
  data$r1 = lm(edindep_8m ~ mat.smok.preg + mat.age.birth
               + parity + birthweight.kg, data=data)$resid # and all C0s would be included in all these
  data$r2 = lm(edindep_21m ~ mat.smok.preg + mat.age.birth
               + parity + birthweight.kg + edindep_8m + sexphysabuse_18m, data=data)$resid
  data$r3 = lm(edindep_33m ~ mat.smok.preg + mat.age.birth
               + parity + birthweight.kg + edindep_8m + edindep_21m
               + sexphysabuse_18m + sexphysabuse_30m, data=data)$resid
  data$r4 = lm(edindep_61m ~ mat.smok.preg + mat.age.birth
               + parity + birthweight.kg + edindep_8m + edindep_21m + edindep_33m
               + sexphysabuse_18m + sexphysabuse_30m + sexphysabuse_42m + 
                 sexphysabuse_57m, data=data)$resid
  data$r5 = lm(edindep_6y ~ mat.smok.preg + mat.age.birth
               + parity + birthweight.kg + edindep_8m + edindep_21m + edindep_33m
               + edindep_61m + sexphysabuse_18m + sexphysabuse_30m + sexphysabuse_42m + 
                 sexphysabuse_57m + sexphysabuse_69m, data=data)$resid
  
  data$tx1 = lm(sexphysabuse_18m ~ r01 + r02 + r03 + r04 + r1 + r2
                + r3 + r4 + r5, data=data)$resid
  data$tx2 = lm(sexphysabuse_30m ~ r01 + r02 + r03 + r04 + r1 + r2
                + r3 + r4 + r5, data=data)$resid
  data$tx3 = lm(sexphysabuse_42m ~ r01 + r02 + r03 + r04 + r1 + r2
                + r3 + r4 + r5, data=data)$resid
  data$tx4 = lm(sexphysabuse_57m ~ r01 + r02 + r03 + r04 + r1 + r2
                + r3 + r4 + r5, data=data)$resid
  data$tx5 = lm(sexphysabuse_69m ~ r01 + r02 + r03 + r04 + r1 + r2
                + r3 + r4 + r5, data=data)$resid
  data$txA = lm(I(sexphysabuse_18m + sexphysabuse_30m + sexphysabuse_42m
                  + sexphysabuse_57m + sexphysabuse_69m) ~ r01 + r02
                + r03 + r04 + r1 + r2 + r3 + r4 + r5, data=data)$resid
  
  te1 = lm(Y ~ tx1, data=data)
  te2 = lm(Y ~ tx2, data=data)
  te3 = lm(Y ~ tx3, data=data)
  te4 = lm(Y ~ tx4, data=data)
  te5 = lm(Y ~ tx5, data=data)
  teA = lm(Y ~ I(tx1 + tx2 + tx3 + tx4 + tx5), data=data)
  
  # return models
  out = list(x1 = te1,
             x2 = te2,
             x3 = te3,
             x4 = te4,
             x5 = te5,
             A = teA)
  
  selection = c(te1$coef[2]*sd(data$tx1),
                te2$coef[2]*sd(data$tx2),
                te3$coef[2]*sd(data$tx3),
                te4$coef[2]*sd(data$tx4),
                te5$coef[2]*sd(data$tx5),
                teA$coef[2]*sd(data$txA))
  
  selection.index = which.max(selection^2)
  
  # print(summary(out[[selection.index]]))
  
  # obtain p-value for selection from SLCMA functions
  myslcma = slcma(Y ~ tx1 + tx2 + tx3 + tx4 + tx5
                  + Accumulation(tx1,tx2,tx3,tx4,tx5), data=data, silent = TRUE)
  
  si = sel_inf(myslcma)
  # selection, coef, cil, ciu, pval
  si.out = si.tc(si)
  
  maxt = max_t(myslcma)
  # selection, coef, sil, ciu, pval
  maxt.out = maxt.tc(maxt)
  
  dme_slcma.ci = confint(out[[selection.index]])[2,1:2] 
  dme_slcma.coefs = summary(out[[selection.index]])$coefficients[2,c(1,4)]
  
  return(list(cpg = cpg.name,
              dme_slcma = c(dme_slcma.coefs[1], dme_slcma.ci, dme_slcma.coefs[2]),
              selection = selection,
              selective.inference = si.out,
              maxt = maxt.out))
}

# SLCMA
SLCMA = function(data, cpg.name){
  
  
  # obtain p-value for selection from SLCMA functions
  myslcma = slcma(Y ~ sexphysabuse_18m + sexphysabuse_30m + sexphysabuse_42m + sexphysabuse_57m + sexphysabuse_69m +
                  + Accumulation(sexphysabuse_18m + sexphysabuse_30m + sexphysabuse_42m + sexphysabuse_57m + sexphysabuse_69m) + 
                    mat.smok.preg + mat.age.birth + parity + birthweight.kg + edindep_8m,
                  data=data, silent = TRUE, adjust = c(7:11))
  
  si = sel_inf(myslcma)
  # selection, coef, cil, ciu, pval
  si.out = si.tc(si)
  
  maxt = max_t(myslcma)
  # selection, coef, sil, ciu, pval
  maxt.out = maxt.tc(maxt)
  
  return(list(cpg = cpg.name,
              selective.inference = si.out,
              maxt = maxt.out))
}

# finds order of samples as they appear in the DNAm matrix
get_sample_indices = function(f, sample.names){
  
  all.samples = read.gdsn(index.gdsn(f, "col.names"))
  indices = which(all.samples %in% sample.names)
  
  return(indices)
}
# returns T/F for if DNAm row is in exposure data
get_sample_selection = function(f, sample.names){
  
  all.samples = read.gdsn(index.gdsn(f, "col.names"))
  selections = all.samples %in% sample.names
  
  return(selections)
}

# applies DME SLCMA to selected individuals in a column
apply_dme_slcma_DNAm_column = function(index, betas, input.data, cpg.names, sample.names){
  
  # get cpg name
  cpg.name = cpg.names[index]
  # cpg.name = row.names(betas)[index]

  # create data.frame <------------ assumes betas in same order as input.data
  to.merge = data.frame(Sample_Name = sample.names,
                        Y = betas)
  input.data = input.data %>% left_join(to.merge, by = join_by(Sample_Name))
  
  #print(to.merge)
  #print(input.data)
  
  #input.data$Y = betas
  #input.data = input.data %>% drop_na()
  
  #print(nrow(input.data))
  input.data = input.data %>% drop_na()
  #print(nrow(input.data))
  
  # apply DME SLCMA
  result = DME_SLCMA(input.data, cpg.name)
  
  # return result list
  return(result)
  
}

# applies DME SLCMA to selected individuals in a column
apply_slcma_DNAm_column = function(index, betas, input.data, cpg.names, sample.names){
  
  # get cpg name
  cpg.name = cpg.names[index]
  # cpg.name = row.names(betas)[index]
  
  # create data.frame <------------ assumes betas in same order as input.data
  to.merge = data.frame(Sample_Name = sample.names,
                        Y = betas)
  input.data = input.data %>% left_join(to.merge, by = join_by(Sample_Name))
  
  #print(to.merge)
  #print(input.data)
  
  #input.data$Y = betas
  #input.data = input.data %>% drop_na()
  
  #print(nrow(input.data))
  input.data = input.data %>% drop_na()
  #print(nrow(input.data))
  
  # apply DME SLCMA
  result = SLCMA(input.data, cpg.name)
  
  # return result list
  return(result)
  
}

################################################################################

# file paths
data.path = "/depot/socgen/data/solomon_beer/epigeneticDMEslcma/sexphy_exp.csv"
gds.path = "/depot/socgen/data/ALSPAC_2024/dnam_epic450_g0_g1/sftp/data/betas/450.gds"
save.name = "/depot/socgen/data/solomon_beer/epigeneticDMEslcma/SLCMA_sexphy_out_3.rds"

# read in exposure and covariate dataframe
input.data = read.csv(data.path)

# sample names
sample.names = input.data$Sample_Name

# open .GDS file
f = openfn.gds(gds.path)

# find order of samples in DNAm matrix
indices = get_sample_indices(f, sample.names)
# get T/F for DNAm samples in input.data
selections = get_sample_selection(f, sample.names)

cpg.names = read.gdsn(index.gdsn(f, "row.names"))
# create selection list to apply DME SLCMA to the correct DNAm measures for each individual
n.cpg = length(cpg.names)
selection.list = list(c(rep(FALSE, 400000), rep(TRUE, n.cpg - 400000)), selections)
# selection.list = list(rep(TRUE, n.cpg), selections)

# arrange input data to same order as samples appear in DNAm
input.data = cbind(input.data, indices) %>% arrange(-indices) %>% dplyr::select(-indices)

# check sample names
samples = read.gdsn(index.gdsn(f, "col.names"))
samples = samples[selections]

#print(indices[1:5])
#print(selections[12:15])
#print(input.data$Sample_Name[1:5])
#print(samples[1:5])
#print(length(samples))
#print(sum(samples %in% input.data$Sample_Name))

# apply DME SLCMA: 
appl.list <- apply.gdsn(index.gdsn(f, "matrix"), margin=1,
                  var.index = "absolute", selection = selection.list,
                  FUN=function(index, x)
                    apply_slcma_DNAm_column(index, x,
                                                input.data = input.data, cpg.names = cpg.names, sample.names = samples))

# close .GDS file
closefn.gds(f)

# save output
saveRDS(appl.list, file = save.name)
