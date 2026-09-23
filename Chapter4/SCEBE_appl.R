################################################################################
# functions
################################################################################

library(gdsfmt)
library(dplyr)
library(tidyr)
library(MASS)
library(nlme)
library(zoo)

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
# returns T/F for if cpg had a hit
get_cpg_selection = function(f, cpg.names){
  
  all.cpgs = read.gdsn(index.gdsn(f, "row.names"))
  selections = all.cpgs %in% cpg.names
  
  return(selections)
}

# scebe function
scebe_function1 = function(cpg.name, SNP.matrix, mod1){
  
  # print times
  time1 = Sys.time()
  
  phenoData = mod1$data[,c(11,1,8)] |> arrange(ID, Age.years)
  #print(phenoData)
  colnames(phenoData) = c("Sample_Name", "ID", "Time")
  # phenoData$Time = phenoData$Time - 1
  SNP.matrix = SNP.matrix[unique(phenoData$Sample_Name),] # 100
  
  #checks
  #print(cpg.name)
  print(SNP.matrix[1:5,1:10])
  print(phenoData$Sample_Name[1:10])
  #print(phenoData[1:5,])
  
  N = length(unique(phenoData$ID)) # nrow(SNP.matrix)
  noCov = ncol(SNP.matrix)
  n_measure = phenoData %>% dplyr::group_by(ID) %>% dplyr::summarize(n=length(ID)) %>% dplyr::arrange(ID)%>%dplyr::select(n)%>%unlist()%>%unname()
  IDs = phenoData$ID
  
  flme0<-summary(mod1)
  # rf=ranef(mod1)$ID # remove if using nlme
  rf=random.effects(mod1) # remove if using lme4
  #r.name=rownames(rf)
  #rf1=data.frame(ID=as.numeric(r.name),rf) %>% arrange(ID)
  #rf = rf1[, -which(names(rf1)=="ID")]
  # ind.parm <- merge(phenoData[!duplicated(phenoData$ID),c("ID",paste0("cov",1:noCov))],rf1,by="ID")
  # names(ind.parm)[c(2+noCov,3+noCov)]=c("Inter","Slope")
  names(rf)[c(1:2)]=c("Inter","Slope")
  
  # y1: rf$Inter,
  # y2: rf$Slope
  # X: as.matrix(genoData)
  
  y1 = rf$Inter
  y2 = rf$Slope
  X = as.matrix(SNP.matrix)
  
  n1 <- length(y1)
  n<-n2<-n1
  
  y1bar <- mean(y1)
  y2bar <- mean(y2)
  
  y1c <- y1 - y1bar
  y2c <- y2 - y2bar
  s1 <- colSums(X)
  s2 <- colSums(X ^ 2)
  sxx <- (s2 - (s1 ^ 2) / n)
  
  b1 <- as.vector(colSums(y1c*X) / sxx)
  b2 <- as.vector(colSums(y2c*X) / sxx)
  y1_sum <- sum(y1)
  y2_sum <- sum(y2)
  y1X <- colSums(y1 * X)
  y2X <- colSums(y2 * X)
  
  #Residuals1 <- (y1 - X %*% diag(b1) - t(replicate(n1, a1))) ## still potential for improvement
  #b1_sd <- sqrt((colSums(Residuals1 ^ 2) / (n1 - 2)) / sxx)
  #Residuals2 <- (y2 - X %*% diag(b2) - t(replicate(n2, a2))) ## still potential for improvement
  #b2_sd <- sqrt((colSums(Residuals2 ^ 2) / (n2 - 2)) / sxx)
  
  sigma11 = sum(y1 ^ 2) - (y1_sum ^ 2 * s2 - 2 * s1 * y1_sum * y1X +
                             n * (y1X) ^ 2) / (n * sxx)
  var.gamma11 = sigma11 / (sxx * (n - 2))
  
  #sigma12 = sum(y1 * y2) - (y1_sum * y2_sum * s2 - y2_sum * s1 * y1X -
  #                           y1_sum * s1 * y2X + n * y1X * y2X) / (n * sxx)
  #var.gamma12 = sigma12 / (sxx * (n - 2))
  
  sigma22 = sum(y2 ^ 2) - (y2_sum ^ 2 * s2 - 2 * s1 * y2_sum * y2X +
                             n * (y2X) ^ 2) / (n * sxx)
  var.gamma22 = sigma22 / (sxx * (n - 2))
  
  b1_sd=sqrt(var.gamma11)
  b2_sd=sqrt(var.gamma22)
  
  
  t1_value <- b1 / b1_sd
  p1_value <- pt(q = -abs(t1_value), df =  n1 - 2) * 2 # Pr(>|t|)
  
  t2_value <- b2 / b2_sd
  p2_value <- pt(q = -abs(t2_value), df =  n2 - 2) * 2 # Pr(>|t|)
  
  
  result <- cbind(b1, b1_sd, p1_value,b2, b2_sd,  p2_value)
  colnames(result) <- c("b1_Est", "b1_Std",  "Pr1(>|t|)","b2_Est", "b2_Std",  "Pr2(>|t|)")
  rownames(result) <- colnames(X)
  
  # spr = result
  
  Data=phenoData
  # R=matrix(as.numeric(VarCorr(mod1)$ID),ncol=2) # lmer
  # R = diag(as.numeric(nlme::VarCorr(mod1)[1:2,1])) # nlme uncorrelated
  R = matrix(as.numeric(nlme::getVarCov(mod1)), ncol=2) # nlme correlated
  Ip=diag(1,2)
  
  #allS=allW=list()
  #allM=allISMIS=list()
  alls11=alls12=alls21=alls22=NULL
  allw11=allw12=allw21=allw22=NULL
  allsms11=allsms12=allsms21=allsms22=NULL
  
  Time.lst = split(Data[["Time"]], Data[["ID"]])
  
  for(i in 1:N){
    nt = n_measure[i]
    time <- Time.lst[[i]] #Data[Data$ID ==i, Time]##
    Zi <- cbind(1, time)
    Gi = diag(sigma(mod1) ^ 2, nt)
    ZinvGZ <- crossprod(Zi, solve(Gi, Zi)) ##when n measure==1, return error!
    Mi=R+ginv(ZinvGZ)
    
    SIGMAi=Zi%*%R%*%t(Zi)+Gi
    Wi <- crossprod(Zi, solve(SIGMAi, Zi))
    #allW[[i]]=Wi
    allw11=c(allw11,Wi[1,1])
    allw12=c(allw12,Wi[1,2])
    allw21=c(allw21,Wi[2,1])
    allw22=c(allw22,Wi[2,2])
    
    
    
    invSi <- R %*% ZinvGZ + Ip
    Si=solve(invSi)
    alls11<-c(alls11,Si[1,1])
    alls12<-c(alls12,Si[1,2])
    alls21<-c(alls21,Si[2,1])
    alls22<-c(alls22,Si[2,2])
    
    Ipi=diag(1,2)
    SMS=(Ipi-Si)%*%Mi%*%t(Ipi-Si)
    
    allsms11<-c(allsms11,SMS[1,1])
    allsms12<-c(allsms12,SMS[1,2])
    allsms21<-c(allsms21,SMS[2,1])
    allsms22<-c(allsms22,SMS[2,2])
    
  }
  
  C<-matrix(c(sum(allw11),sum(allw12),sum(allw21),sum(allw22)),2,byrow=T)
  
  inv_C<-solve(C)
  
  invc11<-inv_C[1,1]
  invc12<-inv_C[1,2]
  invc21<-inv_C[2,1]
  invc22<-inv_C[2,2]
  
  ebe<-result[,c(1,4)] # spr[,c(1,4)]
  
  #  myX<-SNP.matrix #unique(phenoData[,5:(noCov+4)])
  Xc<-apply(SNP.matrix,2,function(x){scale(x,scale=F)}) # apply(myX,2,function(x){scale(x,scale=F)})
  
  XcX<-Xc*SNP.matrix  # Xc*myX
  Xc2<-Xc^2
  Sxc2<-colSums(Xc2)
  
  xstar11<-colSums(allw11*Xc)
  xstar12<-colSums(allw12*Xc)
  xstar21<-colSums(allw21*Xc)
  xstar22<-colSums(allw22*Xc)
  
  A11<-colSums(Xc2*alls11)/Sxc2
  A12<-colSums(Xc2*alls12)/Sxc2
  A21<-colSums(Xc2*alls21)/Sxc2
  A22<-colSums(Xc2*alls22)/Sxc2
  
  B11<-colSums(Xc*alls11)/Sxc2
  B12<-colSums(Xc*alls12)/Sxc2
  B21<-colSums(Xc*alls21)/Sxc2
  B22<-colSums(Xc*alls22)/Sxc2
  
  Sxc22<-(Sxc2)^2
  
  D11=colSums(Xc2*allsms11)/Sxc22
  D12=colSums(Xc2*allsms12)/Sxc22
  D21=colSums(Xc2*allsms21)/Sxc22
  D22=colSums(Xc2*allsms22)/Sxc22
  
  ssh211<-B11*invc11*xstar11+B12*invc21*xstar11+B11*invc12*xstar21+B12*invc22*xstar21
  ssh212<-B11*invc11*xstar12+B12*invc21*xstar12+B11*invc12*xstar22+B12*invc22*xstar22
  ssh221<-B21*invc11*xstar11+B22*invc21*xstar11+B21*invc12*xstar21+B22*invc22*xstar21
  ssh222<-B21*invc11*xstar12+B22*invc21*xstar12+B21*invc12*xstar22+B22*invc22*xstar22
  
  ssh11=1-A11+ssh211
  ssh12=0-A12+ssh212
  ssh21=0-A21+ssh221
  ssh22=1-A22+ssh222
  
  detssh=ssh11*ssh22-ssh12*ssh21
  
  invssh11<-ssh22/detssh
  invssh12<-(-ssh12/detssh)
  invssh21<-(-ssh21/detssh)
  invssh22<-ssh11/detssh
  
  Vgamma211<-B11*invc11*B11+B12*invc21*B11+B11*invc12*B12+B12*invc22*B12
  #Vgamma211<-Vgamma211/SXc2
  Vgamma212<-B11*invc11*B21+B12*invc21*B21+B11*invc12*B22+B12*invc22*B22
  #Vgamma212<-Vgamma212/SXc2
  Vgamma221<-B21*invc11*B11+B22*invc21*B11+B21*invc12*B12+B22*invc22*B12
  #Vgamma221<-Vgamma221/SXc2
  Vgamma222<-B21*invc11*B21+B22*invc21*B21+B21*invc12*B22+B22*invc22*B22
  #Vgamma222<-Vgamma222/SXc2
  
  Vgamma11=D11-Vgamma211
  Vgamma12=D12-Vgamma212
  Vgamma21=D21-Vgamma221
  Vgamma22=D22-Vgamma222
  
  Vssh11<-invssh11*Vgamma11*invssh11+invssh12*Vgamma21*invssh11+invssh11*Vgamma12*invssh12+invssh12*Vgamma22*invssh12
  Vssh12<-invssh11*Vgamma11*invssh21+invssh12*Vgamma21*invssh21+invssh11*Vgamma12*invssh22+invssh12*Vgamma22*invssh22
  Vssh21<-invssh21*Vgamma11*invssh11+invssh22*Vgamma21*invssh11+invssh21*Vgamma12*invssh12+invssh22*Vgamma22*invssh12
  Vssh22<-invssh21*Vgamma11*invssh21+invssh22*Vgamma21*invssh21+invssh21*Vgamma12*invssh22+invssh22*Vgamma22*invssh22
  
  est.ssh.inter<-invssh11*ebe[,1]+invssh12*ebe[,2]
  est.ssh.slope<-invssh21*ebe[,1]+invssh22*ebe[,2]
  
  t.ssh.inter<-est.ssh.inter/sqrt(Vssh11)
  t.ssh.slope<-est.ssh.slope/sqrt(Vssh22)
  
  p.ssh.inter<-2*pnorm(abs(t.ssh.inter),lower.tail=F)
  p.ssh.slope<-2*pnorm(abs(t.ssh.slope),lower.tail=F)
  
  output<-cbind(est.ssh.inter,sqrt(Vssh11),p.ssh.inter,t.ssh.inter,est.ssh.slope,sqrt(Vssh22),p.ssh.slope,t.ssh.slope)
  
  colnames(output)= c("est.sc.b1", "se.sc.b1", "p.sc.b1", "t.sc.b1",
                      "est.sc.b2", "se.sc.b2", "p.sc.b2", "t.sc.b2") 
    
  time2 = Sys.time()
  print(paste0(cpg.name, ":   ", time2 - time1))
  
  return(output)
  
}

################################################################################
# file paths
## longitudinal outcome and covariate data
data.path = "/depot/socgen/data/solomon_beer/epiSCEBE/YC.df.csv"
## beta value DNAm predictors
gds.path = "/depot/socgen/data/ALSPAC_2024/dnam_epic450_g0_g1/sftp/data/betas/450.gds"
## mod1 path
mod1.path = "/depot/socgen/data/solomon_beer/epiSCEBE/mod1.rds"
## CpG keep list path
cpg.path = "/depot/socgen/data/solomon_beer/epiSCEBE/CpG.keep.list.csv"
## save name
save.name = "/depot/socgen/data/solomon_beer/epiSCEBE/epiSCEBE.rds"

# read in exposure and covariate dataframe
input.data = read.csv(data.path)
# read in CpG omitted model
mod1 = readRDS("mod1.rds")

# read in CpG keep list
cpgs = read.csv(cpg.path)$x

# sample names
sample.names = input.data$Sample_Name

# open .GDS file
f = openfn.gds(gds.path)

# find order of samples in DNAm matrix
indices = get_sample_indices(f, sample.names)
# get T/F for DNAm samples in input.data
selections = get_sample_selection(f, sample.names)
# get T/F for CpG getting a hit
cpg.selections = get_cpg_selection(f, cpgs)

cpg.names = read.gdsn(index.gdsn(f, "row.names"))
cpg.names = cpg.names[cpg.selections]
# create selection list to apply SCEBE to the correct DNAm measures for each individual
n.cpg = length(cpg.names)
# selection.list = list(c(rep(TRUE, 100), rep(FALSE, n.cpg-100)), selections)
# selection.list = list(rep(TRUE, n.cpg), selections)
selection.list = list(cpg.selections, selections)

# arrange input data to same order as samples appear in DNAm
# input.data = cbind(input.data, indices) %>% arrange(-indices) %>% dplyr::select(-indices)

# check sample names
samples = read.gdsn(index.gdsn(f, "col.names"))
samples = samples[selections]

#print(indices[1:5])
#print(selections[12:15])
#print(input.data$Sample_Name[1:5])
#print(samples[1:5])
#print(length(samples))
#print(sum(samples %in% input.data$Sample_Name))

# read in DNAm as SNP.matrix
CpG.matrix = readex.gdsn(index.gdsn(f, "matrix"), selection.list)
colnames(CpG.matrix) = samples
rownames(CpG.matrix) = cpg.names #[1:100]
CpG.matrix = t(CpG.matrix)
CpG.matrix = 100*CpG.matrix
CpG.matrix = na.aggregate(CpG.matrix, FUN = median)
print(nrow(CpG.matrix))
print(ncol(CpG.matrix))
print(CpG.matrix[1:5,1:5])

# close .GDS file
closefn.gds(f)

# obtain summaries of all CpG sites
# sum.cpg.df = data.frame(ID = rownames(CpG.matrix), CpG.matrix)
# print(sum.cpg.df[1:5,1:10])
# sum.cpg.df = sum.cpg.df %>% pivot_longer(!ID, names_to = "CpG", values_to = "beta") %>% group_by(CpG) %>% summarise(n = n(), n.not.na = sum(!is.na(beta)),mean = mean(beta,na.rm=T), median = median(beta, na.rm=T),sd = sd(beta, na.rm=T), min = min(beta, na.rm=T),max = max(beta, na.rm=T))
# print(sum.cpg.df[1:5,])
# write.csv(sum.cpg.df, "cpg.summaries.csv")

# apply SCEBE: 
output = scebe_function1("weight.adj.height", CpG.matrix, mod1)
output = data.frame(output) %>% drop_na()

# save output
# saveRDS(appl.list, file = save.name)
write.csv(output, "output.csv", row.names=T)
