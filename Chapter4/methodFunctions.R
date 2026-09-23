########################## script for LMM methods ##############################

################################################################################
# functions to return mod1
################################################################################

# baseline covariates
# uncorrelated randome effects
baseline_mod1 = function(phenoData){
  
  colnames(phenoData) = c("ID", "Sex", "Time", "Height", "Weight")
  phenoData$Time = phenoData$Time - mean(phenoData$Time, na.rm=T)
  mod1 = lme(Weight ~ Time + Sex + Height, random = list(ID = pdDiag(~Time)),
                data = phenoData,
                na.action = "na.omit",
                control=c(msMaxIter=100, opt="optim")
)
  
  return(mod1)
  
}

# tvc
tvc_mod1 = function(phenoData){
  
  colnames(phenoData) = c("ID", "Sex", "tvc", "Time", "Height", "Weight")
  phenoData$Time = phenoData$Time - mean(phenoData$Time, na.rm=T)
  mod1 = lme(Weight ~ Sex + tvc*Time + Height, random = list(ID = pdDiag(~Time)),
                data = phenoData,
                na.action = "na.omit",
                control=c(msMaxIter=100, opt="optim")
  )
  
  return(mod1)
  
}

# non lin
quadratic_mod1 = function(phenoData){
  
  colnames(phenoData) = c("ID", "Sex", "Time", "Height", "Weight")
  phenoData$Time = phenoData$Time - mean(phenoData$Time, na.rm=T)
  mod1 = lme(Weight ~ Sex + poly(Time, 2) + Height, random = list(ID = pdDiag(~Time)), data = phenoData, na.action = "na.omit")
  
  return(mod1)
  
}

################################################################################
# gallop functions
################################################################################

# baseline
gallop_function1 = function(phenoData, SNP.matrix, SNP.names, mod1, k = 6){
  
  colnames(phenoData) = c("ID", "Sex", "Time", "Height", "Weight")
  
  phenoData$Time = phenoData$Time - mean(phenoData$Time, na.rm=T)  # centre time to "To deal with potential multicolinearity caused by adding interaction effect of SNP x time"
  phenoData = phenoData[,-1] # remove ID column
  phenoData[is.na(phenoData)] <- 0 # set NAs as 0, as per GALLOP paper
  
  ## start time here
  start.time = Sys.time()
  
  # Extract variance components and compute penalty matrix (P)
  varcor = nlme::VarCorr(mod1)
  G = diag(as.numeric(varcor[1:2,1]))
  sig = as.numeric(varcor[3,2])
  P = solve(G / sig ^ 2)
  
  # Put data in convenient arrays and vector
  TT =  cbind(1, phenoData[, "Time"]) # matrix with column of ones and time at each measurement
  nxy = ncol(phenoData)
  X = as.matrix(cbind(1, phenoData[ , -nxy])) # data matrix with no ID or outcome column
  nx = ncol(X)
  y = phenoData[ , nxy] # outcome
  
  # Compute components of block-diagonal system with covariates  (without SNP)
  # Additionally compute and store object SS = Rot %*% Si, which is used later on
  ncov = 2 # number of covariates
  n = nrow(SNP.matrix) # number of individuals
  k = k # maximum number of outcome measurements for an individual
  
  A21 = matrix(NA, 2 * n , 2 + ncov)
  q2  = rep(NA, 2 * n)
  SS = matrix(NA, 2 * n , 2)
  
  for (i in 1 : n ) {
    
    uk = (i - 1) * k + (1 : k)
    u2 = (i - 1) * 2 + (1 : 2)
    Ti = TT[uk, ]
    Si = crossprod(Ti, Ti)
    sv = svd(Si + P)
    Rot = sqrt(1 /sv$d) * sv$u
    Q = Rot %*% t(Ti)
    SS[u2, ] = Rot %*% Si
    A21[u2, ] = Q %*% X[uk,  ]
    q2[u2] = Q %*% y[uk]
  }
  
  q1 = crossprod(X, y)
  A11 = crossprod(X)
  # Solve the system (20)
  Q = A11 - crossprod(A21)
  q = q1 - crossprod(A21, q2)
  sol = solve(Q, q)
  blups = q2 - A21 %*% sol
  
  # Compute sums of products per subject involved in the crossprod(X, G), crossprod(G) and crossprod(G, y). We use row-wise Kronecker product to avoid repeating SNP vector k times
  ex = matrix(1, 1, ncov + 2)
  et = matrix(1, 1, 2)
  XTk = kronecker(et, X) * kronecker(TT, ex) ## check here for how GALLOP handles TVC
  TTk = kronecker(et, TT) * kronecker(TT, et)
  Tyk = y * TT
  XTs = matrix(0, n, ncol(XTk))
  TTs = matrix(0, n, ncol(TTk))
  Tys = matrix(0, n, 2)
  AtS = matrix(0, n, 2 * nx)
  for (i in 1:n) {
    uk = (i - 1) * k + (1 : k)
    XTs[i, ] = apply(XTk[uk, ], 2, sum)
    TTs[i, ] = apply(TTk[uk, ], 2, sum)
    Tys[i, ] = apply(Tyk[uk, ], 2, sum)
    u2 = (i - 1) * 2 + (1 : 2)
    AtS[i, ] = c(crossprod(A21[u2, ], SS[u2, ]))
  }
  
  ## SNPs start here
  end.oneoff.time = Sys.time()
  
  SNPS = SNP.matrix
  ns = ncol(SNPS)
  
  # Add SNPs one by one and solve
  Theta = D = matrix(NA, ns, 2)
  for (i in 1 : ns) {
    si = SNPS[, i]
    snp2 = rep(si, each = 2)
    H1 = matrix(crossprod(si, XTs), nx, 2)
    H2 = snp2 * SS
    AtH = matrix(crossprod(si, AtS), nx, 2)
    R = H1 - AtH
    Cfix = solve(Q, R)
    Cran = H2 - A21 %*% Cfix
    GtG = matrix(crossprod(si ^ 2, TTs), 2, 2)
    Gty = matrix(crossprod(si, Tys), 2, 1)
    V = GtG - crossprod(H1, Cfix) - crossprod(H2, Cran)
    v = Gty - crossprod(H1, sol) - crossprod(H2, blups)
    Theta[i , ] = solve(V, v)
    D[i, ] = diag(solve(V))
  }
  
  SE = sqrt(sig**2) * sqrt(D)
  Pval = 2 * pnorm(-abs(Theta / SE))
  
  ## end time here
  end.time = Sys.time()
  
  # create output dataframe
  gallop.out.df = data.frame(
    SNP = SNP.names,
    cs.est = Theta[,1],
    cs.se = SE[,1],
    cs.t = Theta[,1] / SE[,1],
    cs.pval = Pval[,1],
    lon.est = Theta[,2],
    lon.se = SE[,2],
    lon.t = Theta[,2] / SE[,2],
    lon.pval = Pval[,2],
    method = rep("GALLOP", length(SNP.names))
  )
  
  out.list = list(
    output = gallop.out.df,
    time.taken.oneoff = end.oneoff.time - start.time,
    time.taken.SNPs = end.time - end.oneoff.time,
    time.taken.total = end.time - start.time
  )
  return(out.list)
  
}

# missing data
gallop_function1_missing_data = function(phenoData, SNP.matrix, SNP.names, mod1, k = 6){
  
  colnames(phenoData) = c("ID", "Sex", "Time", "Height", "Weight")
  
  phenoData$Time = phenoData$Time - mean(phenoData$Time, na.rm=T)  # centre time to "To deal with potential multicolinearity caused by adding interaction effect of SNP x time"
  phenoData$Sex = phenoData$Sex - mean(phenoData$Sex, na.rm=T)
  phenoData$Height = phenoData$Height - mean(phenoData$Height, na.rm=T)
  phenoData = phenoData[,-1] # remove ID column
  
  ## start time here
  start.time = Sys.time()
  
  # Extract variance components and compute penalty matrix (P)
  varcor = nlme::VarCorr(mod1)
  G = diag(as.numeric(varcor[1:2,1]))
  sig = as.numeric(varcor[3,2])
  P = solve(G / sig ^ 2)
  
  ## introduced to handle missing data
  missing.array = 1 - is.na(phenoData$Time)*1
  phenoData[is.na(phenoData)] <- 0 # set NAs as 0, as per GALLOP paper
  
  # Put data in convenient arrays and vector
  TT =  cbind(missing.array, phenoData[, "Time"]) # matrix with column of ones and time at each measurement
  nxy = ncol(phenoData)
  X = as.matrix(cbind(missing.array, phenoData[ , -nxy])) # data matrix with no ID or outcome column
  nx = ncol(X)
  y = phenoData[ , nxy] # outcome
  
  # Compute components of block-diagonal system with covariates  (without SNP)
  # Additionally compute and store object SS = Rot %*% Si, which is used later on
  ncov = 2 # number of covariates
  n = nrow(SNP.matrix) # number of individuals
  k = k # maximum number of outcome measurements for an individual
  
  A21 = matrix(NA, 2 * n , 2 + ncov)
  q2  = rep(NA, 2 * n)
  SS = matrix(NA, 2 * n , 2)
  
  for (i in 1 : n ) {
    
    uk = (i - 1) * k + (1 : k)
    u2 = (i - 1) * 2 + (1 : 2)
    Ti = TT[uk, ]
    Si = crossprod(Ti, Ti)
    sv = svd(Si + P)
    Rot = sqrt(1 /sv$d) * sv$u
    Q = Rot %*% t(Ti)
    SS[u2, ] = Rot %*% Si
    A21[u2, ] = Q %*% X[uk,  ]
    q2[u2] = Q %*% y[uk]
  }
  
  q1 = crossprod(X, y)
  A11 = crossprod(X)
  # Solve the system (20)
  Q = A11 - crossprod(A21)
  q = q1 - crossprod(A21, q2)
  sol = solve(Q, q)
  blups = q2 - A21 %*% sol
  
  # Compute sums of products per subject involved in the crossprod(X, G), crossprod(G) and crossprod(G, y). We use row-wise Kronecker product to avoid repeating SNP vector k times
  ex = matrix(1, 1, ncov + 2)
  et = matrix(1, 1, 2)
  XTk = kronecker(et, X) * kronecker(TT, ex) ## check here for how GALLOP handles TVC
  TTk = kronecker(et, TT) * kronecker(TT, et)
  Tyk = y * TT
  XTs = matrix(0, n, ncol(XTk))
  TTs = matrix(0, n, ncol(TTk))
  Tys = matrix(0, n, 2)
  AtS = matrix(0, n, 2 * nx)
  for (i in 1:n) {
    uk = (i - 1) * k + (1 : k)
    XTs[i, ] = apply(XTk[uk, ], 2, sum)
    TTs[i, ] = apply(TTk[uk, ], 2, sum)
    Tys[i, ] = apply(Tyk[uk, ], 2, sum)
    u2 = (i - 1) * 2 + (1 : 2)
    AtS[i, ] = c(crossprod(A21[u2, ], SS[u2, ]))
  }
  
  ## SNPs start here
  end.oneoff.time = Sys.time()
  
  SNPS = SNP.matrix
  ns = ncol(SNPS)
  
  # Add SNPs one by one and solve
  Theta = D = matrix(NA, ns, 2)
  for (i in 1 : ns) {
    si = SNPS[, i]
    snp2 = rep(si, each = 2)
    H1 = matrix(crossprod(si, XTs), nx, 2)
    H2 = snp2 * SS
    AtH = matrix(crossprod(si, AtS), nx, 2)
    R = H1 - AtH
    Cfix = solve(Q, R)
    Cran = H2 - A21 %*% Cfix
    GtG = matrix(crossprod(si ^ 2, TTs), 2, 2)
    Gty = matrix(crossprod(si, Tys), 2, 1)
    V = GtG - crossprod(H1, Cfix) - crossprod(H2, Cran)
    v = Gty - crossprod(H1, sol) - crossprod(H2, blups)
    Theta[i , ] = solve(V, v)
    D[i, ] = diag(solve(V))
  }
  
  # sig2 is sigma squared
  
  SE = sqrt(sig**2) * sqrt(D)
  Pval = 2 * pnorm(-abs(Theta / SE))
  
  ## end time here
  end.time = Sys.time()
  
  # create output dataframe
  gallop.out.df = data.frame(
    SNP = SNP.names,
    cs.est = Theta[,1],
    cs.se = SE[,1],
    cs.t = Theta[,1] / SE[,1],
    cs.pval = Pval[,1],
    lon.est = Theta[,2],
    lon.se = SE[,2],
    lon.t = Theta[,2] / SE[,2],
    lon.pval = Pval[,2],
    method = rep("GALLOP", length(SNP.names))
  )
  
  out.list = list(
    output = gallop.out.df,
    time.taken.oneoff = end.oneoff.time - start.time,
    time.taken.SNPs = end.time - end.oneoff.time,
    time.taken.total = end.time - start.time
  )
  return(out.list)
  
}

# tvc
gallop_function2 = function(phenoData, SNP.matrix, SNP.names, mod1, k = 6){
  
  colnames(phenoData) = c("ID", "Sex", "tvc", "Time", "Height", "Weight")
  phenoData$Time = phenoData$Time - mean(phenoData$Time, na.rm=T)  # centre time to "To deal with potential multicolinearity caused by adding interaction effect of SNP x time"
  phenoData = phenoData[,-1] # remove ID column
  
  ## start time here
  start.time = Sys.time()
  
  # Extract variance components and compute penalty matrix (P)
  varcor = nlme::VarCorr(mod1)
  G = diag(as.numeric(varcor[1:2,1]))
  sig = as.numeric(varcor[3,2])
  P = solve(G / sig ^ 2)
  
  # Put data in convenient arrays and vector
  TT =  cbind(1, phenoData[, "Time"]) # matrix with column of ones and time at each measurement
  nxy = ncol(phenoData)
  X = as.matrix(cbind(1, phenoData[ , -nxy])) # data matrix with no ID or outcome column
  nx = ncol(X)
  y = phenoData[ , nxy] # outcome
  
  # Compute components of block-diagonal system with covariates  (without SNP)
  # Additionally compute and store object SS = Rot %*% Si, which is used later on
  ncov = 3 # number of covariates
  n = nrow(SNP.matrix) # number of individuals
  k = k # maximum number of outcome measurements for an individual
  
  A21 = matrix(NA, 2 * n , 2 + ncov)
  q2  = rep(NA, 2 * n)
  SS = matrix(NA, 2 * n , 2)
  
  for (i in 1 : n ) {
    
    uk = (i - 1) * k + (1 : k)
    u2 = (i - 1) * 2 + (1 : 2)
    Ti = TT[uk, ]
    Si = crossprod(Ti, Ti)
    sv = svd(Si + P)
    Rot = sqrt(1 /sv$d) * sv$u
    Q = Rot %*% t(Ti)
    SS[u2, ] = Rot %*% Si
    A21[u2, ] = Q %*% X[uk,  ]
    q2[u2] = Q %*% y[uk]
  }
  
  q1 = crossprod(X, y)
  A11 = crossprod(X)
  # Solve the system (20)
  Q = A11 - crossprod(A21)
  q = q1 - crossprod(A21, q2)
  sol = solve(Q, q)
  blups = q2 - A21 %*% sol
  
  # Compute sums of products per subject involved in the crossprod(X, G), crossprod(G) and crossprod(G, y). We use row-wise Kronecker product to avoid repeating SNP vector k times
  ex = matrix(1, 1, ncov + 2)
  et = matrix(1, 1, 2)
  XTk = kronecker(et, X) * kronecker(TT, ex) ## check here for how GALLOP handles TVC
  TTk = kronecker(et, TT) * kronecker(TT, et)
  Tyk = y * TT
  XTs = matrix(0, n, ncol(XTk))
  TTs = matrix(0, n, ncol(TTk))
  Tys = matrix(0, n, 2)
  AtS = matrix(0, n, 2 * nx)
  for (i in 1:n) {
    uk = (i - 1) * k + (1 : k)
    XTs[i, ] = apply(XTk[uk, ], 2, sum)
    TTs[i, ] = apply(TTk[uk, ], 2, sum)
    Tys[i, ] = apply(Tyk[uk, ], 2, sum)
    u2 = (i - 1) * 2 + (1 : 2)
    AtS[i, ] = c(crossprod(A21[u2, ], SS[u2, ]))
  }
  
  # SNPs start here
  end.oneoff.time = Sys.time()
  
  SNPS = SNP.matrix
  ns = ncol(SNPS)
  
  # Add SNPs one by one and solve
  Theta = D = matrix(NA, ns, 2)
  for (i in 1 : ns) {
    si = SNPS[, i]
    snp2 = rep(si, each = 2)
    H1 = matrix(crossprod(si, XTs), nx, 2)
    H2 = snp2 * SS
    AtH = matrix(crossprod(si, AtS), nx, 2)
    R = H1 - AtH
    Cfix = solve(Q, R)
    Cran = H2 - A21 %*% Cfix
    GtG = matrix(crossprod(si ^ 2, TTs), 2, 2)
    Gty = matrix(crossprod(si, Tys), 2, 1)
    V = GtG - crossprod(H1, Cfix) - crossprod(H2, Cran)
    v = Gty - crossprod(H1, sol) - crossprod(H2, blups)
    Theta[i , ] = solve(V, v)
    D[i, ] = diag(solve(V))
  }
  
  # sig2 is sigma squared
  
  SE = sqrt(sig**2) * sqrt(D)
  Pval = 2 * pnorm(-abs(Theta / SE))
  
  ## end time here
  end.time = Sys.time()
  
  # create output dataframe
  gallop.out.df = data.frame(
    SNP = SNP.names,
    cs.est = Theta[,1],
    cs.se = SE[,1],
    cs.t = Theta[,1] / SE[,1],
    cs.pval = Pval[,1],
    lon.est = Theta[,2],
    lon.se = SE[,2],
    lon.t = Theta[,2] / SE[,2],
    lon.pval = Pval[,2],
    method = rep("GALLOP", length(SNP.names))
  )
  
  out.list = list(
    output = gallop.out.df,
    time.taken.oneoff = end.oneoff.time - start.time,
    time.taken.SNPs = end.time - end.oneoff.time,
    time.taken.total = end.time - start.time
  )
  return(out.list)
  
}

# non linear? Still only looking for linear effects?

################################################################################
# scebe functions
################################################################################

# baseline
scebe_function1 = function(phenoData, SNP.matrix, SNP.names, mod1){

  colnames(phenoData) = c("ID", "Sex", "Time", "Height", "Weight")
  
  phenoData$Time = phenoData$Time - mean(phenoData$Time, na.rm=T)
  phenoData = phenoData %>% dplyr::select(ID,Weight,Height,Sex,Time)
  
  N = nrow(SNP.matrix)
  noCov = ncol(SNP.matrix)
  n_measure = phenoData %>% dplyr::group_by(ID) %>% dplyr::summarize(n=length(ID))%>%dplyr::select(n)%>%unlist()
  IDs = phenoData$ID
  
  ## start time here
  start.time = Sys.time()
  
  flme0<-summary(mod1)
  # rf=ranef(mod1)$ID # remove if using nlme
  rf=random.effects(mod1) # remove if using lme4
  r.name=rownames(rf)
  rf1=data.frame(ID=as.numeric(r.name),rf) %>% arrange(ID)
  rf = rf1[, -which(names(rf1)=="ID")]
  # ind.parm <- merge(phenoData[!duplicated(phenoData$ID),c("ID",paste0("cov",1:noCov))],rf1,by="ID")
  # names(ind.parm)[c(2+noCov,3+noCov)]=c("Inter","Slope")
  names(rf)[c(1:2)]=c("Inter","Slope")
  
  # y1: rf$Inter,
  # y2: rf$Slope
  # X: as.matrix(genoData)
  
  ## NEBE with SNPS starts here <-----------------------------------------------
  int1.time = Sys.time()
  
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
  
  # nebe with SNPs ends here <--------------------------------------------------
  int2.time = Sys.time()
  
  spr = result
  
  Data=phenoData
  # R=matrix(as.numeric(VarCorr(mod1)$ID),ncol=2) # lmer
  R = diag(as.numeric(nlme::VarCorr(mod1)[1:2,1])) # lme4 with uncorrelated random effects
  Ip=diag(1,2)
  
  #allS=allW=list()
  #allM=allISMIS=list()	  
  alls11=alls12=alls21=alls22=NULL
  allw11=allw12=allw21=allw22=NULL
  allsms11=allsms12=allsms21=allsms22=NULL
  
  Time.lst = split(Data[["Time"]], Data[["ID"]])
  
  #print(Time.lst)
  
  for(i in 1:N){ ## this loop does not involve SNPs, I believe
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
  
  ebe<-spr[,c(1,4)]
  
  # SNP section starts here
  end.oneoff.time = Sys.time()
  
  myX<-SNP.matrix #unique(phenoData[,5:(noCov+4)]) ## SNPs start again here
  Xc<-apply(myX,2,function(x){scale(x,scale=F)})
  
  XcX<-Xc*myX
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
  Vgamma22=D22-Vgamma222 ###VAR of NEBE
  
  Vssh11<-invssh11*Vgamma11*invssh11+invssh12*Vgamma21*invssh11+invssh11*Vgamma12*invssh12+invssh12*Vgamma22*invssh12
  Vssh12<-invssh11*Vgamma11*invssh21+invssh12*Vgamma21*invssh21+invssh11*Vgamma12*invssh22+invssh12*Vgamma22*invssh22
  Vssh21<-invssh21*Vgamma11*invssh11+invssh22*Vgamma21*invssh11+invssh21*Vgamma12*invssh12+invssh22*Vgamma22*invssh12
  Vssh22<-invssh21*Vgamma11*invssh21+invssh22*Vgamma21*invssh21+invssh21*Vgamma12*invssh22+invssh22*Vgamma22*invssh22
  
  # addition for single SNP input ---------------------------------------------<
  #  ebe = matrix(ebe, ncol = 2)
  
  est.ssh.inter<-invssh11*ebe[,1]+invssh12*ebe[,2]
  est.ssh.slope<-invssh21*ebe[,1]+invssh22*ebe[,2]
  
  t.ssh.inter<-est.ssh.inter/sqrt(Vssh11)
  t.ssh.slope<-est.ssh.slope/sqrt(Vssh22)
  
  p.ssh.inter<-2*pnorm(abs(t.ssh.inter),lower.tail=F)
  p.ssh.slope<-2*pnorm(abs(t.ssh.slope),lower.tail=F)
  
  ## end time here
  end.time = Sys.time()
  
  output<-cbind(est.ssh.inter,sqrt(Vssh11),p.ssh.inter,t.ssh.inter,est.ssh.slope,sqrt(Vssh22),p.ssh.slope,t.ssh.slope)
  
  colnames(output)= c("est.sc.b1", "se.sc.b1", "p.sc.b1", "t.sc.b1",
                      "est.sc.b2", "se.sc.b2", "p.sc.b2", "t.sc.b2")
  
  out.scebe = data.frame(SNP = SNP.names, matrix(output[,c(1:2,4,3,5:6,8,7)], ncol = 8), method = rep("SCEBE", nrow(output)))
  row.names(out.scebe) = seq(1,nrow(out.scebe),1)
  colnames(out.scebe) = c("SNP","cs.est", "cs.se", "cs.t", "cs.pval", "lon.est", "lon.se", "lon.t", "lon.pval", "method")
  
  out.list = list(
    output = out.scebe,
    time.taken.oneoff = int1.time - start.time + end.oneoff.time - int2.time,
    time.taken.SNPs = end.time - end.oneoff.time + int2.time - int1.time,
    time.taken.total = end.time - start.time
  )
  
  return(out.list)
  
}

# tvc
scebe_function2 = function(phenoData, SNP.matrix, SNP.names, mod1){
  
  # phenoData = data.snp
  colnames(phenoData) = c("ID", "Sex", "tvc", "Time", "Height", "Weight")
  # mod.lm = phenoData %>% lm(formula = Weight ~ Sex) # move out of function to compare
  # phenoData$Weight = mod.lm$residuals 
  
  phenoData$Time = phenoData$Time - mean(phenoData$Time, na.rm=T)
  #phenoData = phenoData %>% dplyr::select(ID,Weight,Time)
  
  N = nrow(SNP.matrix)
  noCov = ncol(SNP.matrix)
  n_measure = phenoData %>% dplyr::group_by(ID) %>% dplyr::summarize(n=length(ID))%>%dplyr::select(n)%>%unlist()
  IDs = phenoData$ID
  
  ## start time here
  start.time = Sys.time()
  
  flme0<-summary(mod1)
  # rf=ranef(mod1)$ID # remove if using nlme
  rf=random.effects(mod1) # remove if using lme4
  r.name=rownames(rf)
  rf1=data.frame(ID=as.numeric(r.name),rf) %>% arrange(ID)
  rf = rf1[, -which(names(rf1)=="ID")]
  # ind.parm <- merge(phenoData[!duplicated(phenoData$ID),c("ID",paste0("cov",1:noCov))],rf1,by="ID")
  # names(ind.parm)[c(2+noCov,3+noCov)]=c("Inter","Slope")
  names(rf)[c(1:2)]=c("Inter","Slope")
  
  # y1: rf$Inter,
  # y2: rf$Slope
  # X: as.matrix(genoData)
  
  ## NEBE with SNPS starts here <-----------------------------------------------
  int1.time = Sys.time()
  
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
  
  # nebe with SNPs ends here <--------------------------------------------------
  int2.time = Sys.time()
  
  spr = result
  
  Data=phenoData
  # R=matrix(as.numeric(VarCorr(mod1)$ID),ncol=2) # lmer
  R = diag(as.numeric(nlme::VarCorr(mod1)[1:2,1])) # lme4 with uncorrelated random effects
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
  
  ebe<-spr[,c(1,4)]
  
  # SNP section starts here
  end.oneoff.time = Sys.time()
  
  myX<-SNP.matrix #unique(phenoData[,5:(noCov+4)])
  Xc<-apply(myX,2,function(x){scale(x,scale=F)})
  
  XcX<-Xc*myX
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
  Vgamma22=D22-Vgamma222 ###VAR of NEBE
  
  Vssh11<-invssh11*Vgamma11*invssh11+invssh12*Vgamma21*invssh11+invssh11*Vgamma12*invssh12+invssh12*Vgamma22*invssh12
  Vssh12<-invssh11*Vgamma11*invssh21+invssh12*Vgamma21*invssh21+invssh11*Vgamma12*invssh22+invssh12*Vgamma22*invssh22
  Vssh21<-invssh21*Vgamma11*invssh11+invssh22*Vgamma21*invssh11+invssh21*Vgamma12*invssh12+invssh22*Vgamma22*invssh12
  Vssh22<-invssh21*Vgamma11*invssh21+invssh22*Vgamma21*invssh21+invssh21*Vgamma12*invssh22+invssh22*Vgamma22*invssh22
  
  # addition for single SNP input ---------------------------------------------<
  #  ebe = matrix(ebe, ncol = 2)
  
  est.ssh.inter<-invssh11*ebe[,1]+invssh12*ebe[,2]
  est.ssh.slope<-invssh21*ebe[,1]+invssh22*ebe[,2]
  
  t.ssh.inter<-est.ssh.inter/sqrt(Vssh11)
  t.ssh.slope<-est.ssh.slope/sqrt(Vssh22)
  
  p.ssh.inter<-2*pnorm(abs(t.ssh.inter),lower.tail=F)
  p.ssh.slope<-2*pnorm(abs(t.ssh.slope),lower.tail=F)
  
  ## end time here
  end.time = Sys.time()
  
  output<-cbind(est.ssh.inter,sqrt(Vssh11),p.ssh.inter,t.ssh.inter,est.ssh.slope,sqrt(Vssh22),p.ssh.slope,t.ssh.slope)
  
  colnames(output)= c("est.sc.b1", "se.sc.b1", "p.sc.b1", "t.sc.b1",
                      "est.sc.b2", "se.sc.b2", "p.sc.b2", "t.sc.b2")
  
  out.scebe = data.frame(SNP = SNP.names, matrix(output[,c(1:2,4,3,5:6,8,7)], ncol = 8), method = rep("SCEBE", nrow(output)))
  row.names(out.scebe) = seq(1,nrow(out.scebe),1)
  colnames(out.scebe) = c("SNP","cs.est", "cs.se", "cs.t", "cs.pval", "lon.est", "lon.se", "lon.t", "lon.pval", "method")
  
  out.list = list(
    output = out.scebe,
    time.taken.oneoff = int1.time - start.time + end.oneoff.time - int2.time,
    time.taken.SNPs = end.time - end.oneoff.time + int2.time - int1.time,
    time.taken.total = end.time - start.time
  )
  
  return(out.list)
  
}

# non linear?

################################################################################
# lme functions
################################################################################

# baseline
lmm_function1 = function(phenoData, SNP.matrix, SNP.names){
  
  #phenoData = data.snp
  colnames(phenoData) = c("ID", "Sex", "Time", "Height", "Weight")
  phenoData$Time = phenoData$Time - mean(phenoData$Time, na.rm=T)  # centre time to "To deal with potential multicolinearity caused by adding interaction effect of SNP x time"
  
  ns = ncol(SNP.matrix)
  out.lmm = matrix(rep(NA, ns*8), ncol = 8) 

  ## start time here
  start.time = Sys.time()
  
  for(i in 1:ncol(SNP.matrix)){
    
    #print(i)

    j.df = left_join(phenoData, data.frame(SNP.matrix[,i], ID = as.numeric(row.names(SNP.matrix))), by = "ID")
    colnames(j.df) = c("ID", "Sex", "Time", "Height", "Weight", "SNP")
    mod1 = suppressWarnings(lme(Weight ~ SNP*Time + Sex + Height,
                                random = list(ID = pdDiag(~Time)),
                                data = j.df,
                                control=c(msMaxIter=100, opt="optim"),
                                na.action="na.omit"))
    out.lmm[i,] = c(coef(summary(mod1))[2,c(1:2,4:5)],coef(summary(mod1))[6,c(1:2,4:5)])
  }
  
  ## end time here
  end.time = Sys.time()
  
  out.lmm = data.frame(SNP = SNP.names, out.lmm, method = rep("LMM", nrow(out.lmm)))
  colnames(out.lmm) = c("SNP","cs.est", "cs.se", "cs.t", "cs.pval", "lon.est", "lon.se", "lon.t", "lon.pval", "method")
  
  out.list = list(
    output = out.lmm,
    time.taken.total = end.time - start.time
  )
  
  return(out.list)
  
}

# tvc 
lmm_function2 = function(phenoData, SNP.matrix, SNP.names){
  
  #phenoData = data.snp
  colnames(phenoData) = c("ID", "Sex", "tvc", "Time", "Height", "Weight")
  phenoData$Time = phenoData$Time - mean(phenoData$Time, na.rm=T)  # centre time to "To deal with potential multicolinearity caused by adding interaction effect of SNP x time"
  
  ns = ncol(SNP.matrix)
  out.lmm = matrix(rep(NA, ns*8), ncol = 8) 
  
  ## start time here
  start.time = Sys.time()
  
  for(i in 1:ncol(SNP.matrix)){
    
    j.df = left_join(phenoData, data.frame(SNP.matrix[,i], ID = as.numeric(row.names(SNP.matrix))), by = "ID")
    colnames(j.df) = c("ID", "Sex", "tvc", "Time", "Height", "Weight", "SNP")
    mod1 = suppressWarnings(lme(Weight ~ SNP*Time + tvc*Time + Sex + Height,
                                random = list(ID = pdDiag(~Time)),
                                data = j.df, control=c(msMaxIter=100, opt="optim")))
    out.lmm[i,] = c(coef(summary(mod1))[2,c(1:2,4:5)],coef(summary(mod1))[7,c(1:2,4:5)])
  }
  
  ## end time here
  end.time = Sys.time()
  
  out.lmm = data.frame(SNP = SNP.names, out.lmm, method = rep("LMM", nrow(out.lmm)))
  colnames(out.lmm) = c("SNP","cs.est", "cs.se", "cs.t", "cs.pval", "lon.est", "lon.se", "lon.t", "lon.pval", "method")
  
  out.list = list(
    output = out.lmm,
    time.taken.total = end.time - start.time
  )
  
  return(out.list)
  
}

# non linear
lmm_function1.nonlin = function(phenoData, SNP.matrix, SNP.names){
  
  #phenoData = data.snp
  colnames(phenoData) = c("ID", "Sex", "Time", "Height", "Weight")
  phenoData$Time = phenoData$Time - mean(phenoData$Time, na.rm=T)  # centre time to "To deal with potential multicolinearity caused by adding interaction effect of SNP x time"
  
  ns = ncol(SNP.matrix)
  out.lmm = matrix(rep(NA, ns*8), ncol = 8) 
  
  ## start time here
  start.time = Sys.time()
  
  for(i in 1:ncol(SNP.matrix)){
    
    j.df = left_join(phenoData, data.frame(SNP.matrix[,i], ID = as.numeric(row.names(SNP.matrix))), by = "ID")
    colnames(j.df) = c("ID", "Sex", "Time", "Height", "Weight", "SNP")
    mod1 = suppressWarnings(lme(Weight ~ SNP*Time + I(Time**2) + Sex + Height,
                                random = list(ID = pdDiag(~Time)),
                                data = j.df, control=c(msMaxIter=100, opt="optim")))
    out.lmm[i,] = c(coef(summary(mod1))[2,c(1:2,4:5)],coef(summary(mod1))[7,c(1:2,4:5)])
  }
  
  ## end time here
  end.time = Sys.time()
  
  out.lmm = data.frame(SNP = SNP.names, out.lmm, method = rep("LMM", nrow(out.lmm)))
  colnames(out.lmm) = c("SNP","cs.est", "cs.se", "cs.t", "cs.pval", "lon.est", "lon.se", "lon.t", "lon.pval", "method")
  
  out.list = list(
    output = out.lmm,
    time.taken.total = end.time - start.time
  )
  
  return(out.list)
  
}

################################################################################
# semi factor functions
################################################################################

# non-linear
add_nonlinear_time_interaction = function(data, es, order=2, xshift=2, yshift=0){
  
  data = data %>% mutate(Weight = Weight + es * (Age.years+xshift)**order + yshift)
  
  return(data)
}



