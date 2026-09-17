## script to loop through all simulation combinations

# set parameter options
Cor.structure = c("Observed(ish)", "X to C greater", "C to X greater", "Increasing", "Decreasing")
Gamma = c(0, 1, 2, 25)
X.hypothesis = c("x1", "x2", "x3", "A", "NULL")
X.direct.effect = c(0.1)
C.pathway = c("c1", "c2", "c3", "All", "None")
C.direct.effect = c(0.05, 0.1, 0.2, 0.5)
Sample.size = 1000

setwd("Documents/...")

################################# functions ####################################

# function to simulate x and c for a given gamma
sim_xc = function(cov.matrix, gamma, sample.size = 100){
  
  n = sample.size
  
  Z <- matrix(rnorm(n*dim(cov.matrix)[1]), nrow=n)
  # Ensure these data are standardised and uncorrelated
  Z <- scale(Z)
  Z <- scale(Z)
  chol_decomp1 <- chol(var(Z))
  Z <-  Z %*% solve(chol_decomp1)
  
  # Generate data with the desired covariance structure
  chol_decomp2 <- chol(cov.matrix)
  gendata <- Z %*% chol_decomp2
  
  # Extract exposures and covariates
  X <- gendata[,c(1,3,5)]
  C <- gendata[,c(2,4,6)]
  
  # Generate new data with different correlation structure
  pred <- X %*% solve(t(X) %*% X) %*% t(X) %*% C # how does R read matrix multiplication?
  res <- C - pred
  chol_decomp3 <- chol(var(C))
  gamma_comb <- gamma * pred + res
  chol_decomp4 <- chol(var(gamma_comb))
  newC <- gamma_comb %*% solve(chol_decomp4) %*% chol_decomp3
  newdata <- cbind(X,newC)
  newdata <- newdata[,c(1,4,2,5,3,6)] # reorder columns to (x1,c1,x2,c2,x3,c3)
  
  return(data.frame(newdata))
}

# function to simulate outcome variable
simulate_outcome = function(xcdata, x.hypo = c(0,0,0,0,1), x.es, c.hypo = c(0,0,0,0,1), c.des = 1){
  
  beta_x1 = x.hypo[1]*x.es + x.hypo[4]*x.es*(1/3)
  beta_x2 = x.hypo[2]*x.es + x.hypo[4]*x.es*(1/3)
  beta_x3 = x.hypo[3]*x.es + x.hypo[4]*x.es*(1/3)
  
  beta_c1 = c.hypo[1]*c.des + c.hypo[4]*c.des*(1/3)
  beta_c2 = c.hypo[2]*c.des + c.hypo[4]*c.des*(1/3)
  beta_c3 = c.hypo[3]*c.des + c.hypo[4]*c.des*(1/3)
  
  xcdata$y = beta_x1*xcdata$x1 + beta_x2*xcdata$x2 + beta_x3*xcdata$x3 + beta_c1*xcdata$c1 + beta_c2*xcdata$c2 + beta_c3*xcdata$c3
  
  return(xcdata)
}

# regular slcma
reg_SLCMA = function(slcma_df){
  
  e1 <- lm(y~x1, data=slcma_df)$coef[2]
  e2 <- lm(y~x2, data=slcma_df)$coef[2]
  e3 <- lm(y~x3, data=slcma_df)$coef[2]
  eA <- lm(y~I(x1+x2+x3), data=slcma_df)$coef[2]
  
  s1 <- e1 * sd(slcma_df$x1)
  s2 <- e2 * sd(slcma_df$x2)
  s3 <- e3 * sd(slcma_df$x3)
  sA <- eA * sd(slcma_df$x1+slcma_df$x2+slcma_df$x3)
  
  return(list(effects = c(e1, e2, e3, eA), standardised = c(s1, s2, s3, sA)))
}

# direct and mediated effects slcma
dm_SLCMA = function(slcma_df){
  
  slcma_df$r1 <- lm(c1~x1, data=slcma_df)$resid
  slcma_df$r2 <- lm(c2~x1+x2+c1, data=slcma_df)$resid 
  slcma_df$r3 <- lm(c3~x1+x2+x3+c1+c2, data=slcma_df)$resid
  
  slcma_df$tx1 <- lm(x1~r1+r2+r3, data=slcma_df)$resid
  slcma_df$tx2 <- lm(x2~r1+r2+r3, data=slcma_df)$resid
  slcma_df$tx3 <- lm(x3~r1+r2+r3, data=slcma_df)$resid
  slcma_df$txA <- lm(I(x1+x2+x3)~r1+r2+r3, data=slcma_df)$resid
  
  te1 <- lm(y~tx1, data=slcma_df)$coef[2]
  te2 <- lm(y~tx2, data=slcma_df)$coef[2]
  te3 <- lm(y~tx3, data=slcma_df)$coef[2]
  teA <- lm(y~txA, data=slcma_df)$coef[2]
  
  ts1 <- te1 * sd(slcma_df$tx1)
  ts2 <- te2 * sd(slcma_df$tx2)
  ts3 <- te3 * sd(slcma_df$tx3)
  tsA <- teA * sd(slcma_df$txA)
  
  return(list(effects = c(te1, te2, te3, teA), standardised = c(ts1, ts2, ts3, tsA)))
}

# can loop through this
full_sim_function = function(cov.matrix, gamma, sample.size, x.hypo = c(0,0,0,0,1), x.es, c.hypo = c(0,0,0,0,1), c.des = 1){
  
  # sim exposure-covariate data
  xcdata = sim_xc(cov.matrix, gamma, sample.size)
  colnames(xcdata) = c("x1", "c1", "x2", "c2", "x3", "c3")
  
  # simulate outcome
  out.data = simulate_outcome(xcdata, x.hypo, x.es, c.hypo, c.des)
  
  # covariance matrix
  #cov.matrix.fsim = cov(out.data)
  
  # run SLCMAs
  reg.list = reg_SLCMA(out.data)
  dm.list = dm_SLCMA(out.data)
  
  # SLCMA effects table
  effects.table = t(data.frame(
    regular = unname(reg.list$effects[1:4]),
    dm = unname(dm.list$effects[1:4])#,
  ))
  colnames(effects.table) = c("x1", "x2", "x3", "A")
  
  # SLCMA standardised table
  standardised.table = t(data.frame(
    regular = unname(reg.list$standardised[1:4]),
    dm = unname(dm.list$standardised[1:4])#,
  ))
  colnames(standardised.table) = c("x1", "x2", "x3", "A")
  
  
  return(list(effect.estimates = effects.table, selection = standardised.table))
}

# reads in covariance matrix for each structure
read_cov = function(cor.structure){
  # read in covariance matrix
  if(cor.structure == "Observed(ish)"){
    cov.matrix = read.csv("adj.obs.cov.matrix.csv", row.names = 1)
  }
  if(cor.structure == "X to C greater"){
    cov.matrix = read.csv("XgC.cov.matrix.csv", row.names = 1)
  }
  if(cor.structure == "C to X greater"){
    cov.matrix = read.csv("CgX.cov.matrix.csv", row.names = 1)
  }
  if(cor.structure == "Increasing"){
    cov.matrix = read.csv("Inc.cov.matrix.csv", row.names = 1)
  }
  if(cor.structure == "Decreasing"){
    cov.matrix = read.csv("Dec.cov.matrix.csv", row.names = 1)
  }
  
  return(cov.matrix)
}

# set up x.hypo
set_x.hypo = function(x.hypothesis){
  
  if(x.hypothesis == "x1"){
    x.hypo = c(1,0,0,0,0)
  }else{
    if(x.hypothesis == "x2"){
      x.hypo = c(0,1,0,0,0)
    }else{
      if(x.hypothesis == "x3"){
        x.hypo = c(0,0,1,0,0)
      }else{
        if(x.hypothesis == "A"){
          x.hypo = c(0,0,0,1,0)
        }else{
          x.hypo = c(0,0,0,0,1)
        }
      }
    }
  }
  
  return(x.hypo)
}

# set up c.hypo
set_c.hypo = function(c.pathway){
  
  if(c.pathway == "c1"){
    c.hypo = c(1,0,0,0,0)
  }else{
    if(c.pathway == "c2"){
      c.hypo = c(0,1,0,0,0)
    }else{
      if(c.pathway == "c3"){
        c.hypo = c(0,0,1,0,0)
      }else{
        if(c.pathway == "All"){
          c.hypo = c(0,0,0,1,0)
        }else{
          c.hypo = c(0,0,0,0,1)
        }
      }
    }
  }
  
  return(c.hypo)
}

################################# loops ########################################

# make empty output dataframe
out.df = data.frame(
  cor.structure = character(),
  x.hypothesis = character(),
  x.direct.effect = numeric(),
  c.pathway = character(),
  c.direct.effect = numeric(),
  gamma = numeric(),
  sample.size = numeric(),
  agree = character(),
  reg.choice = character(),
  dm.choice = character(),
  reg.sta.x1 = numeric(),
  reg.sta.x2 = numeric(),
  reg.sta.x3 = numeric(),
  reg.sta.A = numeric(),
  dm.sta.x1 = numeric(),
  dm.sta.x2 = numeric(),
  dm.sta.x3 = numeric(),
  dm.sta.A = numeric(),
  reg.estimate.x1 = numeric(),
  reg.estimate.x2 = numeric(),
  reg.estimate.x3 = numeric(),
  reg.estimate.A = numeric(),
  dm.estimate.x1 = numeric(),
  dm.estimate.x2 = numeric(),
  dm.estimate.x3 = numeric(),
  dm.estimate.A = numeric()
)

# top level cor structure
for(strct in Cor.structure){
  
  # read in required cov matrix
  cov.matrix = read_cov(strct)
  
  # next level x.hypo
  for(item in X.hypothesis){
    x.hypo = set_x.hypo(item)
    
    # check progress
    print(paste0(strct, "    -    ", item))
    
    # next level x effect size (only one level)
    for(x.effect in X.direct.effect){
      
        # next level c.hypo
      for(jtem in C.pathway){
        c.hypo = set_c.hypo(jtem)
        
        # next level c.direct.effect
        for(c.effect in C.direct.effect){
          
          # last level gamma
          for(gtem in Gamma){
            
            # run full sim function
            sim.out = full_sim_function(cov.matrix, gtem, Sample.size, x.hypo, x.effect, c.hypo, c.effect)
            
            # get selections
            choices = colnames(sim.out$selection)[max.col(sim.out$selection**2,ties.method="first")]
            
            
            # make output row
            out.df[nrow(out.df) + 1,] = c(strct, item, x.effect, jtem, c.effect, gtem, Sample.size,
                                          ifelse(choices[1] == choices[2], "yes", "no"),
                                          choices[1], choices[2],
                                          sim.out$selection[1,1], sim.out$selection[1,2], sim.out$selection[1,3], sim.out$selection[1,4],
                                          sim.out$selection[2,1], sim.out$selection[2,2], sim.out$selection[2,3], sim.out$selection[2,4],
                                          sim.out$effect.estimates[1,1], sim.out$effect.estimates[1,2], sim.out$effect.estimates[1,3], sim.out$effect.estimates[1,4],
                                          sim.out$effect.estimates[2,1], sim.out$effect.estimates[2,2], sim.out$effect.estimates[2,3], sim.out$effect.estimates[2,4]
                                          
                                          )
            
          }
        }
      }
    }
  }
}

write.csv(out.df, file = "out.df.csv", row.names = F)

out.df[out.df$agree == "no" & out.df$gamma != 25 & out.df$x.hypothesis != "NULL",1:10]

