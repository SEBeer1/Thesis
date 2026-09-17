## noisy simulation script

# cor, x.hypo, c.hypo, c.des, gamma, sample size, r^2
baseline = c("C to X greater", "x2", "c1", 0.2, 2, 500, 0.1) # numbers converted to string here
x.des = c(0.1) # stays constant for all

# delete baseline from options
cor.structure.options = c("Observed(ish)", "X to C greater", "Increasing", "Decreasing")
x.hypo.options = c("x1", "x3", "A", "NULL")
c.hypo.options = c("c2", "c3", "All", "None")
c.hypo.des = c(0.05, 0.1, 0.5)
gamma.options = c(1,25)
sample.size.options = c(100,1000)
rsquared.options = c( 0.01)

setwd("/Documents/...")

### functions

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
# c.scale is actually c direct effect size
simulate_outcome_with_noise = function(xcdata, x.hypo = c(0,0,0,0,1), x.es, c.hypo = c(0,0,0,0,1), c.scale = 1, trsq = 0.1){
  
  beta_x1 = x.hypo[1]*x.es + x.hypo[4]*x.es*(1/3)
  beta_x2 = x.hypo[2]*x.es + x.hypo[4]*x.es*(1/3)
  beta_x3 = x.hypo[3]*x.es + x.hypo[4]*x.es*(1/3)
  
  beta_c1 = c.hypo[1]*c.scale + c.hypo[4]*c.scale*(1/3)
  beta_c2 = c.hypo[2]*c.scale + c.hypo[4]*c.scale*(1/3)
  beta_c3 = c.hypo[3]*c.scale + c.hypo[4]*c.scale*(1/3)
  
  xcdata$y = beta_x1*xcdata$x1 + beta_x2*xcdata$x2 + beta_x3*xcdata$x3 + beta_c1*xcdata$c1 + beta_c2*xcdata$c2 + beta_c3*xcdata$c3
  
  xcdata = scale(xcdata)
  xcdata = data.frame(xcdata)
  xcdata$y = xcdata$y + rnorm(length(xcdata$y))*sqrt((1-trsq)/trsq)
  
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
full_sim_function_noisy = function(cov.matrix, gamma, sample.size, x.hypo = c(0,0,0,0,1), x.es, c.hypo = c(0,0,0,0,1), c.scale = 1, trsq = 0.1){
  
  # sim exposure-covariate data
  xcdata = sim_xc(cov.matrix, gamma, sample.size)
  colnames(xcdata) = c("x1", "c1", "x2", "c2", "x3", "c3")
  
  # simulate outcome
  out.data = simulate_outcome_with_noise(xcdata, x.hypo, x.es, c.hypo, c.scale, trsq)
  
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

# make table from full sim output
run_nsim_for_one_combo = function(cov.matrix, gtem, sample.size, x.hypo, x.des, c.hypo, c.des, trsq, nsim, strct, item, jtem){ ## NEED TO SORT OUT ITEM JTEM!!
  
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
    dm.estimate.A = numeric(),
    iteration = numeric()
  )
  
  seed.store = vector("list", nsim)
  
    
  
  
  
  for(i in 1:nsim){
    
    # store starting seed
    seed1 = .Random.seed
    
    # simulate data and apply SLCMAs
    sim.out = full_sim_function_noisy(cov.matrix, gtem, sample.size, x.hypo, x.des, c.hypo, c.des, trsq)
    
    # store ending seed
    seed2 = .Random.seed
  
  
    # get choices
    choices = colnames(sim.out$selection)[max.col(sim.out$selection**2,ties.method="first")]
    
    
    # make output row
    out.df[nrow(out.df) + 1,] = c(strct, item, x.des, jtem, c.des, gtem, sample.size,
                                  ifelse(choices[1] == choices[2], "yes", "no"),
                                  choices[1], choices[2],
                                  sim.out$selection[1,1], sim.out$selection[1,2], sim.out$selection[1,3], sim.out$selection[1,4],
                                  sim.out$selection[2,1], sim.out$selection[2,2], sim.out$selection[2,3], sim.out$selection[2,4],
                                  sim.out$effect.estimates[1,1], sim.out$effect.estimates[1,2], sim.out$effect.estimates[1,3], sim.out$effect.estimates[1,4],
                                  sim.out$effect.estimates[2,1], sim.out$effect.estimates[2,2], sim.out$effect.estimates[2,3], sim.out$effect.estimates[2,4],
                                  i
                                  
    )
    
    Ps = list(iteration = i,
              start.seed = seed1,
              end.seed = seed2)
    
    seed.store[[i]] <- Ps
  
  }
    
    setNames(seed.store, paste0("iteration", 1:nsim))
    
    return(list(out.df = out.df, seed.list = seed.store))
}

############################ observed ish baseline #############################
# cor, x.hypo, c.hypo, c.des, gamma, sample size, r^2
baseline = c("Observed(ish)", "x2", "c2", 0.2, 2, 500, 0.1) # numbers converted to string here - be carefuk when using later
x.des = c(0.1) # stays constant for all

# delete baseline from options
cor.structure.options = c("X to C greater", "C to X greater", "Increasing", "Decreasing")
x.hypo.options = c("x1", "x3", "A", "NULL")
c.hypo.options = c("c1", "c3", "All") 
c.hypo.des = c(0.05, 0.1, 0.5)
gamma.options = c(1,25)
sample.size.options = c(100,1000)
rsquared.options = c( 0.01)


run_noisy_sims_baseline3 = function(baseline, alt.cor, alt.x.hypo, alt.c.hypo, alt.c.hypo.des, alt.gamma, alt.sample.size, alt.rsquared, x.des, nsim){
  
  # run baseline
  
  # set up hypos
  x.hypo = set_x.hypo(baseline[2])
  c.hypo = set_c.hypo(baseline[3])
  
  c1.out = run_nsim_for_one_combo(read_cov(baseline[1]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                  as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, baseline[1], baseline[2], baseline[3])
  
  # loop through cor.structures
  c2.out = run_nsim_for_one_combo(read_cov(alt.cor[1]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                  as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, alt.cor[1], baseline[2], baseline[3])
  
  c3.out = run_nsim_for_one_combo(read_cov(alt.cor[2]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                  as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, alt.cor[2], baseline[2], baseline[3])
  
  c4.out = run_nsim_for_one_combo(read_cov(alt.cor[3]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                  as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, alt.cor[3], baseline[2], baseline[3])
  
  c5.out = run_nsim_for_one_combo(read_cov(alt.cor[4]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                  as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, alt.cor[4], baseline[2], baseline[3])
  
  print("5 done...")
  
  # loop through x.hypos
  x.hypo = set_x.hypo(alt.x.hypo[1])
  c6.out = run_nsim_for_one_combo(read_cov(baseline[1]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                  as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, baseline[1], alt.x.hypo[1], baseline[3])
  
  x.hypo = set_x.hypo(alt.x.hypo[2])
  c7.out = run_nsim_for_one_combo(read_cov(baseline[1]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                  as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, baseline[1], alt.x.hypo[2], baseline[3])
  
  x.hypo = set_x.hypo(alt.x.hypo[3])
  c8.out = run_nsim_for_one_combo(read_cov(baseline[1]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                  as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, baseline[1], alt.x.hypo[3], baseline[3])
  
  x.hypo = set_x.hypo(alt.x.hypo[4])
  c9.out = run_nsim_for_one_combo(read_cov(baseline[1]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                  as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, baseline[1], alt.x.hypo[4], baseline[3])
  
  x.hypo = set_x.hypo(baseline[2])
  
  # loop through c.hypos
  c.hypo = set_c.hypo(alt.c.hypo[1])
  c10.out = run_nsim_for_one_combo(read_cov(baseline[1]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                   as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, baseline[1], baseline[2], alt.c.hypo[1])
  
  print("10 done...")
  
  c.hypo = set_c.hypo(alt.c.hypo[2])
  c11.out = run_nsim_for_one_combo(read_cov(baseline[1]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                   as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, baseline[1], baseline[2], alt.c.hypo[2])
  
  c.hypo = set_c.hypo(alt.c.hypo[3])
  c12.out = run_nsim_for_one_combo(read_cov(baseline[1]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                   as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, baseline[1], baseline[2], alt.c.hypo[3])
  
  c.hypo = set_c.hypo(baseline[3])
  
  # loop through c.des
  c14.out = run_nsim_for_one_combo(read_cov(baseline[1]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                   alt.c.hypo.des[1], as.numeric(baseline[7]), nsim, baseline[1], baseline[2], baseline[3])
  
  c15.out = run_nsim_for_one_combo(read_cov(baseline[1]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                   alt.c.hypo.des[2], as.numeric(baseline[7]), nsim, baseline[1], baseline[2], baseline[3])
  
  print("15 done...")
  
  c16.out = run_nsim_for_one_combo(read_cov(baseline[1]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                   alt.c.hypo.des[3], as.numeric(baseline[7]), nsim, baseline[1], baseline[2], baseline[3])
  
  # loop through other gammas
  c17.out = run_nsim_for_one_combo(read_cov(baseline[1]), alt.gamma[1], as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                   as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, baseline[1], baseline[2], baseline[3])
  
  c18.out = run_nsim_for_one_combo(read_cov(baseline[1]), alt.gamma[2], as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                   as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, baseline[1], baseline[2], baseline[3])
  
  # loop through sample sizes
  c19.out = run_nsim_for_one_combo(read_cov(baseline[1]), as.numeric(baseline[5]), alt.sample.size[1], x.hypo, x.des, c.hypo,
                                   as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, baseline[1], baseline[2], baseline[3])
  
  c20.out = run_nsim_for_one_combo(read_cov(baseline[1]), as.numeric(baseline[5]), alt.sample.size[2], x.hypo, x.des, c.hypo,
                                   as.numeric(baseline[4]), as.numeric(baseline[7]), nsim, baseline[1], baseline[2], baseline[3])
  
  # loop through rsquared
  c21.out = run_nsim_for_one_combo(read_cov(baseline[1]), as.numeric(baseline[5]), as.numeric(baseline[6]), x.hypo, x.des, c.hypo,
                                   as.numeric(baseline[4]), alt.rsquared[1], nsim, baseline[1], baseline[2], baseline[3])
  
  # output: list of 1000 row dataframes
  return(list(
    baseline = c1.out, 
    XtC = c2.out, CtX = c3.out, Inc = c4.out, Dec = c5.out,
    x1 = c6.out, x2 = c7.out, x3 = c8.out, all.x = c9.out,
    c1 = c10.out, c3 = c11.out, all.cc = c12.out,
    cdes.0.05 = c14.out, cdes.0.1 = c15.out, cdes.0.5 = c16.out,
    gamma1 = c17.out, gamma25 = c18.out,
    samples100 = c19.out, samples1000 = c20.out,
    rsq0.001 = c21.out
  ))
}

set.seed(2512)
full_run = run_noisy_sims_baseline3(baseline, cor.structure.options, x.hypo.options, c.hypo.options, c.hypo.des, gamma.options, sample.size.options, rsquared.options, 0.1, 1000)
saveRDS(full_run, file = "baselineFullRun.rds")
