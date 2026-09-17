# slcma functions
## slcma_df a data.frame with columns x1, x2, x3, c1, c2, c3, y
## output effects are the effect estimates for each hypothesis
## output standardised are standardised between hypotheses for selection
## regular SLCMA
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

## direct and mediated effects SLCMA
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

# empirical application
## selected accumulation model
slcma_df$r1 <- lm(c1~x1, data=slcma_df)$resid
slcma_df$r2 <- lm(c2~x1+x2+c1, data=slcma_df)$resid 
slcma_df$r3 <- lm(c3~x1+x2+x3+c1+c2, data=slcma_df)$resid

slcma_df$tx1 <- lm(x1~r1+r2+r3, data=slcma_df)$resid
slcma_df$tx2 <- lm(x2~r1+r2+r3, data=slcma_df)$resid
slcma_df$tx3 <- lm(x3~r1+r2+r3, data=slcma_df)$resid
slcma_df$txA <- lm(I(x1+x2+x3)~r1+r2+r3, data=slcma_df)$resid

selected.model = lm(y~txA, data=slcma_df)
summary(selected.model)

## post selection inference
library(slcma)

tx1 = slcma_df$tx1
tx2 = slcma_df$tx2
tx3 = slcma_df$tx3
y = slcma_df$y

myslcma = slcma(y ~ tx1 + tx2 + tx3 + Accumulation(tx1,tx2,tx3))
slcmaInfer(myslcma, 1) # selective inference
slcmaInfer(myslcma, 1, method = "maxt", do.maxtCI=TRUE) # max-|t|
