## fit mod1 longitudinal mQTL script
#!/usr/bin Rscript


# ---- Packages ----
suppressMessages(library(nlme))
suppressMessages(library(lme4))
suppressMessages(library(lmerTest))
suppressMessages(library(nloptr))
suppressMessages(library(estimatr))
suppressMessages(library(argparse))
suppressMessages(library(data.table))
suppressMessages(library(foreach))
suppressMessages(library(parallel))
suppressMessages(library(doMC))
suppressMessages(library(dplyr))
suppressMessages(library(genio))
suppressMessages(library(miceadds))
suppressMessages(library(performance))

registerDoMC(5)


# ---- Set Chromosome for Analysis ----
parser <- ArgumentParser()
parser$add_argument("--chr", type="character")
Args <- parser$parse_args()
chr <- Args$chr

# ---- Start ----
timer <- system.time({
  
  print(paste0("Fit mod1 per CpG: Chromosome ", chr))
  
  print("Read Input")
  
  # update to my candidate list
  candidates <- fread(paste0("/data2/sbeer/lmQTL/DNAm/DNAmSubsets/CpGbyChr/cpgChr", chr, ".csv"))
  candidates <- candidates[,1] 
  
  # ---- Function for Analysis ----
  
  model_one <- function(cpg.name, dnam) {
    
    dnam <- data.frame(Sample_Name = names(dnam),
                       cpg = dnam)
    
    cov <- fread("/data2/sbeer/lmQTL/covariates/DCHS_covariates_longitudinal.csv")
    
    # scebe requires 2+ outcome measures per individual
    nt = cov %>% dplyr::group_by(ID) %>% dplyr::summarize(n=length(ID)) %>% dplyr::arrange(ID)
    twoplus = nt %>% dplyr::filter(n > 1) %>% dplyr::select(ID) %>% unlist()
    n_measure = nt %>% dplyr::filter(n > 1) %>% dplyr::select(n) %>% unlist()

    data <- merge(cov, dnam, by = "Sample_Name")

    data = data %>% dplyr::filter(ID %in% twoplus)
    
    # Format Data
    data$ID <- as.factor(data$ID)
    data$Age_Years <- as.numeric(data$Age_Years)
    data$sex <- as.factor(data$sex)
    data$site <- as.factor(data$site)
    data$genotype_batch <- as.factor(data$genotype_batch)
    data$Sample_Plate <- as.factor(data$Sample_Plate)
    data$Sample_col <- as.factor(data$Sample_col)
    data$Slide <- as.factor(data$Slide)
    data$sentrix_row <- as.factor(data$sentrix_row)
    #data$snp <- as.numeric(data$snp)
    data$cpg <- as.numeric(data$cpg)
    
    # Centre time: set intercept to age 1
    data$Age_Years <- data$Age_Years - 1 # mean(data$Age_Years)  
        
    # Residualizing
    m1 <- lm(as.formula("cpg ~ Bcell + CD4T + CD8T + Eos + Mono + Neu + NK + Slide + Sample_Plate + Sample_col + sentrix_row"),
             data = data)
    
    data$resid <- resid(m1)
  
    # Mixed Effect Model # 1-lmer_cor_rf, 2-lmer_uncor_rf, 3-nlme
    # f <- as.formula("resid ~ Age_Years + sex + site + genotype_batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + (Age_Years|ID)")
    # f <- as.formula("resid ~ Age_Years + sex + site + genotype_batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + (Age_Years||ID)")
    f <- as.formula("resid ~ Age_Years + sex + site + genotype_batch + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10")
    
    m2 <- NULL
    
    # 1-lmer, 2-nlme_uncor_rf, 3-nlme_cor_rf
    # m2 <- lmer(f, data = data, control = lmerControl(optimizer = "nloptwrap", calc.derivs = FALSE))
    # m2 <- lme(f, random = list(ID = pdDiag(~Age_Years)), data = data, control=c(msMaxIter=100, opt="optim"), na.action = "na.omit")
    m2 <- lme(f, random = ~Age_Years|ID, data = data, control=c(msMaxIter=100, opt="optim"), na.action = "na.omit")

    # Singular flag
    singular = ifelse(check_singularity(m2, tolerance = 1e-6), T, F)

    # remove once working
    # print(cpg.name[[1]])
	
    if(singular == F){
      # Output
    return(list(cpg.name = cpg.name,
                mod1 = m2,
                res.cpg = resid(m1),
                singular = singular))
    }
    else{
      return(list(cpg.name = cpg.name,
                  singular = singular))
    }
    
  }
  
  
  # ---- Main Analysis ----
  print("Run mod1 fitting")
  
  # Load Data
  load.Rdata(paste0("/data2/sbeer/lmQTL/DNAm/DCHS_betas_longitudinal_chr", chr, ".Rdata"),
             "dnam")

  # Run in Parallel
  output <- foreach(i = 1:nrow(candidates), .export = 'model_one') %dopar% {
    
    tryCatch({
    	# Run Analysis
    	model_one(candidates[i,], dnam[candidates$CpG[i],])
    }, error = function(msg){
        print(paste0("Failed to fit CpG: ", candidates[i,]))
	return(list(cpg.name = candidates[i,],
		    singular = TRUE))
    })
  }
  
  # Write Output
  save(output, file = paste0("/data2/sbeer/lmQTL/mod1andresid/Chr", chr, "Output.Rdata"))
  
})


elapsed_time <- timer["elapsed"]

print(paste0("Models fitted in ",
             elapsed_time/60,
             " minutes"))
