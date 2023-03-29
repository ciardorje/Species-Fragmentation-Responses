rm(list=ls()); gc()

library(pacman)
p_load(tidyverse, nimble, wiqid, parallel, doSNOW, abind, MCMCvis)

setwd('C:/Users/ciarn/Desktop/PhD/Dung_Beetles/')
set.seed(56)

load('AF_DB_HMSOM_Data.RData')
dbs_wide <- dbs_wide_incidence

#-----Prepare Data-----
Ymat <- dbs_array
SppBinoms <- colnames(dbs_wide)[36:ncol(dbs_wide)]
nTraps <- as.numeric(dbs_wide$Fragment_Samples)
nSpObs <- dim(Ymat)[2]
nPatches <- dim(Ymat)[1]

Trap_Habitat[Trap_Habitat == 'Core'] <- 1
Trap_Habitat[Trap_Habitat == 'Edge'] <- 2
Trap_Habitat <- apply(Trap_Habitat, 2, as.integer)

z <- (Ymat > 0) * 1
z <- apply(z, 1:2, sum, na.rm = T)
z <- ifelse(z > 0, 1, NA)

#Variables for comparison 
covars <- dbs_wide[,15:22]
covars <- covars %>% select(Area_Ha, Shape, 
                            Prox_100m, Prox_250m, Prox_500m, 
                            Prox_750m, Prox_1000m, Prox_2500m)
covars$Area_Ha <- log10(covars$Area_Ha)
covars[,3:8] <- log10(covars[,3:8] + 1)
covars_std <- as.data.frame(apply(covars, 2, standardize))

#-----Model Structure-----
candidate_model <- nimbleCode({
  
  ###Species-Level Model###
  for(i in 1:nSpObs){
    for(j in 1:nPatches){
      
      #Occurence Model
      logit(psi[j, i]) <- lpsi[i] + bCoVar[1, i] * Area[j] + bCoVar[2, i] * Shape[j] + bCoVar[3, i] * Prox[j]
      z[j, i] ~ dbern(psi[j, i])
      
      #Detection Model
      for(k in 1:nTraps[j]){
        
        logit(p[j, i, k]) <- lp[i] + aHabitat[i, Habitat[j, k]] #Inc. habitat as fixed effect
        y[j, i, k] ~ dbern(p[j, i, k] * z[j, i])
      
      }
    }  
    
    ####Species Level Priors###
    lpsi[i] ~ dnorm(mu.lpsi, tau.lpsi)
    
    for(x in 1:3){
      bCoVar[x, i] ~ dnorm(mu.CoVar[x], tau.CoVar[x])
    }
    
    mu.eta[i] <- mu.lp + rho * sd.lp/sd.lpsi * (lpsi[i] - mu.lpsi)
    lp[i] ~ dnorm(mu.eta[i], tau.eta)
    
    aHabitat[i, 1] <- 0 #Core (Habitat = 1) set as reference habitat (lp = detection probability in core)
    aHabitat[i, 2] ~ dnorm(mu.Habitat, tau.Habitat) #Prior for difference in detection at patch edge
    
  }
  
  ###Community-Level Hyperpriors###
  
  #Occurence
  psi.mean ~ dbeta(1, 1)
  mu.lpsi <- logit(psi.mean)
  sd.lpsi ~ dunif(0, 5)
  tau.lpsi <- 1/sd.lpsi^2
  
  for(x in 1:3){
    mu.CoVar[x] ~ dnorm(0,0.1) #Covariates are standardised
    sd.CoVar[x] ~ dunif(0, 5)
    tau.CoVar[x] <- 1/sd.CoVar[x]^2
  }
  
  p.mean ~ dbeta(1, 1)
  mu.lp <- logit(p.mean)
  sd.lp ~ dunif(0, 5)
  tau.lp <- 1/sd.lp^2
  
  mu.Habitat ~ dnorm(0,0.1)
  sd.Habitat ~ dunif(0, 5)
  tau.Habitat <- 1/sd.Habitat^2
    
  #Relationship between occurrence and detection 
  #(more common species are [generally] more likely to be detected)
  rho ~ dunif(-1, 1) 
  tau.eta <- tau.lp/(1 - rho^2)
    
})

#-----Run Buffer Size Selection Using WAIC-----

prox_buffers <- c('Prox100', 'Prox250', 'Prox500', 'Prox750', 'Prox1000', 'Prox2500')

model_inputs <- list()

for(i in 1:length(prox_buffers)){
  
  data <- list(y = Ymat, nSpObs = nSpObs, nTraps = nTraps, nPatches = nPatches,
               z = z, Habitat = Trap_Habitat, Area = covars_std$Area_Ha, 
               Shape = covars_std$Shape, Prox = covars_std[,(2+i)])
  model_inputs[[i]] <- data
  names(model_inputs)[i] <- prox_buffers[i]
  
}

WAICmodel <- function(x){
  
  model <- nimbleModel(candidate_model, constants = x)
  MCMCconf <- configureMCMC(model, enableWAIC = T, 
                            monitors = c('lpsi', 'lp', 'aHabitat', 'bCoVar'))
  MCMC <- buildMCMC(MCMCconf)
  compModel <- compileNimble(model)
  compMCMC <- compileNimble(MCMC, project = compModel)
  results <- runMCMC(compMCMC, niter = 60000, nburnin = 10000, 
                     thin = 10, nchains = 4, WAIC = T, samplesAsCodaMCMC = T)
  return(results)
  
}

my.cluster <- makeCluster(6)
registerDoSNOW(my.cluster)
ntasks <- length(model_inputs) #set up progress bar
pb <- txtProgressBar(max = ntasks, style = 3)
progress <- function(x) setTxtProgressBar(pb, x)
opts <- list(progress = progress) 

results <- foreach(i = 1:length(model_inputs),
                   .options.snow = opts,
                   .packages = 'nimble'
) %dopar% {
  
  WAICmodel(model_inputs[[i]])
  
}

stopCluster(my.cluster)

names(results) <- prox_buffers
WAIC_results <- data.frame(Variable = prox_buffers, WAIC = NA)
WAIC_results$WAIC <- as.numeric(lapply(results, function(x) x[['WAIC']][['WAIC']]))

#-----Find best model including predetermined best Prox buffer-----

load('AF_DB_HMSOM_WAIC_Results.RData')

#Full Model
full_model_waic <- WAIC_results[WAIC_results$WAIC == min(WAIC_results$WAIC),]

#Univariate models
univariate_model <- nimbleCode({
  
  ###Species-Level Model###
  for(i in 1:nSpObs){
    for(j in 1:nPatches){
      
      #Occurence Model
      logit(psi[j, i]) <- lpsi[i] + bCoVar[i] * CoVar[j]
      z[j, i] ~ dbern(psi[j, i])
      
      #Detection Model
      for(k in 1:nTraps[j]){
        
        logit(p[j, i, k]) <- lp[i] + aHabitat[i, Habitat[j, k]] #Inc. habitat as fixed effect
        y[j, i, k] ~ dbern(p[j, i, k] * z[j, i])
        
      }
    }  
    
    ####Species Level Priors###
    lpsi[i] ~ dnorm(mu.lpsi, tau.lpsi)
    
    bCoVar[i] ~ dnorm(mu.CoVar, tau.CoVar)
    
    mu.eta[i] <- mu.lp + rho * sd.lp/sd.lpsi * (lpsi[i] - mu.lpsi)
    lp[i] ~ dnorm(mu.eta[i], tau.eta)
    
    aHabitat[i, 1] <- 0 #Core (Habitat = 1) set as reference habitat (lp = detection probability in core)
    aHabitat[i, 2] ~ dnorm(mu.Habitat, tau.Habitat) #Prior for difference in detection at patch edge
    
  }
  
  ###Community-Level Hyperpriors###
  
  #Occurence
  psi.mean ~ dbeta(1, 1)
  mu.lpsi <- logit(psi.mean)
  sd.lpsi ~ dunif(0, 5)
  tau.lpsi <- 1/sd.lpsi^2
  
  
  mu.CoVar ~ dnorm(0,0.1) #Covariates are standardised
  sd.CoVar ~ dunif(0, 5)
  tau.CoVar <- 1/sd.CoVar^2
  
  p.mean ~ dbeta(1, 1)
  mu.lp <- logit(p.mean)
  sd.lp ~ dunif(0, 5)
  tau.lp <- 1/sd.lp^2
  
  mu.Habitat ~ dnorm(0,0.1)
  sd.Habitat ~ dunif(0, 5)
  tau.Habitat <- 1/sd.Habitat^2
  
  #Relationship between occurrence and detection 
  #(more common species are [generally] more likely to be detected)
  rho ~ dunif(-1, 1) 
  tau.eta <- tau.lp/(1 - rho^2)
  
})

single_variables <- list(Area = covars_std$Area_Ha, Prox = covars_std$Prox_500m,
                         Shape = covars_std$Shape)

univariate_inputs <- list()
for(i in 1:length(single_variables)){
  univariate_inputs[[i]] <- list(y = Ymat, nSpObs = nSpObs, nTraps = nTraps, nPatches = nPatches,
                                 z = z, Habitat = Trap_Habitat, CoVar = single_variables[[i]])
  names(univariate_inputs)[i] <- names(single_variables)[i]
}

#Bivariate models
bivariate_model <- nimbleCode({
  
  ###Species-Level Model###
  for(i in 1:nSpObs){
    for(j in 1:nPatches){
      
      #Occurence Model
      logit(psi[j, i]) <- lpsi[i] + bCoVar[1, i] * CoVar1[j] + bCoVar[2, i] * CoVar2[j]
      z[j, i] ~ dbern(psi[j, i])
      
      #Detection Model
      for(k in 1:nTraps[j]){
        
        logit(p[j, i, k]) <- lp[i] + aHabitat[i, Habitat[j, k]] #Inc. habitat as fixed effect
        y[j, i, k] ~ dbern(p[j, i, k] * z[j, i])
        
      }
    }  
    
    ####Species Level Priors###
    lpsi[i] ~ dnorm(mu.lpsi, tau.lpsi)
    
    for(x in 1:2){
      bCoVar[x, i] ~ dnorm(mu.CoVar[x], tau.CoVar[x])
    }
    
    mu.eta[i] <- mu.lp + rho * sd.lp/sd.lpsi * (lpsi[i] - mu.lpsi)
    lp[i] ~ dnorm(mu.eta[i], tau.eta)
    
    aHabitat[i, 1] <- 0 #Core (Habitat = 1) set as reference habitat (lp = detection probability in core)
    aHabitat[i, 2] ~ dnorm(mu.Habitat, tau.Habitat) #Prior for difference in detection at patch edge
    
  }
  
  ###Community-Level Hyperpriors###
  
  #Occurence
  psi.mean ~ dbeta(1, 1)
  mu.lpsi <- logit(psi.mean)
  sd.lpsi ~ dunif(0, 5)
  tau.lpsi <- 1/sd.lpsi^2
  
  for(x in 1:2){
    mu.CoVar[x] ~ dnorm(0,0.1) #Covariates are standardised
    sd.CoVar[x] ~ dunif(0, 5)
    tau.CoVar[x] <- 1/sd.CoVar[x]^2
  }
  
  p.mean ~ dbeta(1, 1)
  mu.lp <- logit(p.mean)
  sd.lp ~ dunif(0, 5)
  tau.lp <- 1/sd.lp^2
  
  mu.Habitat ~ dnorm(0,0.1)
  sd.Habitat ~ dunif(0, 5)
  tau.Habitat <- 1/sd.Habitat^2
  
  #Relationship between occurrence and detection 
  #(more common species are [generally] more likely to be detected)
  rho ~ dunif(-1, 1) 
  tau.eta <- tau.lp/(1 - rho^2)
  
})

bivariate_combos <- list(
  Area_Shape = list(y = Ymat, nSpObs = nSpObs, nTraps = nTraps, nPatches = nPatches,
                    z = z, Habitat = Trap_Habitat, CoVar1 = covars_std$Area_Ha, CoVar2 = covars_std$Shape),
  Area_Prox = list(y = Ymat, nSpObs = nSpObs, nTraps = nTraps, nPatches = nPatches,
                   z = z, Habitat = Trap_Habitat, CoVar1 = covars_std$Area_Ha, CoVar2 = covars_std$Prox_500m),
  Prox_Shape = list(y = Ymat, nSpObs = nSpObs, nTraps = nTraps, nPatches = nPatches,
                    z = z, Habitat = Trap_Habitat, CoVar1 = covars_std$Prox_500m, CoVar2 = covars_std$Shape)
)

#Run all models
all_models <- c(univariate_inputs, bivariate_combos)

my.cluster <- makeCluster(6)
registerDoSNOW(my.cluster)
ntasks <- length(all_models) #set up progress bar
pb <- txtProgressBar(max = ntasks, style = 3)
progress <- function(x) setTxtProgressBar(pb, x)
opts <- list(progress = progress) 

results2 <- foreach(i = 1:length(all_models),
                   .options.snow = opts,
                   .packages = 'nimble'
) %dopar% {
  
  if(i %in% 1:3){
    
    model <- nimbleModel(univariate_model, constants = all_models[[i]])
    MCMCconf <- configureMCMC(model, enableWAIC = T, 
                              monitors = c('lpsi', 'lp', 'aHabitat', 'bCoVar'))
    MCMC <- buildMCMC(MCMCconf)
    compModel <- compileNimble(model)
    compMCMC <- compileNimble(MCMC, project = compModel)
    results <- runMCMC(compMCMC, niter = 60000, nburnin = 10000, 
                       thin = 10, nchains = 4, WAIC = T, samplesAsCodaMCMC = T)
    return(results)
    
  } else {
    
    model <- nimbleModel(bivariate_model, constants = all_models[[i]])
    MCMCconf <- configureMCMC(model, enableWAIC = T, 
                              monitors = c('lpsi', 'lp', 'aHabitat', 'bCoVar'))
    MCMC <- buildMCMC(MCMCconf)
    compModel <- compileNimble(model)
    compMCMC <- compileNimble(MCMC, project = compModel)
    results <- runMCMC(compMCMC, niter = 60000, nburnin = 10000, 
                       thin = 10, nchains = 4, WAIC = T, samplesAsCodaMCMC = T)
    return(results)
    
  }
  
}

stopCluster(my.cluster)

names(results2) <- c('Area', 'Prox', 'Shape', 'Area_Shape', 'Area_Prox', 'Prox_Shape')
WAIC_results2 <- data.frame(Variable = names(results2), WAIC = NA)
WAIC_results2$WAIC <- as.numeric(lapply(results2, function(x) x[['WAIC']][['WAIC']]))
full_model_waic$Variable <- 'Full'
WAIC_results2 <- rbind(WAIC_results2, full_model_waic)

save(WAIC_results, results, WAIC_results2, results2, file = 'AF_DB_HMSOM_WAIC_Results.RData')
load('AF_DB_HMSOM_WAIC_Results.RData')

#The below code takes too long to run and I am not aiming to predict, which is where LOO is most needed
#-----Run Model Selection Using Leave-One-Out Cross Validation-----

#0-1 loss function
score01 <- function(preds, bool){
  sum(abs(preds - bool), na.rm = T)
}

#Area under ROC curve
auroc <- function(probs, bool) {
  ok <- is.finite(bool) & is.finite(probs)  # remove NAs
  bool <- bool[ok]
  probs <- probs[ok]
  (mean(rank(probs)[bool==1]) - (sum(bool)+1)/2) / sum(bool==0)
}

# Brier score
brier <- function(probs, bool) {
  sum(bool * (1 - probs)^2 + (1 - bool) * probs^2, na.rm=TRUE)
}

# Log(likelihood) score for Bernoulli trials
loglikBern <- function(probs, bool) {
  sum(dbinom(bool, 1, probs, log=TRUE), na.rm=TRUE)
}

#Define NIMBLE model to generate LOO parameters
loo_model <- function(x){
  
  model <- nimbleModel(candidate_model, constants = x)
  MCMCconf <- configureMCMC(model, 
                            monitors = c('lpsi', 'lp', 'bCoVar', 'aHabitat'))
  MCMC <- buildMCMC(MCMCconf)
  compModel <- compileNimble(model)
  compMCMC <- compileNimble(MCMC, project = compModel)
  results <- runMCMC(compMCMC, niter = 60000, nburnin = 10000, 
                     thin = 20, nchains = 4, samplesAsCodaMCMC = T)
  return(results)
  
}

#Parallel processing 
my.cluster <- makeCluster(detectCores()-1)
registerDoSNOW(my.cluster)

for(i in 1:ncol(prox_buffers)){
  
  cat(paste0('\n', i, '\n'))
  
  #Set up folds
  nFolds <- nPatches
  foldID <- 1:nPatches
  
  #Covariates
  area <- covars$Area_Ha
  shape <- covars$Shape
  Prox <- Prox_buffers[,i]
  buffer_size <- colnames(Prox_buffers)[i]
  
  CVscores <- foreach(fold = 1:nFolds,
                      .packages = c('nimble', 'wiqid', 'MCMCvis'),
                      .inorder = F,
                      .combine = 'rbind'
                      
  ) %dopar% {
    
    #Get position (row index) of test and train sites
    test <- foldID == fold
    train <- foldID != fold
    
    #Re-standardize covars only using training data to avoid data leakage
    trainArea <- wiqid::standardize(area[train]) 
    trainShape <- wiqid::standardize(shape[train])
    trainProx <- wiqid::standardize(Prox[train])
    
    #Standardize test CoVar values to match standardized training data
    testArea <- wiqid::standardize2match(area[test], trainArea)
    testShape <- wiqid::standardize2match(shape[test], trainShape)
    testProx <- wiqid::standardize2match(Prox[test], trainProx)
    
    #Create NIMBLE data for just training sites
    data <- list(y = Ymat[train,,], nTraps = nTraps[train], 
                 nPatches = (nPatches-1), z = z[train,], nSpObs = nSpObs,
                 Habitat = Trap_Habitat[train,], 
                 Area = trainArea, Shape = trainShape, Prox = trainProx)
    
    #Apply model to training data
    out <- loo_model(data)
    
    #Manipulate results
    MCMCresults <- do.call(rbind.data.frame, out)
    nIter <- nrow(MCMCresults)
    
    #Create empty dfs for each metric in each MCMC iteration
    sc01 <- llk <- AUC <- br <- numeric(nIter)
    psi <- Zpred <- rep(NA, nSpObs)
    detect <- matrix(NA, nrow = nSpObs, ncol = nTraps[test])
    true.detect <- Ypred <- ucpd <- detect
    
    #Extract MCMC chains for intercepts and slopes
    psiB0 <- MCMCpstr(out, params = 'lpsi', type = 'chains')[[1]]
    psiCoV <- MCMCpstr(out, params = 'bCoVar', type = 'chains')[[1]]
    pB0 <- MCMCpstr(out, params = 'lp', type = 'chains')[[1]]
    pCoV <- MCMCpstr(out, params = 'aHabitat', type = 'chains')[[1]]
    
    for(n in 1:nIter){
      
      #Extract parameter estimate from iteration n
      psiB0iter <- psiB0[,n]
      psiCoV1iter <- psiCoV[1,,n]
      psiCoV2iter <- psiCoV[2,,n]
      psiCoV3iter <- psiCoV[3,,n]
      pB0iter <- pB0[,n]
      pCoViter <- pCoV[,,n]
      
      for(i in 1:nSpObs){
        
        #Derive occupancy probability estimate for test site in iteration n
        psi[i] <- plogis(psiB0iter[i] + (psiCoV1iter[i] * testArea) + (psiCoV2iter[i] * testShape) + (psiCoV3iter[i] * testProx))
        Zpred[i] <- rbinom(1, 1, psi[i])
        
        for(k in 1:nTraps[test]){
          
          #Extract detection probability estimate for trap j in test site, in iteration n
          detect[i, k] <- plogis(pB0iter[i] + pCoViter[i, Trap_Habitat[test, k]])
          true.detect[i, k] <- detect[i, k] * Zpred[i]
          
          #Derive probability of detection given presence
          #and unconditional probability of detection
          Ypred[i, k] <- rbinom(1, 1, true.detect[i, k])
          ucpd[i, k] <- detect[i, k] * psi[i]
          
        }
      }
      
      Y <- t(na.omit(t(Ymat[test,,])))
      
      sc01[n] <- score01(Ypred, Y)
      llk[n] <- loglikBern(ucpd, Y)
      AUC[n] <- auroc(ucpd, Y) 
      br[n] <- brier(ucpd, Y)
      
    }
    
    #Calculate mean scores across all iterations
    mn_sc01 <- mean(sc01)
    mn_llk <- -2 * mean(llk)
    mn_AUC <- mean(AUC)
    mn_br <- mean(br)
    
    fold_results <- c(fold, mn_sc01, mn_llk, mn_AUC, mn_br)
    return(fold_results)
    
  }
  
  colnames(CVscores) <- c("fold", "score01", "Log_score", "AUC", "Brier")
  var_results <- c(covar_name, colMeans(CVscores[,2:5]))
  final <- list(Summary = var_results, Fold_Results = CVscores)
  
  save(final, file = paste0('./DB_LOOCV/', covar_name, '_LOOCV_Results.RData'))
  gc()
  
}

stopCluster(my.cluster)

cv_files <- list.files(path = './DB_LOOCV', full.names = T)
cv_results <- as.data.frame(matrix(nrow = length(cv_files), ncol = 5))
colnames(cv_results) <- c('Variable', 'Score_01', 'Log_Score', 'AUC', 'Brier')

for(i in 1:length(cv_files)){
  load(cv_files[[i]])
  cv_results[i,] <- final[["Summary"]]
}

cv_results[,2:5] <- sapply(cv_results[,2:5], as.numeric)
