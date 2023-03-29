rm(list=ls()); gc()

library(pacman)
p_load(tidyverse, nimble, wiqid, parallel, nlist,
       doSNOW, abind, MCMCvis, coda, ragg, tidyquant)

setwd('C:/Users/ciarn/Desktop/PhD/Dung_Beetles/')
set.seed(56)
options(memory.limit = 56000)

load('AF_DB_HMSOM_Data.RData')
dbs_wide <- dbs_wide_incidence

#-----Prepare Data-----
Ymat <- dbs_array
SppBinoms <- colnames(dbs_wide)[26:ncol(dbs_wide)]
nTraps <- as.numeric(dbs_wide$Fragment_Samples)
nSpObs <- dim(Ymat)[2]
nPatches <- dim(Ymat)[1]

Trap_Habitat[Trap_Habitat == 'Core'] <- 1
Trap_Habitat[Trap_Habitat == 'Edge'] <- 2
Trap_Habitat <- apply(Trap_Habitat, 2, as.integer)

z <- (Ymat > 0) * 1
z <- apply(z, 1:2, sum, na.rm = T)
z <- ifelse(z > 0, 1, NA)

logArea <- log(dbs_wide$Area_Ha)
AreaStd <- standardize(logArea)
Shape <- dbs_wide$Shape
ShapeStd <- standardize(Shape)
Prox <- log(dbs_wide$Prox_100m + 1)
ProxStd <- standardize(Prox)

#Check distribution of cattle effected sites
#range(logArea)
#breaks <- seq(0, 12, by = 1)
#par(mfrow=2:1)
#hist(log(dbs_wide$Area_Ha[dbs_wide$Cattle == 1]), breaks = breaks, xlim = c(0, 12),
 #    xlab = "log(Area(Ha))", main = "Cattle")
#hist(log(dbs_wide$Area_Ha[dbs_wide$Cattle == 0]), breaks = breaks, xlim = c(0, 12),
 #    xlab = "log(Area(Ha))", main = "No Cattle") 
#hist(dbs_wide$Cattle) 
#no medium or large sites w/ cattle
#many more non-cattle than cattle sites
# not sufficient data to assess impact over all sites, maybe within smaller sites

#Format species Habitat preferences for model inclusion
sps_FHP <- read.csv('AF_DBs_Fragment-Related_Habitat_Preferences.csv', row.names = 1)
sps_order <- gsub('_', ' ', colnames(dbs_wide)[26:108])
sps_FHP <- filter(sps_FHP, Binomial %in% sps_order)
sps_FHP$Binomial == sps_order
  
#Model inputs
data <- list(y = Ymat, nSpObs = nSpObs, nTraps = nTraps, nPatches = nPatches,
             z = z, Habitat = Trap_Habitat, Area = AreaStd,
             Shape = ShapeStd, Prox = ProxStd, FHP = sps_FHP$FHP, FE = sps_FHP$FE)


#-----Structure Nimble HMSOM model w/ data augmentation-----
dbModel_code <-  nimbleCode({
  
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
    mu.CoVar[x] ~ dnorm(0, 0.1) #Covariates are standardised
    sd.CoVar[x] ~ dunif(0, 5)
    tau.CoVar[x] <- 1/sd.CoVar[x]^2
  }
  
  p.mean ~ dbeta(1, 1)
  mu.lp <- logit(p.mean)
  sd.lp ~ dunif(0, 5)
  tau.lp <- 1/sd.lp^2
  
  mu.Habitat ~ dnorm(0, 0.1)
  sd.Habitat ~ dunif(0, 5)
  tau.Habitat <- 1/sd.Habitat^2
  
  #Relationship between occurrence and detection 
  #(more common species are [generally] more likely to be detected)
  rho ~ dunif(-1, 1) 
  tau.eta <- tau.lp/(1 - rho^2)
  
  
  ###Derived quantities###
  
  #Alpha diversity
  for(j in 1:nPatches){
    
    #Richness
    Rich[j] <- sum(z[j,])
    
    #Community Habitat Preference (weighted by occurrence probability, conditional on presence)
    cFHP[j] <- (sum((psi[j, 1:nSpObs] * z[j,]) * FHP[1:nSpObs]) / sum((psi[j, 1:nSpObs] * z[j,]))) 
    cFHPsd[j] <- sqrt(sum((psi[j, 1:nSpObs] * z[j,]) * (FHP[1:nSpObs] - cFHP[j])^2) / 
                        (sum((psi[j, 1:nSpObs] * z[j,]))-1/sum(psi[j, 1:nSpObs] * z[j,])) * sum(psi[j, 1:nSpObs] * z[j,]))
    
    #Forest Endemic Richness
    RichFE[j] <- sum(z[j,] * FE[1:nSpObs]) 
    
  }
  
  #Pairwise Beta Diversity
  for(j in 1:nPatches){
    for(k in 1:nPatches){
      
      #Total community
      a[j,k] <- sum((z[j,] * z[k,]))
      b[j,k] <- sum((z[j,] > z[k,]))
      c[j,k] <- sum((z[k,] > z[j,]))
      
      SorB[j,k] <- ((b[j,k] + c[j,k]) / (2 * a[j,k] + b[j,k] + c[j,k])) #Sorensen Dissimilarity vs major continuous forest
      Turn[j,k] <- ((min(b[j,k], c[j,k])) / (a[j,k] + min(b[j,k], c[j,k]))) #Simpson dissimilarity
      Nest[j,k] <- (SorB[j,k] - Turn[j,k]) #Baselga's nestedness
      
      #Forest Endemics
      aFE[j,k] <- sum(((z[j,] * FE[1:nSpObs]) * (z[k,] * FE[1:nSpObs])))
      bFE[j,k] <- sum(((z[j,] * FE[1:nSpObs]) > (z[k,] * FE[1:nSpObs])))
      cFE[j,k] <- sum(((z[k,] * FE[1:nSpObs]) > (z[j,] * FE[1:nSpObs])))
      
      SorBFE[j,k] <- ((bFE[j,k] + cFE[j,k]) / (2 * aFE[j,k] + bFE[j,k] + cFE[j,k]))
      TurnFE[j,k] <- ((min(bFE[j,k], cFE[j,k])) / (aFE[j,k] + min(bFE[j,k], cFE[j,k])))
      NestFE[j,k] <- (SorBFE[j,k] - TurnFE[j,k])
      
    }
  }
  
  #Number of patches occupied by each species
  for(i in 1:nSpObs){ nOcc[i] <- sum(z[,i]) }
  
})


#-----Run Model-----

outputs <- c('rho', 'p.mean', 'psi.mean',
             'mu.lp', 'sd.lp', 'mu.lpsi', 'sd.lpsi',
             'mu.CoVar', 'sd.CoVar', 
             'mu.Habitat', 'sd.Habitat',
             'lpsi', 'lp', 'psi', 'z',
             'bCoVar', 'aHabitat', 
             'SorB', 'Rich', 'Turn', 'Nest', 'cFHP', 'cFHPsd',
             'RichFE', 'SorBFE', 'TurnFE', 'NestFE', 'nOcc')

#Parallel processing
ncores <- 4 #4 chains
seeds <- 1:ncores #different inits for each chain
my.cluster <- makeCluster(ncores)
registerDoSNOW(my.cluster)

dbModel_chains <- foreach(i = seeds,
                          .packages = 'nimble',
                          .inorder = F
  ) %dopar% {
  
  set.seed(i)
  dbModel <- nimbleModel(dbModel_code, constants = data)
  MCMCconf <- configureMCMC(dbModel, monitors = outputs)
  MCMC <- buildMCMC(MCMCconf)
  comp_dbModel <- compileNimble(dbModel)
  comp_dbMCMC <- compileNimble(MCMC, project = comp_dbModel)
  MCMCsamples <- runMCMC(comp_dbMCMC, niter = 150000, nburnin = 50000, 
                         thin = 20, nchains = 1, samplesAsCodaMCMC = T)
  
  return(MCMCsamples)
  
                          }

stopCluster(my.cluster)

db_results <- mcmc.list(dbModel_chains)

#Inspect results
MCMCtrace(db_results, excl = 'z')
dbModel_summ <- MCMCsummary(db_results, excl = 'z')
#diagPlot(db_results)
db_results <- collapse_chains(db_results)
save(db_results, dbModel_summ, file = 'AF_DB_HMSOM_MCMC.RData')

#-----Assess Goodness of fit of model for each species-----
load('AF_DB_HMSOM_MCMC.RData'); rm(dbModel_summ)

zOut <- MCMCpstr(db_results, params = 'z', type = 'chains')[[1]]
psiOut <- MCMCpstr(db_results, params = 'psi', type = 'chains')[[1]]
lpOut <- MCMCpstr(db_results, params = 'lp', type = 'chains')[[1]]
aHabitatOut <- MCMCpstr(db_results, params = 'aHabitat', type = 'chains')[[1]]
nIter <- dim(db_results[[1]])[1]

sp_bayes_p <- data.frame(ID = 1:nSpObs, Species = SppBinoms, p = NA)
Ysim <- Tobs <- Tsim <- vector('list', nSpObs)
for(i in 1:nSpObs){
  
  Ysim[[i]] <- array(NA, dim = c(nPatches, max(nTraps), nIter))
  Tobs[[i]] <- array(NA, dim = c(nPatches, max(nTraps), nIter))#
  Tsim[[i]] <- array(NA, dim = c(nPatches, max(nTraps), nIter))
  
}  

gft_logic <- vector('list', nSpObs)

for(i in 1:nSpObs){
  for(iter in 1:nIter){
    for(j in 1:nPatches){
      psi <- psiOut[j, i, iter]
      for(k in 1:nTraps[j]){
        p <- as.numeric(plogis(lpOut[i, iter] + aHabitatOut[i, Trap_Habitat[j, k], iter]))
        Tobs[[i]][j, k, iter] <- sum((sqrt(Ymat[j, i, k]) - sqrt(p * zOut[j, i, iter])) ^2)
        Ysim[[i]][j, k, iter] <- rbinom(1, 1, p * zOut[j, i, iter])
        Tsim[[i]][j, k, iter] <- sum((sqrt(Ysim[[i]][j, k, iter]) - sqrt(p * zOut[j, i, iter])) ^2)
      }
    }
  }
  sp_bayes_p$p[i] <- mean(Tsim[[i]] > Tobs[[i]], na.rm = T)
  gft_logic[[i]] <- Tsim[[i]] > Tobs[[i]]
  cat(i, '\n')
}

save(sp_bayes_p, file = 'AF_DB_HMSOM_Sps_Bayes_Ps.RData')
gc()

obj <- ls()
obj <- obj[obj !='gft_logic']
rm(list = obj)
t <- do.call('rbind', gft_logic)
mean(t, na.rm = T)


#R2
nIter <- dim(db_results[[1]])[1]
observed <- Ymat
intercept_prob <- MCMCpstr(db_results, params = 'lp', type = 'chains')[[1]]
aHabitat <- MCMCpstr(db_results, params = 'aHabitat', type = 'chains')[[1]]
fitted <- array(NA, c(nPatches, nSpObs, max(nTraps), nIter))

for(i in 1:nSpObs){
  for(iter in 1:nIter){
    for(j in 1:nPatches){
      for(k in 1:nTraps[j]){
        fitted[j, i, k, iter] <- as.numeric(plogis(intercept_prob[i, iter] + aHabitat[i, Trap_Habitat[j, k], iter]))
      }
    }
  }
}

fitted_mean <- apply(fitted, c(1,2,3), mean, na.rm =T)
intercept_prob <- plogis(intercept_prob) %>% apply(1, mean, na.rm = T)

loglik_full <- sum(log(fitted_mean) * observed + log(1 - fitted_mean) * (1 - observed), na.rm = TRUE)
loglik_null <- vector('numeric', length = nSpObs)
loglik_full <- vector('numeric', length = nSpObs)

for(i in 1:nSpObs){
  
  loglik_full[i] <- sum(log(fitted_mean[,i,]) * observed[,i,] + log(1 - fitted_mean[,i,]) * (1 - observed[,i,]), na.rm = TRUE)
  loglik_null[i] <- log(intercept_prob[i]) * sum(observed[,i,], na.rm = TRUE) + log(1 - intercept_prob[i]) * sum(1 - observed[,i,], na.rm = TRUE)
  
}

loglik_null <- sum(loglik_null)

1 - loglik_full / loglik_null

