rm(list=ls()); gc()

library(pacman)
p_load(tidyverse, nimble, wiqid, MCMCvis, coda, ragg, tidyquant, BayesFactor, Hmisc,
       reshape2, shades, ggpubr, tidybayes, grid, pBrackets, bayestestR, betapart)

setwd('C:/Users/ciarn/Desktop/PhD/Dung_Beetles/')
set.seed(56)

#Plot HMSOM Results
load('AF_DB_HMSOM_MCMC.RData')
load('AF_DB_HMSOM_Data.RData')
options(scipen = 1000000000)

dbs_wide <- dbs_wide_incidence

nSpObs <- 83
nPatches <- 23
SppBinoms <- gsub('_', ' ', colnames(dbs_wide[,26:ncol(dbs_wide)]))

Ymat <- dbs_array
z <- (Ymat > 0) * 1
z <- apply(z, 1:2, sum, na.rm = T)
z <- ifelse(z > 0, 1, NA)

#-----Ocupancy and Detection Probabilities-----
#Get species occupancy probabilities (intercept)
psi0 <- MCMCpstr(db_results, params = 'lpsi', type = 'chains')[[1]]
psi0 <- plogis(psi0)

psi0.median <- apply(psi0, 1, median)
psi0.CRI <- t(apply(psi0, 1, quantile, probs = c(0.05, 0.95)))

#Detection
p0 <- MCMCpstr(db_results, params = 'lp', type = 'chains')[[1]]
p0 <- plogis(p0)

p0.median <- apply(p0, 1, median)
p0.CRI <- t(apply(p0, 1, quantile, probs = c(0.05, 0.95)))

sp.PsiP <- data.frame(Species = colnames(dbs_wide[,26:ncol(dbs_wide)]),
                      psi0 = psi0.median,
                      psi0.LCI = psi0.CRI[,1],
                      psi0.UCI = psi0.CRI[,2],
                      p0 = p0.median,
                      p0.LCI = p0.CRI[,1],
                      p0.UCI = p0.CRI[,2])
sp.PsiP$Species <- gsub('_', ' ', sp.PsiP$Species)
sp.PsiP$Species <- factor(sp.PsiP$Species, levels = sp.PsiP[order(sp.PsiP$psi0), 'Species'])

#Plot species occupancy and detection probabilities
#Occupancy
mu.lpsi <- MCMCpstr(db_results, params = 'mu.lpsi', type = 'chains')[[1]]
sd.lpsi <- MCMCpstr(db_results, params = 'sd.lpsi', type = 'chains')[[1]]

psi0_plot <- 
  ggplot(sp.PsiP, aes(x = psi0, y = Species, colour = Species)) +
    geom_point(stat = 'identity') +
    geom_errorbarh(aes(xmin = psi0.LCI, xmax = psi0.UCI), size = 0.5, height = 0) +
    theme_bw() +
    xlab('Occurence Probability Intercept Estimate') +
    geom_vline(aes(xintercept = plogis(median(mu.lpsi))), 
               colour = 'black', size = 1, linetype = 'dashed') +
    geom_vline(aes(xintercept = plogis(median(mu.lpsi) + 2 * median(sd.lpsi))), 
               colour = 'black', size = 0.5, linetype = 'dashed') +
    geom_vline(aes(xintercept = plogis(median(mu.lpsi) - 2 * median(sd.lpsi))), 
               colour = 'black', size = 0.5, linetype = 'dashed') +
    theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
          axis.line = element_line(colour = 'black', size = 1),
          legend.position = 'none',
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_text(vjust = -0.5, size = 18),
          axis.title.y = element_blank(),
          axis.ticks.x = element_line(size = 1),
          axis.text.x = element_text(colour = 'black', face = 'bold', size = 14),
          axis.text.y = element_text(angle = 15, size = 11, colour = 'black'),
          panel.grid.minor.y = element_line(size = 0.5))

#Detection
mu.lp <- MCMCpstr(db_results, params = 'mu.lp', type = 'chains')[[1]]
sd.lp <- MCMCpstr(db_results, params = 'sd.lp', type = 'chains')[[1]]

p0_plot <- 
  ggplot(sp.PsiP, aes(x = p0, y = Species, colour = Species)) +
    geom_point(stat = 'identity') +
    geom_errorbarh(aes(xmin = p0.LCI, xmax = p0.UCI), size = 0.5, height = 0) +
    theme_bw() +
    xlab('Detection Probability Intercept Estimate') +
    geom_vline(aes(xintercept = plogis(median(mu.lp))), 
               colour = 'black', size = 1, linetype = 'dashed') +
    geom_vline(aes(xintercept = plogis(median(mu.lp) + 2 * median(sd.lp))), 
               colour = 'black', size = 0.5, linetype = 'dashed') +
    geom_vline(aes(xintercept = plogis(median(mu.lp) - 2 * median(sd.lp))), 
               colour = 'black', size = 0.5, linetype = 'dashed') +
    theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
          axis.line = element_line(colour = 'black', size = 1),
          legend.position = 'none',
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_text(vjust = -0.5, size = 18),
          axis.title.y = element_blank(),
          axis.ticks.x = element_line(size = 1),
          axis.text.x = element_text(colour = 'black', face = 'bold', size = 14),
          axis.text.y = element_text(angle = 15, size = 11, colour = 'black'),
          panel.grid.minor.y = element_line(size = 0.5))

p0psi0file <- fs::path('./Figures/',  "Occ_and_Detec_Probs_Plot.png")
agg_png(p0psi0file, width = 50, height = 30, units = "cm", res = 300) 

ggarrange(psi0_plot, p0_plot, 
          ncol = 2, labels = c('A)', 'B)'))

invisible(dev.off())
knitr::include_graphics(p0psi0file)

#Detection core vs edge
core_detect <- MCMCpstr(db_results, params = 'lp', type = 'chains')[[1]]
edge_effect <- MCMCpstr(db_results, params = 'aHabitat', type = 'chains')[[1]]

core_detect <- apply(core_detect, 1, median)
edge_effect <- apply(edge_effect, 1, mean)
edge_detect <- plogis(core_detect + edge_effect)
core_detect <- plogis(core_detect)
p0_change <- edge_detect-core_detect
mean(p0_change)
range(p0_change)
sd(p0_change)

habitat_detection <- as.data.frame(cbind(core_detect, edge_detect)) %>%
  mutate(Species = sub('_', ' ', SppBinoms)) %>%
  pivot_longer(cols = c(core_detect, edge_detect),
               names_to = 'Habitat', 
               values_to = 'Detection')

sum(edge_effect > 0)
sum(edge_effect < 0)

edgefile <- fs::path('./Figures/',  "Edge_Effects_Detection_Plot.png")
agg_png(edgefile, width = 10, height = 15, units = "cm", res = 300) 

ggplot(habitat_detection, aes(x = Habitat, y = Detection, group = Species, colour = Species)) +
  geom_point(stat = 'identity') +
  geom_line() +
  theme_bw() +
  scale_x_discrete(labels = c('Core', 'Edge')) +
  scale_y_continuous(breaks = seq(0,1,.2), labels = seq(0,1,.2), limits = c(0,1)) +
  ylab('Species-Level Median Detection Probability') +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_text(face = 'bold', vjust = 3.3),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'none',
        plot.title = element_text(face = 'bold', vjust = 2.5))

invisible(dev.off())
knitr::include_graphics(edgefile)

#-----Alpha Diversity------
######Species-Level Area Effects#####
par(mfrow=c(1,1))

logArea <- log(dbs_wide$Area_Ha)

areaRange <- range(logArea)
logxx <- seq(areaRange[1], areaRange[2], , 101) #generate hypothetical areas within range of real values
Axx <- exp(logxx)
logxxS <- standardize2match(logxx, logArea) #standardise with mean and SD of real area values
#get species specific area slopes
bArea <- MCMCpstr(db_results, params = 'bCoVar', type = 'chains')[[1]][1,,]
bArea_mean <- apply(bArea, 1, mean)
bArea_sd <- apply(bArea, 1, sd)

#Get model predictions for each generated area
predA <-  matrix(NA, 101, nSpObs)
for(i in 1:101){
  for(x in 1:nSpObs){
    predA[i,x] <- plogis(psi0[x] + bArea_mean[x] * logxxS[i])
  }
}
matplot(Axx, predA, type = 'l', log = 'x', lty = 1,
        ylab = 'Median Occupancy Probability', xlab = 'Area (Ha)', xaxt = 'n')
axis(side = 1, at = c(10, 100, 1000, 10000, 100000), labels = c(10, 100, 1000, 10000, 100000))


######Species Specific Responses to Fragmentation Metrics######
bArea_5 <- apply(bArea, 1, quantile, 0.05)
bArea_95 <- apply(bArea, 1, quantile, 0.95)
bProx <- MCMCpstr(db_results, params = 'bCoVar', type = 'chains')[[1]][3,,]
bProx_mean <- apply(bProx, 1, mean)
bProx_5 <- apply(bProx, 1, quantile, 0.05)
bProx_95 <- apply(bProx, 1, quantile, 0.95)
bShape <-  MCMCpstr(db_results, params = 'bCoVar', type = 'chains')[[1]][2,,]
bShape_mean <- apply(bShape, 1, mean)
bShape_5 <- apply(bShape, 1, quantile, 0.05)
bShape_95 <- apply(bShape, 1, quantile, 0.95)

mu.bs <- MCMCpstr(db_results, params = 'mu.CoVar', type = 'chains')[[1]]
sd.bs <- MCMCpstr(db_results, params = 'sd.CoVar', type = 'chains')[[1]]
mu.bs <- apply(mu.bs, 1, mean)
sd.bs <- apply(sd.bs, 1, mean)

sp_coefs <- data.frame(Species = SppBinoms,
                       psi0 = psi0.median,
                       bArea = bArea_mean,
                       bArea_5 = bArea_5,
                       bArea_95 = bArea_95,
                       bProx = bProx_mean,
                       bProx_5 = bProx_5,
                       bProx_95 = bProx_95,
                       bShape = bShape_mean,
                       bShape_5 = bShape_5,
                       bShape_95 = bShape_95)

sp_coefs$Species <- factor(sp_coefs$Species, levels = sp_coefs[order(sp_coefs$psi0), 'Species'])

bArea_plot <- 
  ggplot(sp_coefs, aes(x = bArea, y = Species, colour = Species)) +
  geom_vline(aes(xintercept = 0), 
             colour = 'black', size = 1) +
  geom_point(stat = 'identity') +
  geom_errorbarh(aes(xmin = bArea_5, xmax = bArea_95), size = 0.5, height = 0) +
  theme_bw() +
  scale_x_continuous(breaks = seq(-3,3,1), labels = seq(-3,3,1), limits = c(-3.1,3.1)) +
  xlab('Area Response (\u03b2 Coefficient)') +
  geom_vline(aes(xintercept = mu.bs[1]), 
             colour = 'black', size = 1, linetype = 'dashed') +
  geom_vline(aes(xintercept = mu.bs[1] + 2 * sd.bs[1]), 
             colour = 'black', size = 0.5, linetype = 'dashed') +
  geom_vline(aes(xintercept = mu.bs[1] - 2 * sd.bs[1]), 
             colour = 'black', size = 0.5, linetype = 'dashed') +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black', size = 1),
        legend.position = 'none',
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(vjust = -0.5, size = 18),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(size = 1),
        axis.text.x = element_text(colour = 'black', face = 'bold', size = 14),
        axis.text.y = element_text(angle = 15, size = 11, colour = 'black'),
        panel.grid.minor.y = element_line(size = 0.5))
bProx_plot <- 
  ggplot(sp_coefs, aes(x = bProx, y = Species, colour = Species)) +
  geom_vline(aes(xintercept = 0), 
             colour = 'black', size = 1) +
  geom_point(stat = 'identity') +
  geom_errorbarh(aes(xmin = bProx_5, xmax = bProx_95), size = 0.5, height = 0) +
  theme_bw() +
  scale_x_continuous(breaks = seq(-3,3,1), labels = seq(-3,3,1), limits = c(-3.1,3.1)) +
  xlab('Mean Habitat Amount Index Response (\u03b2 Coefficient)') +
  geom_vline(aes(xintercept = mu.bs[3]), 
             colour = 'black', size = 1, linetype = 'dashed') +
  geom_vline(aes(xintercept = mu.bs[3] + 2 * sd.bs[3]), 
             colour = 'black', size = 0.5, linetype = 'dashed') +
  geom_vline(aes(xintercept = mu.bs[3] - 2 * sd.bs[3]), 
             colour = 'black', size = 0.5, linetype = 'dashed') +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black', size = 1),
        legend.position = 'none',
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(vjust = -0.5, size = 18),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(size = 1),
        axis.text.x = element_text(colour = 'black', face = 'bold', size = 14),
        axis.text.y = element_text(angle = 15, size = 11, colour = 'black'),
        panel.grid.minor.y = element_line(size = 0.5))
bShape_plot <- 
  ggplot(sp_coefs, aes(x = bShape, y = Species, colour = Species)) +
  geom_vline(aes(xintercept = 0), 
             colour = 'black', size = 1) +
  geom_point(stat = 'identity') +
  geom_errorbarh(aes(xmin = bShape_5, xmax = bShape_95), size = 0.5, height = 0) +
  theme_bw() +
  scale_x_continuous(breaks = seq(-3,3,1), labels = seq(-3,3,1), limits = c(-3.1,3.1)) +
  xlab('Mean Shape Response (\u03b2 Coefficient)') +
  geom_vline(aes(xintercept = mu.bs[2]), 
             colour = 'black', size = 1, linetype = 'dashed') +
  geom_vline(aes(xintercept = mu.bs[2] + 2 * sd.bs[2]), 
             colour = 'black', size = 0.5, linetype = 'dashed') +
  geom_vline(aes(xintercept = mu.bs[2] - 2 * sd.bs[2]), 
             colour = 'black', size = 0.5, linetype = 'dashed') +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black', size = 1),
        legend.position = 'none',
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(vjust = -0.5, size = 18),
        axis.title.y = element_blank(),
        axis.ticks.x = element_line(size = 1),
        axis.text.x = element_text(colour = 'black', face = 'bold', size = 14),
        axis.text.y = element_text(angle = 15, size = 11, colour = 'black'),
        panel.grid.minor.y = element_line(size = 0.5))

spsBfile <- fs::path('./Figures/',  "Sps_Metric_Responses_Plot.png")
agg_png(spsBfile, width = 65, height = 30, units = "cm", res = 300) 

ggarrange(bArea_plot, bProx_plot, bShape_plot, 
          ncol = 3, labels = c('A)', 'B)'))

invisible(dev.off())
knitr::include_graphics(p0psi0file)

#Do coefs overlap 0?
sp_coefs$Area_Sig <- as.integer(mapply(FUN = dplyr::between, 
                                       x = 0, 
                                       left = sp_coefs$bArea_5, 
                                       right = sp_coefs$bArea_95))
sp_coefs$Prox_Sig <- as.integer(mapply(FUN = dplyr::between, 
                                       x = 0, 
                                       left = sp_coefs$bProx_5, 
                                       right = sp_coefs$bProx_95))
sp_coefs$Shape_Sig <- as.integer(mapply(FUN = dplyr::between, 
                                       x = 0, 
                                       left = sp_coefs$bShape_5, 
                                       right = sp_coefs$bShape_95))

for(i in 1:nSpObs){
  sp_coefs$Area_ROPE[i] <- rope(bArea[i,], ci = 0.9)[1,4]
  sp_coefs$Prox_ROPE[i] <- rope(bProx[i,], ci = 0.9)[1,4]
  sp_coefs$Shape_ROPE[i] <- rope(bShape[i,], ci = 0.9)[1,4]
}
######Plot area responses vs habitat preferences#####
sps_FHP <- read.csv('AF_DBs_Fragment-Related_Habitat_Preferences.csv', row.names = 1)
sps_FHP <- filter(sps_FHP, Binomial %in% SppBinoms)

#Plot occ prob x area w/ lines colour coded by habitat preference
colnames(predA) <- SppBinoms
predA <- as.data.frame(cbind(Axx, predA))
colnames(predA)[1] <- 'Area'
predAlong <- pivot_longer(predA, 
                          cols = 2:ncol(predA),
                          names_to = 'Binomial',
                          values_to = 'Occ_Prob')
predAlong <- merge(predAlong, sps_FHP, by = 'Binomial')

Area_Line_Plot <- 
  ggplot(predAlong, aes(x = log(Area), y = Occ_Prob, group = Binomial, colour = FHP^3)) +
    geom_line() +
    saturation(scale_color_distiller(palette = 'RdYlGn', 
                                     direction = 1, 
                                     name = 'FHP',
                                     limits = c(0, 1),
                                     breaks = seq(0,1,.25),
                                     labels = seq(0,1,.25)), 
               delta(0.5)) +
    scale_y_continuous(breaks = seq(0,1,.1), 
                       limits = c(0,1),
                       expand = c(0,0)) +
    scale_x_continuous(limits = c(0, log(140000)),
                       breaks = c(0, log(c(10, 100, 1000, 10000, 100000))),
                       labels = c(0, 10, 100, 1000, 10000, 100000),
                       expand = c(0,0)) +
    expand_limits(x = 0, y = 0) +
    xlab('Patch Area (Ha)') +
    ylab('Occurence Probability') +
    theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
          axis.line = element_line(colour = 'black'),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_text(face = 'bold', vjust = -1.5),
          axis.title.y = element_text(face = 'bold', vjust = 3.25),
          axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
          axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
          axis.line.x = element_line(size = 0.75),
          axis.line.y = element_line(size = 0.75),
          legend.title = element_text(face = 'bold'),
          legend.text = element_text(face = 'bold'),
          legend.key.height = unit(1.5, 'cm'))
  


#Unstandardize coefficient
bArea <- MCMCpstr(db_results, params = 'bCoVar', type = 'chains')[[1]][1,,]
bProx <- MCMCpstr(db_results, params = 'bCoVar', type = 'chains')[[1]][3,,]
bShape <- MCMCpstr(db_results, params = 'bCoVar', type = 'chains')[[1]][2,,]
psi0 <- plogis(MCMCpstr(db_results, params = 'lpsi', type = 'chains')[[1]])

bArea <- bArea/sd(log(dbs_wide$Area_Ha))
bArea_mean <- apply(bArea, 1, mean)
bArea_sd <- apply(bArea, 1, sd)
bProx <- bProx/sd(log(dbs_wide$Prox_100m + 1))
bProx_mean <- apply(bProx, 1, mean)
bProx_sd <- apply(bProx, 1, sd)
bShape <- bShape/sd(dbs_wide$Shape)
bShape_mean <- apply(bShape, 1, mean)
bShape_sd <- apply(bShape, 1, sd)
psi0_mean <- apply(psi0, 1, mean)
psi0_sd <- apply(psi0, 1, sd)

FHMod_consts <- list(FHP = sps_FHP$FHP,
                     n = nSpObs)
AreaFH_data <- list(x = bArea_mean,
                    sd = bArea_sd)
ProxFH_data <- list(x = bProx_mean,
                    sd = bProx_sd)
ShapeFH_data <- list(x = bShape_mean,
                     sd = bShape_sd)
psi0FH_data <- list(x = psi0_mean,
                    sd = psi0_sd)

SpFH_model <- nimbleCode({
  
  #Priors
  a ~ dnorm(0, 0.0001) #Intercept Prior
  b ~ dnorm(0, 0.0001) #Slope prior
  sd.patch ~ dunif(0, 10) #Error in residuals
  tau.patch <- pow(sd.patch, -2) #Accuracy of residuals
  
  #Likelihood
  for(p in 1:n){
    tau.bsd[p] <- pow(sd[p], -2) #Known residual component - aka measurement error
    eps[p] ~ dnorm(0, tau.patch) #Standard residual component 
    muB[p] <- a + FHP[p] * b + eps[p] #Add residual error to model 
    x[p] ~ dnorm(muB[p], tau.bsd[p]) #Measurement error for area response estimates
  }
  
})

outs <- c('a', 'b') #Monitors
inits <- list(a = rnorm(1), b = rnorm(1)) #Initial values

#Run Model
model <- nimbleModel(SpFH_model, constants = FHMod_consts,
                     data = AreaFH_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

AreaFH_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                      thin = 10, nchains = 4, samplesAsCodaMCMC = T) 
compModel$setData(ProxFH_data)
ProxFH_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                      thin = 10, nchains = 4, samplesAsCodaMCMC = T) 
compModel$setData(ShapeFH_data)
ShapeFH_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                      thin = 10, nchains = 4, samplesAsCodaMCMC = T) 
compModel$setData(psi0FH_data)
psi0FH_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                       thin = 10, nchains = 4, samplesAsCodaMCMC = T) 

AreaFH_chains <- MCMCpstr(AreaFH_out, params = c('a', 'b'), type = 'chains')
AreaFH_chains <- data.frame(a = t(AreaFH_chains[['a']]), b = t(AreaFH_chains[['b']]))
ProxFH_chains <- MCMCpstr(ProxFH_out, params = c('a', 'b'), type = 'chains')
ProxFH_chains <- data.frame(a = t(ProxFH_chains[['a']]), b = t(ProxFH_chains[['b']]))
ShapeFH_chains <- MCMCpstr(ShapeFH_out, params = c('a', 'b'), type = 'chains')
ShapeFH_chains <- data.frame(a = t(ShapeFH_chains[['a']]), b = t(ShapeFH_chains[['b']]))
psi0FH_chains <- MCMCpstr(psi0FH_out, params = c('a', 'b'), type = 'chains')
psi0FH_chains <- data.frame(a = t(psi0FH_chains[['a']]), b = t(psi0FH_chains[['b']]))

AreaFH_lines <- AreaFH_chains[sample(nrow(AreaFH_chains), 250, replace = F),]
ProxFH_lines <- ProxFH_chains[sample(nrow(ProxFH_chains), 250, replace = F),]
ShapeFH_lines <- ShapeFH_chains[sample(nrow(ShapeFH_chains), 250, replace = F),]
psi0FH_lines <- psi0FH_chains[sample(nrow(psi0FH_chains), 250, replace = F),]

AreaFH_summary <- MCMCsummary(AreaFH_out)[,c(1,3,5)]
area_a <- AreaFH_summary[1,1]
area_b <- AreaFH_summary[2,1]
ProxFH_summary <- MCMCsummary(ProxFH_out)[,c(1,3,5)]
prox_a <- ProxFH_summary[1,1]
prox_b <- ProxFH_summary[2,1]
ShapeFH_summary <- MCMCsummary(ShapeFH_out)[,c(1,3,5)]
shape_a <- ShapeFH_summary[1,1]
shape_b <- ShapeFH_summary[2,1]
psi0FH_summary <- MCMCsummary(psi0FH_out)[,c(1,3,5)]
psi0_a <- psi0FH_summary[1,1]
psi0_b <- psi0FH_summary[2,1]

AreaFH_data <- data.frame(Species = SppBinoms,
                          BArea = bArea_mean,
                          BAreaSD = bArea_sd,
                          FHP = sps_FHP$FHP)
ProxFH_data <- data.frame(Species = SppBinoms,
                          BProx = bProx_mean,
                          BProxSD = bProx_sd,
                          FHP = sps_FHP$FHP)
ShapeFH_data <- data.frame(Species = SppBinoms,
                           BShape = bShape_mean,
                           BShapeSD = bShape_sd,
                           FHP = sps_FHP$FHP)
psi0FH_data <- data.frame(Species = SppBinoms,
                           psi0 = psi0_mean,
                           psi0SD = psi0_sd,
                           FHP = sps_FHP$FHP)

Area_FHP_Model_Plot <- 
  ggplot(AreaFH_data, aes(x = FHP, y = BArea)) +
    geom_abline(slope = AreaFH_lines[,2], intercept = AreaFH_lines[,1], size = 0.75, alpha = 0.07, colour = '#F6C65A') +
    geom_abline(slope = area_b, intercept = area_a, colour = '#D89605', size = 1.5, linetype = 'dashed') +
    #geom_errorbar(aes(ymin = (BArea - BAreaSD), ymax = (BArea + BAreaSD)), size = 0.7, width = 0.01, colour = '#D89605') +
    geom_point(color = 'black', size = 1.3) +
    #ggtitle('Relationship Between Species Habitat Preference and Area Response') +
    xlab('Forest Specialism Index') +
    ylab('Species Specific Area Response (\u03b2 Coefficient)') +
    scale_y_continuous(breaks = seq(-0.5, 0.5, 0.1), limits = c(-0.5, 0.5)) +
    scale_x_continuous(breaks = seq(0, 1, 0.1),
                       labels = seq(0, 1, 0.1),
                       limits = c(0,1.02), 
                       expand = c(0,0)) +
    theme_bw() +
    theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
          axis.line = element_line(colour = 'black'),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_text(face = 'bold', vjust = -1.5),
          axis.title.y = element_text(face = 'bold', vjust = 2.75),
          axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
          axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
          axis.line.x = element_line(size = 0.75),
          axis.line.y = element_line(size = 0.75),
          legend.title = element_text(face = 'bold', size = 9.5, hjust = 0.4),
          legend.text = element_text(face = 'bold', size = 8),
          legend.key.size = unit(1.8, 'line'))

spareafile <- fs::path('./Figures/',  "Species_Area_Responses_Plot.png")
agg_png(spareafile, width = 25, height = 15, units = "cm", res = 300) 

ggarrange(Area_FHP_Model_Plot, Area_Line_Plot,
          ncol = 2, labels = c('A)', 'B)'))

invisible(dev.off())
knitr::include_graphics(spareafile)


#Other variables
Prox_FHP_Model_Plot <- 
  ggplot(ProxFH_data, aes(x = FHP, y = BProx)) +
  geom_abline(slope = ProxFH_lines[,2], intercept = ProxFH_lines[,1], size = 0.75, alpha = 0.07, colour = '#F6C65A') +
  geom_abline(slope = prox_b, intercept = prox_a, colour = '#D89605', size = 1.5, linetype = 'dashed') +
  geom_errorbar(aes(ymin = (BProx - BProxSD), ymax = (BProx + BProxSD)), size = 0.5, width = 0.01, colour = '#D89605') +
  geom_point(color = 'black', size = 1.3) +
  #ggtitle('Relationship Between Species Habitat Preference and Prox Response') +
  xlab('Forest Specialism Index') +
  ylab('Species Specific Habitat Amount Response (\u03b2 Coefficient)') +
  scale_y_continuous(breaks = seq(-0.5, 0.5, 0.1), limits = c(-0.5, 0.5)) +
  scale_x_continuous(breaks = seq(0, 1, 0.1),
                     labels = seq(0, 1, 0.1),
                     limits = c(0,1.02), 
                     expand = c(0,0)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.title = element_text(face = 'bold', size = 9.5, hjust = 0.4),
        legend.text = element_text(face = 'bold', size = 8),
        legend.key.size = unit(1.8, 'line'))

Shape_FHP_Model_Plot <- 
  ggplot(ShapeFH_data, aes(x = FHP, y = BShape)) +
  geom_abline(slope = ShapeFH_lines[,2], intercept = ShapeFH_lines[,1], size = 0.75, alpha = 0.07, colour = '#F6C65A') +
  geom_abline(slope = shape_b, intercept = shape_a, colour = '#D89605', size = 1.5, linetype = 'dashed') +
  geom_errorbar(aes(ymin = (BShape - BShapeSD), ymax = (BShape + BShapeSD)), size = 0.5, width = 0.01, colour = '#D89605') +
  geom_point(color = 'black', size = 1.3) +
  #ggtitle('Relationship Between Species Habitat Preference and Shape Response') +
  xlab('Forest Specialism Index') +
  ylab('Species Specific Shape Response (\u03b2 Coefficient)') +
  scale_y_continuous(breaks = seq(-1, 1, 0.2), limits = c(-1, 1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.1),
                     labels = seq(0, 1, 0.1),
                     limits = c(0,1.02), 
                     expand = c(0,0)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.title = element_text(face = 'bold', size = 9.5, hjust = 0.4),
        legend.text = element_text(face = 'bold', size = 8),
        legend.key.size = unit(1.8, 'line'))

psi0_FHP_Model_Plot <- 
  ggplot(psi0FH_data, aes(x = FHP, y = psi0)) +
  geom_abline(slope = psi0FH_lines[,2], intercept = psi0FH_lines[,1], size = 0.75, alpha = 0.07, colour = '#F6C65A') +
  geom_abline(slope = psi0_b, intercept = psi0_a, colour = '#D89605', size = 1.5, linetype = 'dashed') +
  geom_errorbar(aes(ymin = (psi0 - psi0SD), ymax = (psi0 + psi0SD)), size = 0.5, width = 0.01, colour = '#D89605') +
  geom_point(color = 'black', size = 1.3) +
  #ggtitle('Relationship Between Species Habitat Preference and psi0') +
  xlab('Forest Specialism Index') +
  ylab('Species Specific Occurence Model Intercept') +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 1, 0.1),
                     labels = seq(0, 1, 0.1),
                     limits = c(0,1.02), 
                     expand = c(0,0)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.title = element_text(face = 'bold', size = 9.5, hjust = 0.4),
        legend.text = element_text(face = 'bold', size = 8),
        legend.key.size = unit(1.8, 'line'))

FHP_models <- fs::path('./Figures/',  "Other_FHP_Models_Plot.png")
agg_png(FHP_models, width = 30, height = 15, units = "cm", res = 300) 

ggarrange(Prox_FHP_Model_Plot, Shape_FHP_Model_Plot, psi0_FHP_Model_Plot,
          ncol = 3, labels = c('A)', 'B)', 'C)'))

invisible(dev.off())
knitr::include_graphics(FHP_models)

#-----Meta-analytical SAR------
###Use Bayesian meta analysis to propagate uncertainty in HMSOM richness estimates###

#Richness estimates df
patch_richness <- dbs_wide %>% select(Fragment, Area_Ha, Shape, Prox_100m)

patch_chains <- MCMCpstr(db_results, params = 'Rich', type = 'chains')[[1]]
patch_richness$Median <- apply(patch_chains, 1, median)
patch_richness$Mean <- rowMeans(patch_chains)
patch_richness$LCI <- apply(patch_chains, 1, function(x) quantile(x, 0.025))
patch_richness$UCI <- apply(patch_chains, 1, function(x) quantile(x, 0.975))

patch_FE_chains <- MCMCpstr(db_results, params = 'RichFE', type = 'chains')[[1]]
patch_richness$FE_Median <- apply(patch_FE_chains, 1, median)
patch_richness$FE_Mean <- rowMeans(patch_FE_chains)
patch_richness$FE_LCI <- apply(patch_FE_chains, 1, function(x) quantile(x, 0.025))
patch_richness$FE_UCI <- apply(patch_FE_chains, 1, function(x) quantile(x, 0.975))

#Mean estimates with SDs
HMSOM.R.mean <- as.numeric(apply(patch_chains, 1, mean))
HMSOM.R.sd <- as.numeric(apply(patch_chains, 1, sd))

#Prep predictors
area <- log(dbs_wide$Area_Ha)

#Set up data and constant for model
HMSOM_MASAR_data <- list(area = area, n = length(area), N = HMSOM.R.mean, psd = HMSOM.R.sd)

#Model structure
SAR_Model <- nimbleCode({
  
  #Priors
  c ~ dnorm(0, 0.0001) #Intercept Prior
  z ~ dnorm(0, 0.0001) #Area effect (slope) prior
  sd.patch ~ dunif(0, 10) #Error in SAR residuals
  tau.patch <- pow(sd.patch, -2) #Accuracy of SAR residuals
  
  #Likelihood
  for(p in 1:n){
    tau.psd[p] <- pow(psd[p], -2) #Known residual component - aka measurement error
    eps.patch[p] ~ dnorm(0, tau.patch) #Standard residual component - that for SAR parameters
    muN[p] <- c + area[p] * z + eps.patch[p] #Add SAR residual error to model 
    N[p] ~ dnorm(muN[p], tau.psd[p]) #Measurement error for richness estimates
  }
  
})

outs <- c('c', 'z') #Monitors
inits <- list(c = rnorm(1), z = rnorm(1)) #Initial values

#Run Model
model <- nimbleModel(SAR_Model, constants = HMSOM_MASAR_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)
HMSOM_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                     thin = 10, nchains = 4, samplesAsCodaMCMC = T) 

MASAR_chains <- MCMCpstr(HMSOM_out, params = c('c', 'z'), type = 'chains')
MASAR_chains <- data.frame(c = t(MASAR_chains[['c']]), z = t(MASAR_chains[['z']]))

MASAR_lines <- MASAR_chains[sample(nrow(MASAR_chains), 250, replace = F),]

HMSOM_MASAR_summary <- MCMCsummary(HMSOM_out)[,c(1,3,5)]
c <- HMSOM_MASAR_summary[1,1]
z <- HMSOM_MASAR_summary[2,1]


MASAR_data <- data.frame(Area = dbs_wide$Area_Ha,
                         Mean = HMSOM.R.mean,
                         SD = HMSOM.R.sd)

Total_SAR_Plot <- 
  ggplot(MASAR_data, aes(x = log(Area), y = Mean)) +
  geom_abline(slope = MASAR_lines[,2], intercept = MASAR_lines[,1], size = 0.75, alpha = 0.07, colour = '#F6C65A') +
  geom_abline(slope = z, intercept = c, colour = '#D89605', size = 1.5, linetype = 'dashed') +
  geom_errorbar(aes(ymin = (Mean - SD), ymax = (Mean + SD)), size = 0.7, width = 0.15, colour = '#D89605') +
  geom_point(color = 'black', size = 1.3) +
  #ggtitle('Total Species Richness') +
  xlab('Area (Ha)') +
  ylab('Estimated Total Species Richness') +
  scale_y_continuous(breaks = seq(0, 80, 10), limits = c(10, 70)) +
  scale_x_continuous(breaks = log(c(1, 10, 100, 1000, 10000, 100000)),
                     labels = c(1, 10, 100, 1000, 10000, 100000)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.title = element_text(face = 'bold', size = 9.5, hjust = 0.4),
        legend.text = element_text(face = 'bold', size = 8),
        legend.key.size = unit(1.8, 'line'))


#Repeat for Forest Endemics

HMSOM.FER.mean <- as.numeric(apply(patch_FE_chains, 1, mean))
HMSOM.FER.sd <- as.numeric(apply(patch_FE_chains, 1, sd))

#Set up data and constant for model
HMSOM_FE_MASAR_data <- list(area = area, n = length(area),
                            N = HMSOM.FER.mean, psd = HMSOM.FER.sd)

#Run Model
FE_model <- nimbleModel(SAR_Model, constants = HMSOM_FE_MASAR_data, inits = inits)
FE_MCMCconf <- configureMCMC(FE_model, monitors = outs)
FE_MCMC <- buildMCMC(FE_MCMCconf)
FE_compModel <- compileNimble(FE_model)
FE_compMCMC <- compileNimble(FE_MCMC, project = FE_compModel)
FE_HMSOM_out <- runMCMC(FE_compMCMC, niter = 50000, nburnin = 10000, 
                        thin = 10, nchains = 4, samplesAsCodaMCMC = T) 


FE_MASAR_chains <- MCMCpstr(FE_HMSOM_out, params = c('c', 'z'), type = 'chains')
FE_MASAR_chains <- data.frame(c = t(FE_MASAR_chains[['c']]), z = t(FE_MASAR_chains[['z']]))

FE_MASAR_lines <- FE_MASAR_chains[sample(nrow(FE_MASAR_chains), 250, replace = F),]

HMSOM_FE_MASAR_summary <- MCMCsummary(FE_HMSOM_out)[,c(1,3,5)]
FE_c <- HMSOM_FE_MASAR_summary[1,1]
FE_z <- HMSOM_FE_MASAR_summary[2,1]


FE_MASAR_data <- data.frame(Area = dbs_wide$Area_Ha,
                            Mean = HMSOM.FER.mean,
                            SD = HMSOM.FER.sd)

FE_SAR_plot <- 
  ggplot(FE_MASAR_data, aes(x = log(Area), y = Mean)) +
  geom_abline(slope = FE_MASAR_lines[,2], intercept = FE_MASAR_lines[,1], size = 0.75, alpha = 0.07, colour = '#58dff5') +
  geom_abline(slope = FE_z, intercept = FE_c, colour = '#05889e', size = 1.5, linetype = 'dashed') +
  geom_errorbar(aes(ymin = (Mean - SD), ymax = (Mean + SD)), size = 0.7, width = 0.15, colour = '#05889e') +
  geom_point(color = 'black', size = 1.3) +
  #ggtitle('Forest Endemic Species Richness') +
  xlab('Area (Ha)') +
  ylab('Estimated Forest Specialist Species Richness') +
  scale_y_continuous(breaks = seq(0, 80, 10), limits = c(10, 70)) +
  scale_x_continuous(breaks = log(c(1, 10, 100, 1000, 10000, 100000)),
                     labels = c(1, 10, 100, 1000, 10000, 100000)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.title = element_text(face = 'bold', size = 9.5, hjust = 0.4),
        legend.text = element_text(face = 'bold', size = 8),
        legend.key.size = unit(1.8, 'line'))

#Compare slopes of FEs vs total
comp_z <- data.frame(Total = MASAR_chains$z,
                     FE = FE_MASAR_chains$z)
comp_z <- pivot_longer(comp_z, 
                       cols = 1:2,
                       names_to = 'Group',
                       values_to = 'Z_Est') %>%
  group_by(Group) %>%
  mutate(Mean = mean(Z_Est),
         UCI = quantile(Z_Est, 0.975),
         LCI = quantile(Z_Est, 0.025))

Slope_Comparison <- 
  ggplot(comp_z, aes(x = Z_Est, fill = Group, color = Group)) +
  geom_density(alpha = 0.3, position = 'identity', size = 1.2) +
  scale_colour_manual(values = c('#58dff5', '#F6C65A'), labels = c('Forest Specialist', 'Overall'), name = 'Community') +
  scale_fill_manual(values = c('#58dff5', '#F6C65A'), , labels = c('Forest Specialist', 'Overall'), name = 'Community') +
  geom_segment(aes(x = mean(MASAR_chains$z), xend = mean(MASAR_chains$z),
                   y = -0.1, yend = -0.025), 
               size = 2, colour = '#F6C65A') +
  geom_segment(aes(x = quantile(MASAR_chains$z, 0.975), xend = quantile(MASAR_chains$z, 0.975),
                   y = -0.1, yend = -0.025), 
               linetype = 'dotted', size = 1.2, colour = '#F6C65A') +
  geom_segment(aes(x = quantile(MASAR_chains$z, 0.025), xend = quantile(MASAR_chains$z, 0.025),
                   y = -0.1, yend = -0.025), 
               linetype = 'dotted', size = 1.2, colour = '#F6C65A') +
  geom_segment(aes(x = mean(FE_MASAR_chains$z), xend = mean(FE_MASAR_chains$z),
                   y = -0.1, yend = -0.025), 
               size = 2, colour = '#58dff5') +
  geom_segment(aes(x = quantile(FE_MASAR_chains$z, 0.975), xend = quantile(FE_MASAR_chains$z, 0.975),
                   y = -0.1, yend = -0.025), 
               linetype = 'dotted', size = 1.2, colour = '#58dff5') +
  geom_segment(aes(x = quantile(FE_MASAR_chains$z, 0.025), xend = quantile(FE_MASAR_chains$z, 0.025),
                   y = -0.1, yend = -0.025), 
               linetype = 'dotted', size = 1.2, colour = '#58dff5') +
  xlab('SAR \u03B2 Estimate') +
  ylab('Density') +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.25,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'bottom',
        legend.title = element_text(face = 'bold'),
        legend.text = element_text(face = 'bold'))


#Combine SAR plots and slope comparison on one page
sarfile <- fs::path('./Figures/',  "Species_Area_Relationships_Plot.png")
agg_png(sarfile, width = 25, height = 20, units = "cm", res = 300) 

ggarrange(ggarrange(Total_SAR_Plot, FE_SAR_plot,
                    labels = c('A)','B)'),
                    ncol = 2),
          ggarrange(NULL, Slope_Comparison, NULL, 
                    labels = c('', 'C)', ''),
                    ncol = 3),
          nrow = 2,
          heights = c(3,1.2))

invisible(dev.off())
knitr::include_graphics(sarfile)

#Test for difference in slopes between FE and overall community using interaction terms
SAR_Comp_Model <- nimbleCode({
  
  #Priors
  for(i in 1:2){a[i] ~ dnorm(0, 0.0001)}
  for(i in 1:2){b[i] ~ dnorm(0, 0.0001)}
  
  sd.patch ~ dunif(0, 10) #Residual error
  tau.patch <- pow(sd.patch, -2) #Accuracy of residuals
  
  #Likelihood
  for(p in 1:n){
    tau.psd[p] <- pow(psd[p], -2) #Known residual component - aka measurement error
    eps.patch[p] ~ dnorm(0, tau.patch) #Standard residual component - that for SAR parameters
    muN[p] <- a[1] + b[1] * area[p] + a[2] * comm[p] + b[2] * (area[p] * comm[p]) + eps.patch[p] #Add SAR residual error to model 
    N[p] ~ dnorm(muN[p], tau.psd[p]) #Measurement error for richness estimates
  }
  
})

outs <- c('a', 'b') #Monitors
inits <- list(a = rnorm(2), b = rnorm(2)) #Initial values

SAR_comp_data <- list(area = c(area,area), n = 2*nPatches, 
                      N = c(HMSOM.R.mean, HMSOM.FER.mean), 
                      psd = c(HMSOM.R.sd, HMSOM.FER.sd),
                      comm = c(rep(0, nPatches), rep(1, nPatches)))

#Run Model
model <- nimbleModel(SAR_Comp_Model, constants = SAR_comp_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)
SAR_comp_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                     thin = 10, nchains = 4, samplesAsCodaMCMC = T) 

SAR_comp_chains <- MCMCpstr(SAR_comp_out, params = c('a', 'b'), type = 'chains')
SAR_comp_summ <- MCMCsummary(SAR_comp_out)

#-----Model richness with all fragmentation metrics as predictors-----

area <- standardize(log(dbs_wide$Area_Ha))
prox <- standardize(log(dbs_wide$Prox_100m + 1))
shape <- standardize(dbs_wide$Shape)

RichMV_Data <- list(area = area, prox = prox, shape = shape,
                    n = length(area), N = HMSOM.R.mean, psd = HMSOM.R.sd)


#Model structure
RichnessMV_Model <- nimbleCode({
  
  #Priors
  a ~ dnorm(0, 0.0001) #Intercept Prior
  b1 ~ dnorm(0, 0.0001) #Area effect (slope) prior
  b2 ~ dnorm(0, 0.0001)
  b3 ~ dnorm(0, 0.0001)
  sd.patch ~ dunif(0, 10) #Error in SAR residuals
  tau.patch <- pow(sd.patch, -2) #Accuracy of SAR residuals
  
  #Likelihood
  for(p in 1:n){
    tau.psd[p] <- pow(psd[p], -2) #Known residual component - aka measurement error
    eps.patch[p] ~ dnorm(0, tau.patch) #Standard residual component - that for SAR parameters
    muN[p] <- a + b1 * area[p] + b2 * prox[p] + b3 * shape[p] + eps.patch[p] #Add SAR residual error to model 
    N[p] ~ dnorm(muN[p], tau.psd[p]) #Measurement error for richness estimates
  }
  
})

outs <- c('a', 'b1', 'b2', 'b3') #Monitors
inits <- list(a = rnorm(1), b1 = rnorm(1), b2 = rnorm(1), b3 = rnorm(1)) #Initial values

#Run Model
model <- nimbleModel(RichnessMV_Model, constants = RichMV_Data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)
RMV_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                     thin = 10, nchains = 4, samplesAsCodaMCMC = T) 

RMV_chains <- MCMCpstr(RMV_out, params = c('b1', 'b2', 'b3'), type = 'chains')

#Repeat for forest endemics

#Prep data
FERichMV_Data <- list(area = area, prox = prox, shape = shape,
                      n = length(area), N = HMSOM.FER.mean, psd = HMSOM.FER.sd)

#Run Model
model <- nimbleModel(RichnessMV_Model, constants = FERichMV_Data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)
FERMV_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                   thin = 10, nchains = 4, samplesAsCodaMCMC = T) 

FERMV_chains <- MCMCpstr(FERMV_out, params = c('b1', 'b2', 'b3'), type = 'chains')

#Prep results for posterior plots
MVR_df <- data.frame(Total_Area = as.numeric(RMV_chains[[1]]),
                     Total_Prox = as.numeric(RMV_chains[[2]]),
                     Total_Shape = as.numeric(RMV_chains[[3]]),
                     FE_Area = as.numeric(FERMV_chains[[1]]),
                     FE_Prox = as.numeric(FERMV_chains[[2]]),
                     FE_Shape = as.numeric(FERMV_chains[[3]]))
MVR_df <- pivot_longer(MVR_df,
                       cols = 1:6,
                       names_to = 'Variable',
                       values_to = 'Estimate')
MVR_df$Community <- gsub('_.*', '', MVR_df$Variable)
MVR_df$Community <- factor(MVR_df$Community, levels = c('Total', 'FE'))
MVR_df$Variable <- factor(MVR_df$Variable, 
                          levels = rev(c('Total_Area', 'FE_Area', 
                                     'Total_Shape', 'FE_Shape',
                                     'Total_Prox', 'FE_Prox')))

MVRichPlot <- fs::path('./Figures/',  "Multivariate_Richness_Plot.png")
agg_png(MVRichPlot, width = 15, height = 20, units = "cm", res = 300) 

MVR_plot <- 
  ggplot(data = MVR_df, aes(x = Estimate, y = Variable, colour = Community, fill = Community)) +
    geom_vline(xintercept = 0, linetype = 9, colour = 'black', size = .8) +
    stat_halfeye(point_interval = 'mean_qi') +
    scale_fill_manual(values = alpha(c('#F6C65A', '#58dff5'), 0.8), labels = c('Overall', 'Forest Endemics')) +
    scale_color_manual(values = c('#D89605', '#05889e'), labels = c('Overall', 'Forest Endemics')) +
    scale_x_continuous(limits = c(-15, 17), 
                       breaks = seq(-15, 15, 5), 
                       labels = seq(-15, 15, 5)) +
    scale_y_discrete(labels = rev(c('Area', '', 'Shape', '', 'Proximity', ''))) +
    xlab('\u03B2 Estimate') +
    theme_bw() +
    theme(plot.margin = unit(c(.5,.5,.5,.5), 'cm'),
          axis.line = element_line(colour = 'black'),
          panel.border = element_blank(), 
          #panel.grid.major = element_blank(),
          #panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_text(face = 'bold', vjust = -1.5),
          axis.title.y = element_blank(),
          axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
          axis.text.y = element_text(face = 'bold', colour = 'black', 
                                     vjust = 0.4, hjust = 0.5, size = 10),
          axis.line.x = element_line(size = 0.75),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = 'none',
          #legend.title = element_text(face = 'bold', size = 8),
          #legend.text = element_text(size = 7)
          )
MVR_plot
grid.brackets(395, 705, 395, 60, lwd = 2, type = 1, col = 'black')
grid.brackets(395, 1353, 395, 706, lwd = 2, type = 1, col = 'black')
grid.brackets(395, 2000, 395, 1354, lwd = 2, type = 1, col = 'black')

invisible(dev.off())
knitr::include_graphics(MVRichPlot)

#Combine all richness plots
allRplots <- fs::path('./Figures/',  "All_Richness_Plots.png")
agg_png(allRplots, width = 25, height = 25, units = "cm", res = 300) 

ggarrange(Total_SAR_Plot, FE_SAR_plot, Slope_Comparison, MVR_plot,
          labels = c('A)','B)', 'C)', 'D)'),
          ncol = 2, nrow = 2)

invisible(dev.off())
knitr::include_graphics(allRplots)

#Significance of coefficients
MVR_df <- data.frame(Total_Area = as.numeric(RMV_chains[[1]]),
                     Total_Prox = as.numeric(RMV_chains[[2]]),
                     Total_Shape = as.numeric(RMV_chains[[3]]),
                     FE_Area = as.numeric(FERMV_chains[[1]]),
                     FE_Prox = as.numeric(FERMV_chains[[2]]),
                     FE_Shape = as.numeric(FERMV_chains[[3]]))

mean(MVR_df$Total_Area)
quantile(MVR_df$Total_Area, 0.95)
quantile(MVR_df$Total_Area, 0.05)
mean(MVR_df$Total_Prox)
quantile(MVR_df$Total_Prox, 0.95)
quantile(MVR_df$Total_Prox, 0.05)
mean(MVR_df$Total_Shape)
quantile(MVR_df$Total_Shape, 0.95)
quantile(MVR_df$Total_Shape, 0.05)

mean(MVR_df$FE_Area)
quantile(MVR_df$FE_Area, 0.95)
quantile(MVR_df$FE_Area, 0.05)
mean(MVR_df$FE_Prox)
quantile(MVR_df$FE_Prox, 0.95)
quantile(MVR_df$FE_Prox, 0.05)
mean(MVR_df$FE_Shape)
quantile(MVR_df$FE_Shape, 0.95)
quantile(MVR_df$FE_Shape, 0.05)

#Test for difference in coefficients
Comp_MV_Model <- nimbleCode({
  
  #Priors
  for(i in 1:2){a[i] ~ dnorm(0, 0.0001)}
  for(i in 1:6){b[i] ~ dnorm(0, 0.0001)}
  
  sd.patch ~ dunif(0, 10) #Residual error
  tau.patch <- pow(sd.patch, -2) #Accuracy of residuals
  
  #Likelihood
  for(p in 1:n){
    tau.psd[p] <- pow(psd[p], -2) #Known residual component - aka measurement error
    eps.patch[p] ~ dnorm(0, tau.patch) #Standard residual component - that for SAR parameters
    muN[p] <- a[1] + b[1] * area[p] + b[2] * prox[p] + b[3] * shape[p] + 
      a[2] * comm[p] + b[4] * (area[p] * comm[p]) + b[5] * (prox[p] * comm[p]) + b[6] * (shape[p] * comm[p]) + eps.patch[p] #Add SAR residual error to model 
    N[p] ~ dnorm(muN[p], tau.psd[p]) #Measurement error for richness estimates
  }
  
})

outs <- c('a', 'b') #Monitors
inits <- list(a = rnorm(2), b = rnorm(6)) #Initial values

data <- list(area = c(area,area), prox = c(prox,prox), shape = c(shape,shape),
             n = 2*nPatches, 
             N = c(HMSOM.R.mean, HMSOM.FER.mean), psd = c(HMSOM.R.sd, HMSOM.FER.sd),
             comm = c(rep(0, nPatches), rep(1, nPatches)))

#Run Model
model <- nimbleModel(Comp_MV_Model, constants = data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)
Comp_MV_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                     thin = 10, nchains = 4, samplesAsCodaMCMC = T) 

Comp_MV_summ <- MCMCsummary(Comp_MV_out, probs = c(0.05, 0.95), pg0 = T)

#####Community Level Habitat Preference#####

CFHP <- MCMCpstr(db_results, params = 'cFHP', type = 'chains')[[1]]
CFHP_Mean <- apply(CFHP, 1, mean)
CFHP_SD <- apply(CFHP, 1, sd)

#Weighted sd
z <- MCMCpstr(db_results, params = 'z', type = 'chains')[[1]]
cFHP <- MCMCpstr(db_results, params = 'cFHP', type = 'chains')[[1]]
psi <- MCMCpstr(db_results, params = 'psi', type = 'chains')[[1]]

cFHPsd <- matrix(NA, nrow = nPatches, ncol = dim(z)[3])
sp_zpsi <- sp_zcFHP <- z

for(i in 1:dim(z)[3]){
  for(p in 1:dim(z)[1]){
    sp_zcFHP[p,,i] <- z[p,,i] * sps_FHP$FHP 
    sp_zpsi[p,,i] <- z[p,,i] * psi[p,,i]
    cFHPsd[p,i] <-  sqrt(wtd.var(x = sp_zcFHP[p,,i], weights = sp_zpsi[p,,i]))
  }
}

cFHPsd_Mean <- apply(cFHPsd, 1, mean) 
CFHPsd_SD <- apply(cFHPsd, 1, sd)

logit_model <- nimbleCode({
  
  #Priors
  a ~ dnorm(0, .0001) 
  b1 ~ dnorm(0, .0001)
  b2 ~ dnorm(0, .0001)
  b3 ~ dnorm(0, .0001)
  b4 ~ dnorm(0, .0001)
  sd.patch ~ dunif(0, 10) 
  tau.patch <- pow(sd.patch, -2) 
  
  #Likelihood
  for(p in 1:n){
    tau[p] <- pow(sd[p], -2) 
    eps[p] ~ dnorm(0, tau.patch) 
    logit(mu[p]) <- a + b1 * area[p] + b2 * area2[p] + b3 * prox[p] + b4 * shape[p] + eps[p] 
    x[p] ~ dnorm(mu[p], tau[p]) 
  }
  
})

CFHP_Model_Consts <- list(area = standardize(log(dbs_wide$Area_Ha)),
                          area2 = standardize(log(dbs_wide$Area_Ha)^2),
                          prox = standardize(log(dbs_wide$Prox_100m + 1)),
                          shape = standardize(dbs_wide$Shape),
                          n = nPatches)
CFHP_Model_Data <- list(x = CFHP_Mean,
                        sd = CFHP_SD)

CFHPsd_Model_Data <- list(x = cFHPsd_Mean,
                          sd = CFHPsd_SD)

inits <- list(a = rnorm(1), b1 = rnorm(1), b2 = rnorm(1), b3 = rnorm(1), b4 = rnorm(1))
outs <- c('a', 'b1', 'b2', 'b3', 'b4')

model <- nimbleModel(logit_model, constants = CFHP_Model_Consts, 
                     data = CFHP_Model_Data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)
model_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                     thin = 10, nchains = 4, samplesAsCodaMCMC = T) 

compModel$setData(CFHPsd_Model_Data)
sdmodel_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                       thin = 10, nchains = 4, samplesAsCodaMCMC = T)

CFHP_CoEfs <- MCMCpstr(model_out, params = c('a', 'b1', 'b2', 'b3', 'b4'), type = 'chains')
CFHP_CoEfs[[1]] <- CFHP_CoEfs[[1]] - 
  (CFHP_CoEfs[[2]] * mean(log(dbs_wide$Area_Ha))/sd(log(dbs_wide$Area_Ha))) -
    (CFHP_CoEfs[[3]] * mean(log(dbs_wide$Area_Ha)^2)/sd(log(dbs_wide$Area_Ha)^2))
        
CFHP_CoEfs[[2]] <- CFHP_CoEfs[[2]]/sd(log(dbs_wide$Area_Ha))
CFHP_CoEfs[[3]] <- CFHP_CoEfs[[3]]/sd(log(dbs_wide$Area_Ha)^2)

mean_inter <- mean(CFHP_CoEfs[[1]])
mean_slope <- mean(CFHP_CoEfs[[2]])
mean_quad_slope <- mean(CFHP_CoEfs[[3]])

CFHPsd_CoEfs <- MCMCpstr(sdmodel_out, params = c('a', 'b1', 'b2', 'b3', 'b4'), type = 'chains')
CFHPsd_CoEfs[[1]] <- CFHPsd_CoEfs[[1]] - 
  (CFHPsd_CoEfs[[2]] * mean(log(dbs_wide$Area_Ha))/sd(log(dbs_wide$Area_Ha))) -
  (CFHPsd_CoEfs[[3]] * mean(log(dbs_wide$Area_Ha)^2)/sd(log(dbs_wide$Area_Ha)^2))

CFHPsd_CoEfs[[2]] <- CFHPsd_CoEfs[[2]]/sd(log(dbs_wide$Area_Ha))
CFHPsd_CoEfs[[3]] <- CFHPsd_CoEfs[[3]]/sd(log(dbs_wide$Area_Ha)^2)

meansd_inter <- mean(CFHPsd_CoEfs[[1]])
meansd_slope <- mean(CFHPsd_CoEfs[[2]])
meansd_quad_slope <- mean(CFHPsd_CoEfs[[3]])


# #Calculate CFHP based on occurrence (i.e. not weighted by occupancy prob)
# FHPMat <- zMat <- MCMCpstr(db_results, params = 'z', type= 'chains')[[1]]
# FHP <- sps_FHP$FHP
# 
# for(i in 1:dim(zMat)[3]){
#   FHPMat[,,i] <- zMat[,,i] %*% diag(FHP)
# }
# 
# FHPMat[FHPMat == 0] <- NA
# FHPMat_summ <- apply(FHPMat, c(1,3), mean, na.rm = T)
# meanFHP <- apply(FHPMat_summ, 1, mean)

CFHP_Data <- data.frame(Area = log(dbs_wide$Area_Ha), 
                        Mean = CFHP_Mean,
                        SD = CFHP_SD)
CFHP_Lines <- data.frame(Intercept = as.numeric(CFHP_CoEfs[[1]]), 
                         Slope = as.numeric(CFHP_CoEfs[[2]]),
                         Quad_Slope = as.numeric(CFHP_CoEfs[[3]]))
CFHP_Lines <- CFHP_Lines[sample(nrow(CFHP_Lines), 250, replace = F),]

CFHPsd_Data <- data.frame(Area = log(dbs_wide$Area_Ha), 
                        Mean = cFHPsd_Mean,
                        SD = CFHPsd_SD)
CFHPsd_Lines <- data.frame(Intercept = as.numeric(CFHPsd_CoEfs[[1]]), 
                         Slope = as.numeric(CFHPsd_CoEfs[[2]]),
                         Quad_Slope = as.numeric(CFHPsd_CoEfs[[3]]))
CFHPsd_Lines <- CFHPsd_Lines[sample(nrow(CFHPsd_Lines), 250, replace = F),]

x_plot <- seq(0, max(log(dbs_wide$Area_Ha)), 0.05)
ysd_plot <- y_plot <- matrix(NA, length(x_plot), 250)

for(i in 1:250){
  y_plot[,i] <- plogis(CFHP_Lines[i,1] + (CFHP_Lines[i,2] * x_plot) + (CFHP_Lines[i,3] * (x_plot^2)))
  ysd_plot[,i] <- plogis(CFHPsd_Lines[i,1] + (CFHPsd_Lines[i,2] * x_plot) + (CFHPsd_Lines[i,3] * (x_plot^2)))
}

colnames(y_plot) <- 1:250
lines_pred <- as.data.frame(cbind(x_plot, y_plot))
lines_pred <- melt(lines_pred, id.vars = 'x_plot', variable.name = 'line',
                   value.name = 'y_plot')
colnames(ysd_plot) <- 1:250
linessd_pred <- as.data.frame(cbind(x_plot, ysd_plot))
linessd_pred <- melt(linessd_pred, id.vars = 'x_plot', variable.name = 'line',
                   value.name = 'y_plot')

mean_y <- plogis(mean_inter + (mean_slope * x_plot) + (mean_quad_slope * (x_plot^2)))
mean_line <- data.frame(x_plot = x_plot, y_plot = mean_y)
meansd_y <- plogis(meansd_inter + (meansd_slope * x_plot) + (meansd_quad_slope * (x_plot^2)))
meansd_line <- data.frame(x_plot = x_plot, y_plot = meansd_y)

CFHP_Area_Plot <- 
  ggplot(lines_pred) +
    geom_line(aes(x_plot, y_plot, group = line), alpha = 0.3, colour = '#56EF98') +
    geom_line(data = mean_line, aes(x_plot, y_plot), colour = '#29AF63', size = 1.5, linetype = 'dashed') +
    geom_errorbar(data = CFHP_Data, aes(Area, Mean, ymin = (Mean - SD), ymax = (Mean + SD)), 
                  size = 0.7, width = 0.15, colour = '#29AF63') +
    geom_point(data = CFHP_Data, aes(Area, Mean), color = 'black', size = 1.8) +
    xlab('Area (Ha)') +
    ylab('Community Mean Forest Specialism Index') +
    scale_y_continuous(breaks = seq(0.75, 1, 0.05), limits = c(0.77, 0.95)) +
    scale_x_continuous(breaks = log(c(1, 10, 100, 1000, 10000, 100000)),
                       labels = c(1, 10, 100, 1000, 10000, 100000)) +
    theme_bw() +
    theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
          axis.line = element_line(colour = 'black'),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_text(face = 'bold', vjust = -1.5),
          axis.title.y = element_text(face = 'bold', vjust = 2.75),
          axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
          axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
          axis.line.x = element_line(size = 0.75),
          axis.line.y = element_line(size = 0.75),
          legend.position = 'none')

CFHPsd_Area_Plot <- 
  ggplot(linessd_pred) +
  geom_line(aes(x_plot, y_plot, group = line), alpha = 0.3, colour = '#56EF98') +
  geom_line(data = meansd_line, aes(x_plot, y_plot), colour = '#29AF63', size = 1.5, linetype = 'dashed') +
  geom_errorbar(data = CFHPsd_Data, aes(Area, Mean, ymin = (Mean - SD), ymax = (Mean + SD)), 
                size = 0.7, width = 0.15, colour = '#29AF63') +
  geom_point(data = CFHPsd_Data, aes(Area, Mean), color = 'black', size = 1.8) +
  xlab('Area (Ha)') +
  ylab('St. Dev of Community Forest Specialism Index') +
  scale_y_continuous(breaks = seq(0.1, 0.3, 0.05), limits = c(0.1, 0.3)) +
  scale_x_continuous(breaks = log(c(1, 10, 100, 1000, 10000, 100000)),
                     labels = c(1, 10, 100, 1000, 10000, 100000)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'none')

#Plot relative effects of each predictor (posterior distribution of standardized betas)
CFHP_stdCoEfs <- MCMCpstr(model_out, params = c('b1', 'b2', 'b3', 'b4'), type = 'chains')
CFHP_stdCoEfs <- data.frame(Area = as.numeric(CFHP_stdCoEfs[[1]]),
                            Area2 = as.numeric(CFHP_stdCoEfs[[2]]),
                            Prox = as.numeric(CFHP_stdCoEfs[[3]]),
                            Shape = as.numeric(CFHP_stdCoEfs[[4]]))
CFHP_stdCoEfs <- pivot_longer(CFHP_stdCoEfs, 
                              cols = 1:4,
                              names_to = 'Variable', 
                              values_to = 'Estimate')
CFHP_Coefs_Plot <- 
  ggplot(data = CFHP_stdCoEfs, 
         aes(x = Estimate, 
             y = Variable, 
             colour = Variable, fill = Variable)) +
    geom_vline(xintercept = 0, linetype = 9, colour = 'black', size = .8) +
    stat_halfeye(point_interval = 'mean_qi') +
    scale_fill_manual(values = alpha(c('#56EF98', '#56EF98', '#F6C65A', '#58dff5'), 0.7)) +
    scale_colour_manual(values = c('#29AF63', '#29AF63', '#D89605', '#05889e')) +
    theme_bw() +
    xlab('\u03B2 Estimate') +
    scale_y_discrete(labels = c(expression(bold(~Area)), expression(bold(~Area^2)), 
                                expression(bold(~Proximity)), expression(bold(~Shape))),
                     position = 'right') +
    scale_x_continuous(limits = c(-.7, .7), 
                       breaks = seq(-.6, .6, .2), 
                       labels = c(-.6,-.4,-.2,0,.2,.4,.6)) +
    theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
          axis.line = element_line(colour = 'black'),
          panel.border = element_blank(), 
          #panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_text(face = 'bold', vjust = -1.5),
          axis.title.y = element_blank(),
          axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
          axis.text.y = element_text(face = 'bold', colour = 'black', 
                                     hjust = 1, vjust = -1.5, size = 10),
          axis.line.x = element_line(size = 0.75),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = 'none')

CFHPsd_stdCoEfs <- MCMCpstr(sdmodel_out, params = c('b1', 'b2', 'b3', 'b4'), type = 'chains')
CFHPsd_stdCoEfs <- data.frame(Area = as.numeric(CFHPsd_stdCoEfs[[1]]),
                            Area2 = as.numeric(CFHPsd_stdCoEfs[[2]]),
                            Prox = as.numeric(CFHPsd_stdCoEfs[[3]]),
                            Shape = as.numeric(CFHPsd_stdCoEfs[[4]]))
CFHPsd_stdCoEfs <- pivot_longer(CFHPsd_stdCoEfs, 
                              cols = 1:4,
                              names_to = 'Variable', 
                              values_to = 'Estimate')


CFHPsd_Coefs_Plot <- 
  ggplot(data = CFHPsd_stdCoEfs, 
         aes(x = Estimate, 
             y = Variable, 
             colour = Variable, fill = Variable)) +
  geom_vline(xintercept = 0, linetype = 9, colour = 'black', size = .8) +
  stat_halfeye(point_interval = 'mean_qi') +
  scale_fill_manual(values = alpha(c('#56EF98', '#56EF98', '#F6C65A', '#58dff5'), 0.7)) +
  scale_colour_manual(values = c('#29AF63', '#29AF63', '#D89605', '#05889e')) +
  theme_bw() +
  xlab('\u03B2 Estimate') +
  scale_y_discrete(labels = c(expression(bold(~Area)), expression(bold(~Area^2)), 
                              expression(bold(~Proximity)), expression(bold(~Shape))),
                   position = 'right') +
  scale_x_continuous(limits = c(-.7, .7), 
                     breaks = seq(-.6, .6, .2), 
                     labels = c(-.6,-.4,-.2,0,.2,.4,.6)) +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_blank(),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', 
                                   hjust = 1, vjust = -1.5, size = 10),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = 'none')


#Combine
CFHPplot <- fs::path('./Figures/',  "Community_Habitat_Preference_Plot.png")
agg_png(CFHPplot, width = 25, height = 20, units = "cm", res = 300) 

ggarrange(CFHP_Area_Plot, CFHP_Coefs_Plot,
          ncol = 2, widths = c(2.1,1),
          labels = c('A)', 'B)'))

invisible(dev.off())
knitr::include_graphics(CFHPplot)

CFHPsdplot <- fs::path('./Figures/',  "SD_of_Community_Habitat_Preference_Plot.png")
agg_png(CFHPsdplot, width = 25, height = 20, units = "cm", res = 300) 

ggarrange(CFHPsd_Area_Plot, CFHPsd_Coefs_Plot,
          ncol = 2, widths = c(2.1,1),
          labels = c('A)', 'B)'))

invisible(dev.off())
knitr::include_graphics(CFHPsdplot)

CFHP_mean_summ <- MCMCsummary(model_out, probs = c(0.05, 0.95), pg0 = T)
CFHP_sd_summ <- MCMCsummary(sdmodel_out, probs = c(0.05, 0.95), pg0 = T)

#------Beta Diversity-------

######Pairwise######
#Load geographic distance matrix
dist_mat <- read.csv('AF_DB_Patch_Distance_Matrix.csv', row.names = 1)
dist_mat <- dist_mat[!(row.names(dist_mat) %in% 'HFA'), !(colnames(dist_mat) %in% 'HFA')]

#Create environmental distance matrices
env_vars <- dbs_wide %>% dplyr::select(Fragment, Area_Ha, Shape, Prox_100m)
prox_dist <- shape_dist <- area_dist <- as.data.frame(matrix(nrow = nPatches, ncol = nPatches))
colnames(area_dist) <- colnames(shape_dist) <- colnames(prox_dist) <- dbs_wide$Fragment
rownames(area_dist) <- rownames(shape_dist) <- rownames(prox_dist) <- dbs_wide$Fragment

for(i in 1:nPatches){
  for(j in 1:nPatches){
    
    area_dist[i,j] <- abs(env_vars$Area_Ha[i] - env_vars$Area_Ha[j])
    shape_dist[i,j] <- abs(env_vars$Shape[i] - env_vars$Shape[j])
    prox_dist[i,j] <- abs(env_vars$Prox_100m[i] - env_vars$Prox_100m[j])
    
    if(i==j){area_dist[i,j] <- shape_dist[i,j] <- prox_dist[i,j] <- NA}
    
  }
}

colnames(dist_mat) <- sub('X', '', colnames(dist_mat))
dist_mat <- dist_mat[order(match(row.names(dist_mat), row.names(area_dist))),]
dist_mat <- dist_mat[, order(match(colnames(dist_mat), colnames(area_dist))),]

area_dist_std <- standardize(as.matrix(log(area_dist + 1)))
area2_dist_std <- standardize(as.matrix(log(area_dist + 1)^2))
prox_dist_std <- standardize(as.matrix(log(prox_dist + 1)))
shape_dist_std <- standardize(as.matrix(shape_dist))
dist_mat_std <- standardize(as.matrix(dist_mat))

#Extract posteriors for pairwise beta estimates
sor <- MCMCpstr(db_results, params = 'SorB', type = 'chains')[[1]]
nest <- MCMCpstr(db_results, params = 'Nest', type = 'chains')[[1]]
turn <- MCMCpstr(db_results, params = 'Turn', type = 'chains')[[1]]
sorFE <- MCMCpstr(db_results, params = 'SorBFE', type = 'chains')[[1]]
nestFE <- MCMCpstr(db_results, params = 'NestFE', type = 'chains')[[1]]
turnFE <- MCMCpstr(db_results, params = 'TurnFE', type = 'chains')[[1]]

sor_mean <- apply(sor, c(1,2), mean, na.rm = T)
sor_sd <- apply(sor, c(1,2), sd, na.rm = T)
nest_mean <- apply(nest, c(1,2), mean, na.rm = T)
nest_sd <- apply(nest, c(1,2), sd, na.rm = T)
turn_mean <- apply(turn, c(1,2), mean, na.rm = T)
turn_sd <- apply(turn, c(1,2), sd, na.rm = T)
sorFE_mean <- apply(sorFE, c(1,2), mean, na.rm = T)
sorFE_sd <- apply(sorFE, c(1,2), sd, na.rm = T)
nestFE_mean <- apply(nestFE, c(1,2), mean, na.rm = T)
nestFE_sd <- apply(nestFE, c(1,2), sd, na.rm = T)
turnFE_mean <- apply(turnFE, c(1,2), mean, na.rm = T)
turnFE_sd <- apply(turnFE, c(1,2), sd, na.rm = T)

#Convert matrices to tables and remove duplicate pairs
distMat2table <- function(x){
  x <- melt(x, varnames = c('Patch1', 'Patch2'))
  x <- x[x$Patch1 != x$Patch2,]
  x <- x[!duplicated(t(apply(x[c("Patch1", "Patch2")], 1, sort))), ]
}

dist <- distMat2table(dist_mat_std)$value
area <- distMat2table(area_dist_std)$value
area2 <- distMat2table(area2_dist_std)$value
prox <- distMat2table(prox_dist_std)$value
shape <- distMat2table(shape_dist_std)$value
sor_mean <- distMat2table(sor_mean)$value
sor_sd <- distMat2table(sor_sd)$value
nest_mean <- distMat2table(nest_mean)$value
nest_sd <- distMat2table(nest_sd)$value
turn_mean <- distMat2table(turn_mean)$value
turn_sd <- distMat2table(turn_sd)$value
sorFE_mean <- distMat2table(sorFE_mean)$value
sorFE_sd <- distMat2table(sorFE_sd)$value
nestFE_mean <- distMat2table(nestFE_mean)$value
nestFE_sd <- distMat2table(nestFE_sd)$value
turnFE_mean <- distMat2table(turnFE_mean)$value
turnFE_sd <- distMat2table(turnFE_sd)$value

area_non_std <- distMat2table(as.matrix(log(area_dist + 1)))$value
area2_non_std <- distMat2table(as.matrix(log(area_dist + 1)^2))$value
prox_non_std <- distMat2table(as.matrix(log(prox_dist + 1)))$value

#####Linear area effect#####
#Structure model
beta_model <- nimbleCode({
  
  #Priors
  for(v in 1:5){b[v] ~ dnorm(0, 0.0001)}
  sd.patch ~ dunif(0, 10) 
  tau.patch <- pow(sd.patch, -2)
  
  #Likelihood
  for(i in 1:n){
        
        tau[i] <- pow(sd[i], -2) 
        eps[i] ~ dnorm(0, tau.patch) 
        logit(mu[i]) <- b[1] + b[2] * area[i] + b[3] * prox[i] + b[4] * shape[i] + b[5] * dist[i] + eps[i] 
        x[i] ~ dnorm(mu[i], tau[i]) 
        
  }
  
})

#Prep inputs
beta_consts <- list(area = area, prox = prox, shape = shape, dist = dist,
                    n = ((nPatches * (nPatches-1))/2))
sor_data <- list(x = sor_mean, sd = sor_sd)
turn_data <- list(x = turn_mean, sd = turn_sd)
nest_data <- list(x = nest_mean, sd = nest_sd)
sorFE_data <- list(x = sorFE_mean, sd = sorFE_sd)
turnFE_data <- list(x = turnFE_mean, sd = turnFE_sd)
nestFE_data <- list(x = nestFE_mean, sd = nestFE_sd)

inits <- list(b = rnorm(5))
outs <- c('b')

#Run Models
model <- nimbleModel(beta_model, constants = beta_consts, 
                     data = sor_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs, enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

sor_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                   thin = 10, nchains = 4, samplesAsCodaMCMC = T) 
compModel$setData(turn_data)
turn_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                   thin = 10, nchains = 4, samplesAsCodaMCMC = T) 
compModel$setData(nest_data)
nest_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                    thin = 10, nchains = 4, samplesAsCodaMCMC = T)
compModel$setData(sorFE_data)
sorFE_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                   thin = 10, nchains = 4, samplesAsCodaMCMC = T) 
compModel$setData(turnFE_data)
turnFE_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                    thin = 10, nchains = 4, samplesAsCodaMCMC = T) 
compModel$setData(nestFE_data)
nestFE_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                    thin = 10, nchains = 4, samplesAsCodaMCMC = T) 


#####Curvilinear Area Effect#####
#Curvilinear beta models
curvilinear_logit_model <- nimbleCode({
  
  #Priors
  for(v in 1:6){ 
    b[v] ~ dnorm(0, 0.0001) }
  sd.patch ~ dunif(0, 10) 
  tau.patch <- pow(sd.patch, -2) 
  
  #Likelihood
  for(p in 1:n){
    tau[p] <- pow(sd[p], -2) 
    eps[p] ~ dnorm(0, tau.patch) 
    logit(mu[p]) <- b[1] + b[2] * area[p] + b[3] * area2[p] + b[4] * prox[p] + b[5] * shape[p] + b[6] * dist[p] + eps[p] 
    x[p] ~ dnorm(mu[p], tau[p]) 
  }
  
})

#Prep inputs
beta_consts <- list(area = area, area2 = area2,
                    prox = prox, shape = shape, dist = dist,
                    n = ((nPatches * (nPatches-1))/2))
sor_data <- list(x = sor_mean, sd = sor_sd)
turn_data <- list(x = turn_mean, sd = turn_sd)
nest_data <- list(x = nest_mean, sd = nest_sd)
sorFE_data <- list(x = sorFE_mean, sd = sorFE_sd)
turnFE_data <- list(x = turnFE_mean, sd = turnFE_sd)
nestFE_data <- list(x = nestFE_mean, sd = nestFE_sd)

inits <- list(b = rnorm(6))
outs <- c('b')

#Run Models
model <- nimbleModel(curvilinear_logit_model, constants = beta_consts, 
                     data = sor_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs, enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

sor_quad_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                   thin = 10, nchains = 4, samplesAsCodaMCMC = T) 
compModel$setData(turn_data)
turn_quad_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                    thin = 10, nchains = 4, samplesAsCodaMCMC = T) 
compModel$setData(nest_data)
nest_quad_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                    thin = 10, nchains = 4, samplesAsCodaMCMC = T)
compModel$setData(sorFE_data)
sorFE_quad_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                     thin = 10, nchains = 4, samplesAsCodaMCMC = T) 
compModel$setData(turnFE_data)
turnFE_quad_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                      thin = 10, nchains = 4, samplesAsCodaMCMC = T) 
compModel$setData(nestFE_data)
nestFE_quad_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                      thin = 10, nchains = 4, samplesAsCodaMCMC = T) 

#Compare WAIC scores
beta_consts <- list(area = area, 
                    prox = prox, shape = shape, dist = dist,
                    n = ((nPatches * (nPatches-1))/2))
sor_data <- list(x = sor_mean, sd = sor_sd)
turn_data <- list(x = turn_mean, sd = turn_sd)
nest_data <- list(x = nest_mean, sd = nest_sd)
sorFE_data <- list(x = sorFE_mean, sd = sorFE_sd)
turnFE_data <- list(x = turnFE_mean, sd = turnFE_sd)
nestFE_data <- list(x = nestFE_mean, sd = nestFE_sd)

inits <- list(b = rnorm(5))
outs <- c('b')

model <- nimbleModel(beta_model, constants = beta_consts, 
                     data = sor_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs, enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

sor_WAIC_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                   thin = 10, nchains = 4, samplesAsCodaMCMC = T, WAIC = T) 

model <- nimbleModel(beta_model, constants = beta_consts, 
                     data = turn_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs, enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

turn_WAIC_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                    thin = 10, nchains = 4, samplesAsCodaMCMC = T, WAIC = T) 

model <- nimbleModel(beta_model, constants = beta_consts, 
                     data = nest_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs, enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

nest_WAIC_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                    thin = 10, nchains = 4, samplesAsCodaMCMC = T, WAIC = T) 

model <- nimbleModel(beta_model, constants = beta_consts, 
                     data = sorFE_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs, enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

sorFE_WAIC_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                     thin = 10, nchains = 4, samplesAsCodaMCMC = T, WAIC = T) 

model <- nimbleModel(beta_model, constants = beta_consts, 
                     data = turnFE_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs, enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

turnFE_WAIC_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                      thin = 10, nchains = 4, samplesAsCodaMCMC = T, WAIC = T) 

model <- nimbleModel(beta_model, constants = beta_consts, 
                     data = nestFE_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs, enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

nestFE_WAIC_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                      thin = 10, nchains = 4, samplesAsCodaMCMC = T, WAIC = T) 



beta_consts <- list(area = area, area2 = area2,
                    prox = prox, shape = shape, dist = dist,
                    n = ((nPatches * (nPatches-1))/2))
sor_data <- list(x = sor_mean, sd = sor_sd)
turn_data <- list(x = turn_mean, sd = turn_sd)
nest_data <- list(x = nest_mean, sd = nest_sd)
sorFE_data <- list(x = sorFE_mean, sd = sorFE_sd)
turnFE_data <- list(x = turnFE_mean, sd = turnFE_sd)
nestFE_data <- list(x = nestFE_mean, sd = nestFE_sd)

inits <- list(b = rnorm(6))
outs <- c('b')


model <- nimbleModel(curvilinear_logit_model, constants = beta_consts, 
                     data = sor_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs, enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

sor_Quad_WAIC_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                   thin = 10, nchains = 4, samplesAsCodaMCMC = T, WAIC = T) 

model <- nimbleModel(curvilinear_logit_model, constants = beta_consts, 
                     data = turn_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs, enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

turn_Quad_WAIC_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                    thin = 10, nchains = 4, samplesAsCodaMCMC = T, WAIC = T) 

model <- nimbleModel(curvilinear_logit_model, constants = beta_consts, 
                     data = nest_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs, enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

nest_Quad_WAIC_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                    thin = 10, nchains = 4, samplesAsCodaMCMC = T, WAIC = T) 

model <- nimbleModel(curvilinear_logit_model, constants = beta_consts, 
                     data = sorFE_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs, enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

sorFE_Quad_WAIC_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                     thin = 10, nchains = 4, samplesAsCodaMCMC = T, WAIC = T) 

model <- nimbleModel(curvilinear_logit_model, constants = beta_consts, 
                     data = turnFE_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs, enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

turnFE_Quad_WAIC_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                      thin = 10, nchains = 4, samplesAsCodaMCMC = T, WAIC = T) 

model <- nimbleModel(curvilinear_logit_model, constants = beta_consts, 
                     data = nestFE_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs, enableWAIC = T)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

nestFE_Quad_WAIC_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                      thin = 10, nchains = 4, samplesAsCodaMCMC = T, WAIC = T)


#Use models including quadratic term for Sorenson and turnover
#####Plot Beta Models#####
#Extract betas
sor_CoEfs <- t(MCMCpstr(sor_out, params = 'b', type = 'chains')[[1]])
nest_CoEfs <- t(MCMCpstr(nest_out, params = 'b', type = 'chains')[[1]])
turn_CoEfs <- t(MCMCpstr(turn_quad_out, params = 'b', type = 'chains')[[1]])
sorFE_CoEfs <- t(MCMCpstr(sorFE_out, params = 'b', type = 'chains')[[1]])
nestFE_CoEfs <- t(MCMCpstr(nestFE_out, params = 'b', type = 'chains')[[1]])
turnFE_CoEfs <- t(MCMCpstr(turnFE_quad_out, params = 'b', type = 'chains')[[1]])

sor_CoEfs[,1] <- sor_CoEfs[,1] - 
  (sor_CoEfs[,2] * mean(area_non_std)/sd(area_non_std))
turn_CoEfs[,1] <- turn_CoEfs[,1] - 
  (turn_CoEfs[,2] * mean(area_non_std)/sd(area_non_std))  -
  (turn_CoEfs[,3] * mean(area2_non_std)/sd(area2_non_std))
nest_CoEfs[,1] <- nest_CoEfs[,1] - 
  (nest_CoEfs[,2] * mean(area_non_std)/sd(area_non_std))
sorFE_CoEfs[,1] <- sorFE_CoEfs[,1] - 
  (sorFE_CoEfs[,2] * mean(area_non_std)/sd(area_non_std))
turnFE_CoEfs[,1] <- turnFE_CoEfs[,1] - 
  (turnFE_CoEfs[,2] * mean(area_non_std)/sd(area_non_std)) -
  (turnFE_CoEfs[,3] * mean(area2_non_std)/sd(area2_non_std))
nestFE_CoEfs[,1] <- nestFE_CoEfs[,1] - 
  (nestFE_CoEfs[,2] * mean(area_non_std)/sd(area_non_std))

sor_CoEfs[,2] <- sor_CoEfs[,2]/sd(area_non_std)
turn_CoEfs[,2] <- turn_CoEfs[,2]/sd(area_non_std)
turn_CoEfs[,3] <- turn_CoEfs[,3]/sd(area2_non_std)
nest_CoEfs[,2] <- nest_CoEfs[,2]/sd(area_non_std)
sorFE_CoEfs[,2] <- sorFE_CoEfs[,2]/sd(area_non_std)
turnFE_CoEfs[,2] <- turnFE_CoEfs[,2]/sd(area_non_std)
turnFE_CoEfs[,3] <- turnFE_CoEfs[,3]/sd(area2_non_std)
nestFE_CoEfs[,2] <- nestFE_CoEfs[,2]/sd(area_non_std)

sor_mean_inter <- mean(sor_CoEfs[,1])
sor_mean_slope <- mean(sor_CoEfs[,2])
turn_mean_inter <- mean(turn_CoEfs[,1])
turn_mean_slope <- mean(turn_CoEfs[,2])
turn_mean_slope2 <- mean(turn_CoEfs[,3])
nest_mean_inter <- mean(nest_CoEfs[,1])
nest_mean_slope <- mean(nest_CoEfs[,2])
sorFE_mean_inter <- mean(sorFE_CoEfs[,1])
sorFE_mean_slope <- mean(sorFE_CoEfs[,2])
turnFE_mean_inter <- mean(turnFE_CoEfs[,1])
turnFE_mean_slope <- mean(turnFE_CoEfs[,2])
turnFE_mean_slope2 <- mean(turnFE_CoEfs[,3])
nestFE_mean_inter <- mean(nestFE_CoEfs[,1])
nestFE_mean_slope <- mean(nestFE_CoEfs[,2])

Sor_Data <- data.frame(Area = area_non_std, 
                       Mean = sor_mean,
                       SD = sor_sd)
Sor_Lines <- data.frame(Intercept = as.numeric(sor_CoEfs[,1]), 
                        Slope = as.numeric(sor_CoEfs[,2]))
Sor_Lines <- Sor_Lines[sample(nrow(Sor_Lines), 250, replace = F),]
Turn_Data <- data.frame(Area = area_non_std, 
                        Mean = turn_mean,
                        SD = turn_sd)
Turn_Lines <- data.frame(Intercept = as.numeric(turn_CoEfs[,1]), 
                         Slope = as.numeric(turn_CoEfs[,2]),
                         Slope2 = as.numeric(turn_CoEfs[,3]))
Turn_Lines <- Turn_Lines[sample(nrow(Turn_Lines), 250, replace = F),]
Nest_Data <- data.frame(Area = area_non_std, 
                        Mean = nest_mean,
                        SD = nest_sd)
Nest_Lines <- data.frame(Intercept = as.numeric(nest_CoEfs[,1]), 
                         Slope = as.numeric(nest_CoEfs[,2]))
Nest_Lines <- Nest_Lines[sample(nrow(Nest_Lines), 250, replace = F),]
SorFE_Data <- data.frame(Area = area_non_std, 
                         Mean = sorFE_mean,
                         SD = sorFE_sd)
SorFE_Lines <- data.frame(Intercept = as.numeric(sorFE_CoEfs[,1]), 
                          Slope = as.numeric(sorFE_CoEfs[,2]))
SorFE_Lines <- SorFE_Lines[sample(nrow(SorFE_Lines), 250, replace = F),]
TurnFE_Data <- data.frame(Area = area_non_std, 
                          Mean = turnFE_mean,
                          SD = turnFE_sd)
TurnFE_Lines <- data.frame(Intercept = as.numeric(turnFE_CoEfs[,1]), 
                           Slope = as.numeric(turnFE_CoEfs[,2]),
                           Slope2 = as.numeric(turnFE_CoEfs[,3]))
TurnFE_Lines <- TurnFE_Lines[sample(nrow(TurnFE_Lines), 250, replace = F),]
NestFE_Data <- data.frame(Area = area_non_std, 
                          Mean = nest_mean,
                          SD = nest_sd)
NestFE_Lines <- data.frame(Intercept = as.numeric(nestFE_CoEfs[,1]), 
                           Slope = as.numeric(nestFE_CoEfs[,2]))
NestFE_Lines <- NestFE_Lines[sample(nrow(NestFE_Lines), 250, replace = F),]

x_plot <- seq(0, max(area_non_std), 0.05)
sor_y_plot <- turn_y_plot <- nest_y_plot <- 
  sorFE_y_plot <- turnFE_y_plot <- nestFE_y_plot <- matrix(NA, length(x_plot), 250)
for(i in 1:250){
  sor_y_plot[,i] <- plogis(Sor_Lines[i,1] + (Sor_Lines[i,2] * x_plot))
  turn_y_plot[,i] <- plogis(Turn_Lines[i,1] + (Turn_Lines[i,2] * x_plot) + (Turn_Lines[i,3] * (x_plot^2)))
  nest_y_plot[,i] <- plogis(Nest_Lines[i,1] + (Nest_Lines[i,2] * x_plot))
  sorFE_y_plot[,i] <- plogis(SorFE_Lines[i,1] + (SorFE_Lines[i,2] * x_plot))
  turnFE_y_plot[,i] <- plogis(TurnFE_Lines[i,1] + (TurnFE_Lines[i,2] * x_plot) + (TurnFE_Lines[i,3] * (x_plot^2)))
  nestFE_y_plot[,i] <- plogis(NestFE_Lines[i,1] + (NestFE_Lines[i,2] * x_plot))
  
}

colnames(sor_y_plot) <- colnames(turn_y_plot) <- colnames(nest_y_plot) <- 
  colnames(sorFE_y_plot) <- colnames(turnFE_y_plot) <- colnames(nestFE_y_plot) <- 1:250

sor_lines_pred <- as.data.frame(cbind(x_plot, sor_y_plot))
sor_lines_pred <- melt(sor_lines_pred, id.vars = 'x_plot', variable.name = 'line',
                       value.name = 'y_plot')
turn_lines_pred <- as.data.frame(cbind(x_plot, turn_y_plot))
turn_lines_pred <- melt(turn_lines_pred, id.vars = 'x_plot', variable.name = 'line',
                        value.name = 'y_plot')
nest_lines_pred <- as.data.frame(cbind(x_plot, nest_y_plot))
nest_lines_pred <- melt(nest_lines_pred, id.vars = 'x_plot', variable.name = 'line',
                        value.name = 'y_plot')
sorFE_lines_pred <- as.data.frame(cbind(x_plot, sorFE_y_plot))
sorFE_lines_pred <- melt(sorFE_lines_pred, id.vars = 'x_plot', variable.name = 'line',
                         value.name = 'y_plot')
turnFE_lines_pred <- as.data.frame(cbind(x_plot, turnFE_y_plot))
turnFE_lines_pred <- melt(turnFE_lines_pred, id.vars = 'x_plot', variable.name = 'line',
                          value.name = 'y_plot')
nestFE_lines_pred <- as.data.frame(cbind(x_plot, nestFE_y_plot))
nestFE_lines_pred <- melt(nestFE_lines_pred, id.vars = 'x_plot', variable.name = 'line',
                          value.name = 'y_plot')

sor_mean_y <- plogis(sor_mean_inter + (sor_mean_slope * x_plot))
sor_mean_line <- data.frame(x_plot = x_plot, y_plot = sor_mean_y)
turn_mean_y <- plogis(turn_mean_inter + (turn_mean_slope * x_plot) + (turn_mean_slope2 * (x_plot^2)))
turn_mean_line <- data.frame(x_plot = x_plot, y_plot = turn_mean_y)
nest_mean_y <- plogis(nest_mean_inter + (nest_mean_slope * x_plot))
nest_mean_line <- data.frame(x_plot = x_plot, y_plot = nest_mean_y)
sorFE_mean_y <- plogis(sorFE_mean_inter + (sorFE_mean_slope * x_plot))
sorFE_mean_line <- data.frame(x_plot = x_plot, y_plot = sorFE_mean_y)
turnFE_mean_y <- plogis(turnFE_mean_inter + (turnFE_mean_slope * x_plot) + (turnFE_mean_slope2 * (x_plot^2))) 
turnFE_mean_line <- data.frame(x_plot = x_plot, y_plot = turnFE_mean_y)
nestFE_mean_y <- plogis(nestFE_mean_inter + (nestFE_mean_slope * x_plot)) 
nestFE_mean_line <- data.frame(x_plot = x_plot, y_plot = nestFE_mean_y)

Sor_Area_Plot <- 
  ggplot(sor_lines_pred) +
  geom_line(aes(x_plot, sor_y_plot, group = line), alpha = 0.3, colour = '#F6C65A') +
  geom_line(data = sor_mean_line, aes(x_plot, y_plot), colour = '#D89605', size = 1.5, linetype = 'dashed') +
  geom_errorbar(data = Sor_Data, aes(Area, Mean, ymin = (Mean - SD), ymax = (Mean + SD)), 
                size = 0.7, width = 0.15, colour = '#D89605') +
  geom_point(data = Sor_Data, aes(Area, Mean), color = 'black', size = 1.3) +
  xlab('Difference in Area (Ha)') +
  ylab('Sorenson Dissimilarity') +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  scale_x_continuous(breaks = log(c(1, 11, 101, 1001, 10001, 100001)),
                     labels = c(0, 10, 100, 1000, 10000, 100000)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'none')
Turn_Area_Plot <- 
  ggplot(turn_lines_pred) +
  geom_line(aes(x_plot, turn_y_plot, group = line), alpha = 0.3, colour = '#F6C65A') +
  geom_line(data = turn_mean_line, aes(x_plot, y_plot), colour = '#D89605', size = 1.5, linetype = 'dashed') +
  geom_errorbar(data = Turn_Data, aes(Area, Mean, ymin = (Mean - SD), ymax = (Mean + SD)), 
                size = 0.7, width = 0.15, colour = '#D89605') +
  geom_point(data = Turn_Data, aes(Area, Mean), color = 'black', size = 1.3) +
  xlab('Difference in Area (Ha)') +
  ylab('Turnover Resultant Dissimilarity') +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  scale_x_continuous(breaks = log(c(1, 11, 101, 1001, 10001, 100001)),
                     labels = c(0, 10, 100, 1000, 10000, 100000)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'none')
Nest_Area_Plot <- 
  ggplot(nest_lines_pred) +
  geom_line(aes(x_plot, nest_y_plot, group = line), alpha = 0.3, colour = '#F6C65A') +
  geom_line(data = nest_mean_line, aes(x_plot, y_plot), colour = '#D89605', size = 1.5, linetype = 'dashed') +
  geom_errorbar(data = Nest_Data, aes(Area, Mean, ymin = (Mean - SD), ymax = (Mean + SD)), 
                size = 0.7, width = 0.15, colour = '#D89605') +
  geom_point(data = Nest_Data, aes(Area, Mean), color = 'black', size = 1.3) +
  xlab('Difference in Area (Ha)') +
  ylab('Nestedness Resultant Dissimilarity') +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  scale_x_continuous(breaks = log(c(1, 11, 101, 1001, 10001, 100001)),
                     labels = c(0, 10, 100, 1000, 10000, 100000)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'none')

SorFE_Area_Plot <- 
  ggplot(sorFE_lines_pred) +
  geom_line(aes(x_plot, sorFE_y_plot, group = line), alpha = 0.3, colour = '#58dff5') +
  geom_line(data = sorFE_mean_line, aes(x_plot, y_plot), colour = '#05889e', size = 1.5, linetype = 'dashed') +
  geom_errorbar(data = SorFE_Data, aes(Area, Mean, ymin = (Mean - SD), ymax = (Mean + SD)), 
                size = 0.7, width = 0.15, colour = '#05889e') +
  geom_point(data = SorFE_Data, aes(Area, Mean), color = 'black', size = 1.3) +
  xlab('Difference in Area (Ha)') +
  ylab('Sorenson Dissimilarity') +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  scale_x_continuous(breaks = log(c(1, 11, 101, 1001, 10001, 100001)),
                     labels = c(0, 10, 100, 1000, 10000, 100000)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'none')
TurnFE_Area_Plot <- 
  ggplot(turnFE_lines_pred) +
  geom_line(aes(x_plot, turnFE_y_plot, group = line), alpha = 0.3, colour = '#58dff5') +
  geom_line(data = turnFE_mean_line, aes(x_plot, y_plot), colour = '#05889e', size = 1.5, linetype = 'dashed') +
  geom_errorbar(data = TurnFE_Data, aes(Area, Mean, ymin = (Mean - SD), ymax = (Mean + SD)), 
                size = 0.7, width = 0.15, colour = '#05889e') +
  geom_point(data = TurnFE_Data, aes(Area, Mean), color = 'black', size = 1.3) +
  xlab('Difference in Area (Ha)') +
  ylab('Turnover Resultant Dissimilarity') +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  scale_x_continuous(breaks = log(c(1, 11, 101, 1001, 10001, 100001)),
                     labels = c(0, 10, 100, 1000, 10000, 100000)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'none')
NestFE_Area_Plot <- 
  ggplot(nestFE_lines_pred) +
  geom_line(aes(x_plot, nestFE_y_plot, group = line), alpha = 0.3, colour = '#58dff5') +
  geom_line(data = nestFE_mean_line, aes(x_plot, y_plot), colour = '#05889e', size = 1.5, linetype = 'dashed') +
  geom_errorbar(data = NestFE_Data, aes(Area, Mean, ymin = (Mean - SD), ymax = (Mean + SD)), 
                size = 0.7, width = 0.15, colour = '#05889e') +
  geom_point(data = NestFE_Data, aes(Area, Mean), color = 'black', size = 1.3) +
  xlab('Difference in Area (Ha)') +
  ylab('Nestedness Resultant Dissimilarity') +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  scale_x_continuous(breaks = log(c(1, 11, 101, 1001, 10001, 100001)),
                     labels = c(0, 10, 100, 1000, 10000, 100000)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'none')

PairwisePlot <- fs::path('./Figures/',  "Pairwise_Dissimilarity_Plot.png")
agg_png(PairwisePlot, width = 20, height = 30, units = "cm", res = 300)

pair_plot <- ggarrange(Sor_Area_Plot, SorFE_Area_Plot, Turn_Area_Plot, TurnFE_Area_Plot, Nest_Area_Plot, NestFE_Area_Plot,
                       nrow = 3, ncol = 2, labels = c('A)', '', 'B)', '','C)', ''))
annotate_figure(pair_plot, fig.lab = '                   Overall Community                                                Forest Specialists',
                fig.lab.pos = 'top.left', fig.lab.face = 'bold', fig.lab.size = 12.5)

invisible(dev.off())
knitr::include_graphics(PairwisePlot)

######Create plots of pairwise difference in sorrounding forest amount#####
#Extract betas#
sor_CoEfs <- t(MCMCpstr(sor_out, params = 'b', type = 'chains')[[1]])
nest_CoEfs <- t(MCMCpstr(nest_out, params = 'b', type = 'chains')[[1]])
turn_CoEfs <- t(MCMCpstr(turn_quad_out, params = 'b', type = 'chains')[[1]])
sorFE_CoEfs <- t(MCMCpstr(sorFE_out, params = 'b', type = 'chains')[[1]])
nestFE_CoEfs <- t(MCMCpstr(nestFE_out, params = 'b', type = 'chains')[[1]])
turnFE_CoEfs <- t(MCMCpstr(turnFE_quad_out, params = 'b', type = 'chains')[[1]])

sor_CoEfs[,1] <- sor_CoEfs[,1] - 
  (sor_CoEfs[,3] * mean(prox_non_std)/sd(prox_non_std)) 
turn_CoEfs[,1] <- turn_CoEfs[,1] - 
  (sor_CoEfs[,4] * mean(prox_non_std)/sd(prox_non_std)) 
nest_CoEfs[,1] <- nest_CoEfs[,1] - 
  (sor_CoEfs[,3] * mean(prox_non_std)/sd(prox_non_std)) 
sorFE_CoEfs[,1] <- sorFE_CoEfs[,1] - 
  (sor_CoEfs[,3] * mean(prox_non_std)/sd(prox_non_std)) 
turnFE_CoEfs[,1] <- turnFE_CoEfs[,1] - 
  (sor_CoEfs[,4] * mean(prox_non_std)/sd(prox_non_std)) 
nestFE_CoEfs[,1] <- nestFE_CoEfs[,1] - 
  (sor_CoEfs[,3] * mean(prox_non_std)/sd(prox_non_std)) 

sor_CoEfs[,3] <- sor_CoEfs[,3]/sd(prox_non_std)
turn_CoEfs[,4] <- turn_CoEfs[,4]/sd(prox_non_std)
nest_CoEfs[,3] <- nest_CoEfs[,3]/sd(prox_non_std)
sorFE_CoEfs[,3] <- sorFE_CoEfs[,3]/sd(prox_non_std)
turnFE_CoEfs[,4] <- turnFE_CoEfs[,4]/sd(prox_non_std)
nestFE_CoEfs[,3] <- nestFE_CoEfs[,3]/sd(prox_non_std)

sor_mean_inter <- mean(sor_CoEfs[,1])
sor_mean_slope <- mean(sor_CoEfs[,3])
turn_mean_inter <- mean(turn_CoEfs[,1])
turn_mean_slope <- mean(turn_CoEfs[,4])
nest_mean_inter <- mean(nest_CoEfs[,1])
nest_mean_slope <- mean(nest_CoEfs[,3])
sorFE_mean_inter <- mean(sorFE_CoEfs[,1])
sorFE_mean_slope <- mean(sorFE_CoEfs[,3])
turnFE_mean_inter <- mean(turnFE_CoEfs[,1])
turnFE_mean_slope <- mean(turnFE_CoEfs[,4])
nestFE_mean_inter <- mean(nestFE_CoEfs[,1])
nestFE_mean_slope <- mean(nestFE_CoEfs[,3])

Sor_Data <- data.frame(Prox = prox_non_std, 
                       Mean = sor_mean,
                       SD = sor_sd)
Sor_Lines <- data.frame(Intercept = as.numeric(sor_CoEfs[,1]), 
                        Slope = as.numeric(sor_CoEfs[,3]))
Sor_Lines <- Sor_Lines[sample(nrow(Sor_Lines), 250, replace = F),]
Turn_Data <- data.frame(Prox = prox_non_std,
                        Mean = turn_mean,
                        SD = turn_sd)
Turn_Lines <- data.frame(Intercept = as.numeric(turn_CoEfs[,1]), 
                         Slope = as.numeric(turn_CoEfs[,4]))
Turn_Lines <- Turn_Lines[sample(nrow(Turn_Lines), 250, replace = F),]
Nest_Data <- data.frame(Prox = prox_non_std, 
                        Mean = nest_mean,
                        SD = nest_sd)
Nest_Lines <- data.frame(Intercept = as.numeric(nest_CoEfs[,1]), 
                         Slope = as.numeric(nest_CoEfs[,3]))
Nest_Lines <- Nest_Lines[sample(nrow(Nest_Lines), 250, replace = F),]
SorFE_Data <- data.frame(Prox = prox_non_std, 
                         Mean = sorFE_mean,
                         SD = sorFE_sd)
SorFE_Lines <- data.frame(Intercept = as.numeric(sorFE_CoEfs[,1]), 
                          Slope = as.numeric(sorFE_CoEfs[,3]))
SorFE_Lines <- SorFE_Lines[sample(nrow(SorFE_Lines), 250, replace = F),]
TurnFE_Data <- data.frame(Prox = prox_non_std, 
                          Mean = turnFE_mean,
                          SD = turnFE_sd)
TurnFE_Lines <- data.frame(Intercept = as.numeric(turnFE_CoEfs[,1]), 
                           Slope = as.numeric(turnFE_CoEfs[,4]))
TurnFE_Lines <- TurnFE_Lines[sample(nrow(TurnFE_Lines), 250, replace = F),]
NestFE_Data <- data.frame(Prox = prox_non_std, 
                          Mean = nest_mean,
                          SD = nest_sd)
NestFE_Lines <- data.frame(Intercept = as.numeric(nestFE_CoEfs[,1]), 
                           Slope = as.numeric(nestFE_CoEfs[,3]))
NestFE_Lines <- NestFE_Lines[sample(nrow(NestFE_Lines), 250, replace = F),]

x_plot <- seq(0, max(prox_non_std), 0.05)
x2_plot <- 
  sor_y_plot <- turn_y_plot <- nest_y_plot <- 
  sorFE_y_plot <- turnFE_y_plot <- nestFE_y_plot <- matrix(NA, length(x_plot), 250)
for(i in 1:250){
  sor_y_plot[,i] <- plogis(Sor_Lines[i,1] + (Sor_Lines[i,2] * x_plot))
  turn_y_plot[,i] <- plogis(Turn_Lines[i,1] + (Turn_Lines[i,2] * x_plot))
  nest_y_plot[,i] <- plogis(Nest_Lines[i,1] + (Nest_Lines[i,2] * x_plot))
  sorFE_y_plot[,i] <- plogis(SorFE_Lines[i,1] + (SorFE_Lines[i,2] * x_plot))
  turnFE_y_plot[,i] <- plogis(TurnFE_Lines[i,1] + (TurnFE_Lines[i,2] * x_plot))
  nestFE_y_plot[,i] <- plogis(NestFE_Lines[i,1] + (NestFE_Lines[i,2] * x_plot))
  
}

colnames(sor_y_plot) <- colnames(turn_y_plot) <- colnames(nest_y_plot) <- 
  colnames(sorFE_y_plot) <- colnames(turnFE_y_plot) <- colnames(nestFE_y_plot) <- 1:250

sor_lines_pred <- as.data.frame(cbind(x_plot, sor_y_plot))
sor_lines_pred <- melt(sor_lines_pred, id.vars = 'x_plot', variable.name = 'line',
                       value.name = 'y_plot')
turn_lines_pred <- as.data.frame(cbind(x_plot, turn_y_plot))
turn_lines_pred <- melt(turn_lines_pred, id.vars = 'x_plot', variable.name = 'line',
                        value.name = 'y_plot')
nest_lines_pred <- as.data.frame(cbind(x_plot, nest_y_plot))
nest_lines_pred <- melt(nest_lines_pred, id.vars = 'x_plot', variable.name = 'line',
                        value.name = 'y_plot')
sorFE_lines_pred <- as.data.frame(cbind(x_plot, sorFE_y_plot))
sorFE_lines_pred <- melt(sorFE_lines_pred, id.vars = 'x_plot', variable.name = 'line',
                         value.name = 'y_plot')
turnFE_lines_pred <- as.data.frame(cbind(x_plot, turnFE_y_plot))
turnFE_lines_pred <- melt(turnFE_lines_pred, id.vars = 'x_plot', variable.name = 'line',
                          value.name = 'y_plot')
nestFE_lines_pred <- as.data.frame(cbind(x_plot, nestFE_y_plot))
nestFE_lines_pred <- melt(nestFE_lines_pred, id.vars = 'x_plot', variable.name = 'line',
                          value.name = 'y_plot')

sor_mean_y <- plogis(sor_mean_inter + (sor_mean_slope * x_plot))
sor_mean_line <- data.frame(x_plot = x_plot, y_plot = sor_mean_y) 
turn_mean_y <- plogis(turn_mean_inter + (turn_mean_slope * x_plot))
turn_mean_line <- data.frame(x_plot = x_plot, y_plot = turn_mean_y)
nest_mean_y <- plogis(nest_mean_inter + (nest_mean_slope * x_plot))
nest_mean_line <- data.frame(x_plot = x_plot, y_plot = nest_mean_y)
sorFE_mean_y <- plogis(sorFE_mean_inter + (sorFE_mean_slope * x_plot))
sorFE_mean_line <- data.frame(x_plot = x_plot, y_plot = sorFE_mean_y)
turnFE_mean_y <- plogis(turnFE_mean_inter + (turnFE_mean_slope * x_plot)) 
turnFE_mean_line <- data.frame(x_plot = x_plot, y_plot = turnFE_mean_y)
nestFE_mean_y <- plogis(nestFE_mean_inter + (nestFE_mean_slope * x_plot)) 
nestFE_mean_line <- data.frame(x_plot = x_plot, y_plot = nestFE_mean_y)

Sor_Prox_Plot <- 
  ggplot(sor_lines_pred) +
  geom_line(aes(x_plot, sor_y_plot, group = line), alpha = 0.3, colour = '#F6C65A') +
  geom_line(data = sor_mean_line, aes(x_plot, y_plot), colour = '#D89605', size = 1.5, linetype = 'dashed') +
  geom_errorbar(data = Sor_Data, aes(Prox, Mean, ymin = (Mean - SD), ymax = (Mean + SD)), 
                size = 0.7, width = 0.15, colour = '#D89605') +
  geom_point(data = Sor_Data, aes(Prox, Mean), color = 'black', size = 1.3) +
  xlab('Difference in Sorrounding Forest Amount (log(x+1))') +
  ylab('Sorenson Dissimilarity') +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'none')
Turn_Prox_Plot <- 
  ggplot(turn_lines_pred) +
  geom_line(aes(x_plot, turn_y_plot, group = line), alpha = 0.3, colour = '#F6C65A') +
  geom_line(data = turn_mean_line, aes(x_plot, y_plot), colour = '#D89605', size = 1.5, linetype = 'dashed') +
  geom_errorbar(data = Turn_Data, aes(Prox, Mean, ymin = (Mean - SD), ymax = (Mean + SD)), 
                size = 0.7, width = 0.15, colour = '#D89605') +
  geom_point(data = Turn_Data, aes(Prox, Mean), color = 'black', size = 1.3) +
  xlab('Difference in Sorrounding Forest Amount (log(x+1))') +
  ylab('Turnover Resultant Dissimilarity') +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'none')
Nest_Prox_Plot <- 
  ggplot(nest_lines_pred) +
  geom_line(aes(x_plot, nest_y_plot, group = line), alpha = 0.3, colour = '#F6C65A') +
  geom_line(data = nest_mean_line, aes(x_plot, y_plot), colour = '#D89605', size = 1.5, linetype = 'dashed') +
  geom_errorbar(data = Nest_Data, aes(Prox, Mean, ymin = (Mean - SD), ymax = (Mean + SD)), 
                size = 0.7, width = 0.15, colour = '#D89605') +
  geom_point(data = Nest_Data, aes(Prox, Mean), color = 'black', size = 1.3) +
  xlab('Difference in Sorrounding Forest Amount (log(x+1))') +
  ylab('Nestedness Resultant Dissimilarity') +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'none')

SorFE_Prox_Plot <- 
  ggplot(sorFE_lines_pred) +
  geom_line(aes(x_plot, sorFE_y_plot, group = line), alpha = 0.3, colour = '#58dff5') +
  geom_line(data = sorFE_mean_line, aes(x_plot, y_plot), colour = '#05889e', size = 1.5, linetype = 'dashed') +
  geom_errorbar(data = SorFE_Data, aes(Prox, Mean, ymin = (Mean - SD), ymax = (Mean + SD)), 
                size = 0.7, width = 0.15, colour = '#05889e') +
  geom_point(data = SorFE_Data, aes(Prox, Mean), color = 'black', size = 1.3) +
  xlab('Difference in Sorrounding Forest Amount (log(x+1))') +
  ylab('Sorenson Dissimilarity') +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'none')
TurnFE_Prox_Plot <- 
  ggplot(turnFE_lines_pred) +
  geom_line(aes(x_plot, turnFE_y_plot, group = line), alpha = 0.3, colour = '#58dff5') +
  geom_line(data = turnFE_mean_line, aes(x_plot, y_plot), colour = '#05889e', size = 1.5, linetype = 'dashed') +
  geom_errorbar(data = TurnFE_Data, aes(Prox, Mean, ymin = (Mean - SD), ymax = (Mean + SD)), 
                size = 0.7, width = 0.15, colour = '#05889e') +
  geom_point(data = TurnFE_Data, aes(Prox, Mean), color = 'black', size = 1.3) +
  xlab('Difference in Sorrounding Forest Amount (log(x+1))') +
  ylab('Turnover Resultant Dissimilarity') +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'none')
NestFE_Prox_Plot <- 
  ggplot(nestFE_lines_pred) +
  geom_line(aes(x_plot, nestFE_y_plot, group = line), alpha = 0.3, colour = '#58dff5') +
  geom_line(data = nestFE_mean_line, aes(x_plot, y_plot), colour = '#05889e', size = 1.5, linetype = 'dashed') +
  geom_errorbar(data = NestFE_Data, aes(Prox, Mean, ymin = (Mean - SD), ymax = (Mean + SD)), 
                size = 0.7, width = 0.15, colour = '#05889e') +
  geom_point(data = NestFE_Data, aes(Prox, Mean), color = 'black', size = 1.3) +
  xlab('Difference in Sorrounding Forest Amount (log(x+1))') +
  ylab('Nestedness Resultant Dissimilarity') +
  scale_y_continuous(breaks = seq(0, 0.8, 0.2), limits = c(0, 0.8)) +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'none')

PairwiseProxPlot <- fs::path('./Figures/',  "Pairwise_Prox_Dissimilarity_Plot.png")
agg_png(PairwiseProxPlot, width = 22, height = 30, units = "cm", res = 300)

pair_prox_plot <- ggarrange(Sor_Prox_Plot, SorFE_Prox_Plot, Turn_Prox_Plot, TurnFE_Prox_Plot, Nest_Prox_Plot, NestFE_Prox_Plot,
                       nrow = 3, ncol = 2, labels = c('A)', '', 'B)', '','C)', ''))
annotate_figure(pair_prox_plot, fig.lab = '                   Overall Community                                                          Forest Specialists',
                fig.lab.pos = 'top.left', fig.lab.face = 'bold', fig.lab.size = 12.5)

invisible(dev.off())
knitr::include_graphics(PairwiseProxPlot)

#####Assess Beta Model coefficients#####

#Do coefficients overlap 0?
sor_summ <- MCMCsummary(sor_quad_out, probs = c(0.05, 0.95), pg0 = T)
turn_summ <- MCMCsummary(turn_quad_out, probs = c(0.05, 0.95), pg0 = T)
nest_summ <- MCMCsummary(nest_out, probs = c(0.05, 0.95), pg0 = T)
sorFE_summ <- MCMCsummary(sorFE_quad_out, probs = c(0.05, 0.95), pg0 = T)
turnFE_summ <- MCMCsummary(turnFE_quad_out, probs = c(0.05, 0.95), pg0 = T)
nestFE_summ <- MCMCsummary(nestFE_out, probs = c(0.05, 0.95), pg0 = T)

#Do coefficients differ between overall and forest endemic community?
#Structure model
beta_comp_model <- nimbleCode({
  
  #Priors
  for(v in 1:10){b[v] ~ dnorm(0, 0.0001)}
  sd.patch ~ dunif(0, 10) 
  tau.patch <- pow(sd.patch, -2)
  
  #Likelihood
  for(i in 1:n){
    
    tau[i] <- pow(sd[i], -2) 
    eps[i] ~ dnorm(0, tau.patch) 
    logit(mu[i]) <- b[1] + b[3] * area[i] + b[4] * prox[i] + b[5] * shape[i] + b[6] * dist[i] + 
      b[2] * comm[i] + b[7] * (area[i] * comm[i]) + b[8] * (prox[i] * comm[i]) + 
      b[9] * (shape[i] * comm[i]) + b[10] * (dist[i] * comm[i]) + eps[i] 
    x[i] ~ dnorm(mu[i], tau[i]) 
    
  }
  
})

beta_quad_comp_model <- nimbleCode({
  
  #Priors
  for(v in 1:12){b[v] ~ dnorm(0, 0.0001)}
  sd.patch ~ dunif(0, 10) 
  tau.patch <- pow(sd.patch, -2)
  
  #Likelihood
  for(i in 1:n){
    
    tau[i] <- pow(sd[i], -2) 
    eps[i] ~ dnorm(0, tau.patch) 
    logit(mu[i]) <- b[1] + b[3] * area[i] + b[4] * area2[i] + b[5] * prox[i] + b[6] * shape[i] + b[7] * dist[i] + 
      b[2] * comm[i] + b[8] * (area[i] * comm[i]) + b[9] * (area2[i] * comm[i]) + b[10] * (prox[i] * comm[i]) + 
      b[11] * (shape[i] * comm[i]) + b[12] * (dist[i] * comm[i]) + eps[i] 
    x[i] ~ dnorm(mu[i], tau[i]) 
    
  }
  
})

beta_comp_consts <- list(area = c(area,area), prox = c(prox,prox), 
                         shape = c(shape,shape), dist = c(dist,dist),
                         n = (((nPatches * (nPatches-1))/2)*2),
                         comm = c(rep(0, ((nPatches * (nPatches-1))/2)), rep(1, ((nPatches * (nPatches-1))/2))))
beta_quad_comp_consts <- list(area = c(area,area), area2 = c(area2, area2), prox = c(prox,prox), 
                         shape = c(shape,shape), dist = c(dist,dist),
                         n = (((nPatches * (nPatches-1))/2)*2),
                         comm = c(rep(0, ((nPatches * (nPatches-1))/2)), rep(1, ((nPatches * (nPatches-1))/2))))
sor_comp_data <- list(x = c(sor_mean, sorFE_mean), sd = c(sor_sd, sorFE_sd))
turn_comp_data <- list(x = c(turn_mean, turnFE_mean), sd = c(turn_sd, turnFE_sd))
nest_comp_data <- list(x = c(nest_mean, nestFE_mean), sd = c(nest_sd, nestFE_sd))

outs <- c('b')
inits <- list(b = rnorm(10))
quad_inits <- list(b = rnorm(12))

#Sorenson and turnover
model <- nimbleModel(beta_quad_comp_model, constants = beta_quad_comp_consts, 
                     data = sor_comp_data, inits = quad_inits)
MCMCconf <- configureMCMC(model, monitors = outs)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)

sor_comp_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                         thin = 10, nchains = 4, samplesAsCodaMCMC = T)

compModel$setData(turn_comp_data)
turn_comp_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                        thin = 10, nchains = 4, samplesAsCodaMCMC = T)

#Nestedness
model <- nimbleModel(beta_comp_model, constants = beta_comp_consts, 
                     data = nest_comp_data, inits = inits)
MCMCconf <- configureMCMC(model, monitors = outs)
MCMC <- buildMCMC(MCMCconf)
compModel <- compileNimble(model)
compMCMC <- compileNimble(MCMC, project = compModel)
nest_comp_out <- runMCMC(compMCMC, niter = 50000, nburnin = 10000, 
                         thin = 10, nchains = 4, samplesAsCodaMCMC = T)

Sor_comp_summ <- MCMCsummary(sor_comp_out, probs = c(0.05, 0.95), pg0 = T)
Turn_comp_summ <- MCMCsummary(turn_comp_out, probs = c(0.05, 0.95), pg0 = T)
Nest_comp_summ <- MCMCsummary(nest_comp_out, probs = c(0.05, 0.95), pg0 = T)


#####Landscape level dissimilarity (i.e. multi site)#####

sps_FHP <- read.csv('AF_DBs_Fragment-Related_Habitat_Preferences.csv')
sps_order <- gsub('_', ' ', colnames(dbs_wide)[26:108])
sps_FHP <- filter(sps_FHP, Binomial %in% sps_order)
sps_FHP$Binomial == sps_order

zmat <- MCMCpstr(db_results, params = 'z', type = 'chains')[[1]]
zmat_fs <- zmat

for(i in 1:dim(zmat)[3]){
  for(n in 1:nrow(dbs_wide)){
  
  zmat_fs[n,,i] <- zmat[n,,i] * sps_FHP$FE 
  
  }
}

sor <- turn <- nest <- sor_fs <- turn_fs <- nest_fs <- vector('numeric', dim(zmat)[3])

for(i in 1:dim(zmat)[3]){
  
  overall <- beta.multi(zmat[,,i])
  fs <- beta.multi(zmat_fs[,,i])
  
  sor[i] <- overall[3]
  turn[i] <- overall[1]
  nest[i] <- overall[2]
  sor_fs[i] <- fs[3]
  turn_fs[i] <- fs[1]
  nest_fs[i] <- fs[2]
  
}

mean(unlist(sor))
quantile(unlist(sor), c(0.05,0.95))
mean(unlist(sor_fs))
quantile(unlist(sor_fs), c(0.05,0.95))

mean(unlist(turn))
quantile(unlist(turn), c(0.05,0.95))
mean(unlist(turn_fs))
quantile(unlist(turn_fs), c(0.05,0.95))


mean(unlist(nest))
quantile(unlist(nest), c(0.05,0.95))
mean(unlist(nest_fs))
quantile(unlist(nest_fs), c(0.05,0.95))

#####Imperfect Detection Plots#####

dbs_inc <- ifelse(dbs_wide[,26:108] > 0, 1, 0)
obsR <- apply(dbs_inc, 1, sum)
estR <- apply(patch_chains, 1, median)
estRuci <- rep(apply(patch_chains, 1, quantile, 0.95), 2)
estRlci <- rep(apply(patch_chains, 1, quantile, 0.05), 2)
mean(estR-obsR)
sd(estR-obsR)

id_richness <- data.frame(patch = dbs_wide$Fragment, 
                          area = dbs_wide$Area_Ha, 
                          obsR, 
                          estR)

estRlci <- estRlci[order(rep(id_richness$area, 2))][c(1:42, 43, 45, 44, 46)]
estRuci <- estRuci[order(rep(id_richness$area, 2))][c(1:42, 43, 45, 44, 46)]
id_richness$estR <- id_richness$estR - id_richness$obsR
id_richness <- id_richness[order(id_richness$area),]
id_richness$patch <- 1:23
id_richness <- pivot_longer(id_richness,
                            cols = c(obsR, estR),
                            names_to = 'Estimator',
                            values_to = 'Richness')

estRplot <- fs::path('./Figures/',  "Richness_Estimates_Plot.png")
agg_png(estRplot, width = 30, height = 20, units = "cm", res = 300) 

ggplot(id_richness, aes(fill = Estimator, y = Richness, x = patch)) +
  geom_bar(position = 'stack', stat = 'identity') +
  geom_errorbar(aes(ymin = estRlci, ymax = estRuci), width = .2) +
  scale_fill_manual(values = c('#58dff5', '#F6C65A'), labels = c( 'Estimated', 'Observed')) +
  scale_y_continuous(breaks = seq(0, 70, 5), labels = seq(0, 70, 5), limits = c(0, 70), expand = c(0,0)) + 
  scale_x_continuous(breaks = 1:23, labels = c(paste0('P', 21:1), 'CF1', 'CF2'), expand = c(0.01,0)) +
  labs(y = 'Estimated Species Richness',
       x = 'Forest Site') +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 0),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(face = 'bold', size = 9))


invisible(dev.off())
knitr::include_graphics(estRplot)

#No. patches occupied by species
obsOcc <- apply(dbs_inc, 2, sum)
estOcc <- MCMCpstr(db_results, params = 'nOcc', type = 'chains')[[1]]

estOccMed <- apply(estOcc, 1, median)
estOccLCI <- apply(estOcc, 1, quantile, 0.05)
estOccUCI <- apply(estOcc, 1, quantile, 0.95)

mean(estOccMed-obsOcc)
sd(estOccMed-obsOcc)

estOccLCI <- rep(estOccLCI, each = 2)
estOccUCI <- rep(estOccUCI, each = 2)

id_occ <- data.frame(Species = SppBinoms, 
                     obsOcc, 
                     estOccMed)

id_occ$estOccMed <- id_occ$estOccMed - id_occ$obsOcc
id_occ <- pivot_longer(id_occ,
                       cols = c(obsOcc, estOccMed),
                       names_to = 'Estimator',
                       values_to = 'nOcc')

estOccplot <- fs::path('./Figures/',  "Occupancy_Estimates_Plot.png")
agg_png(estOccplot, width = 20, height = 30, units = "cm", res = 300) 

ggplot(id_occ, aes(fill = Estimator, y = Species, x = nOcc, width = .75)) +
  geom_bar(position = 'stack', stat = 'identity') +
  geom_errorbar(aes(xmin = estOccLCI, xmax = estOccUCI), width = .5) +
  scale_fill_manual(values = c('#58dff5', '#F6C65A'), labels = c( 'Estimated', 'Observed')) +
  scale_x_continuous(breaks = seq(0, 23, 2), labels = seq(0, 23, 2), limits = c(0, 23), expand = c(0,0)) +
  scale_y_discrete(expand = c(0.02, 0)) +
  labs(x = 'Estimated Number of Occupied Sites',
       y = 'Species') +
  theme_bw() +
  theme(plot.margin = unit(c(.75,.75,.75,.75), 'cm'),
        axis.line = element_line(colour = 'black'),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_text(face = 'bold', vjust = -1.5),
        axis.title.y = element_text(face = 'bold', vjust = 2.75),
        axis.text.x = element_text(face = 'bold', colour = 'black', vjust = -0.9),
        axis.text.y = element_text(face = 'bold', colour = 'black', hjust = 1, size = 8),
        axis.line.x = element_line(size = 0.75),
        axis.line.y = element_line(size = 0.75),
        legend.position = 'top',
        legend.title = element_blank(),
        legend.text = element_text(face = 'bold', size = 9))

invisible(dev.off())
knitr::include_graphics(estOccplot)