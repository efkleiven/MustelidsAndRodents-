###################################################################################################################
##  A dynamic occupancy model for interacting species with two spatial scales model analyzing data from         ##  
##  a long-term monitoring program of small mammals on the arctic tundra                                         ##
##                 by Eivind Flittie Kleiven and Frederic Barraquand                                             ##
##                                                                                          ##
###################################################################################################################

rm(list=ls())

# Call jags(and other packages)
library(jagsUI)
require(rjags)

## import data
setwd("./data") # set wd to where the data is stored

load("occm_var_lem_rmNA.rda")    
#load("case_study_data.RData")
yb <-occm_va_lem # change name of imported object to fit with the rest of the code

dim(yb) # check that dimensions are ok

## import hab cov
load("hab.rda")

## import seasonal cov
load("season_cov.rda")
load("season_cov_block.rda")

seas <- season_cov
seas_block <- season_cov_block

#
seas[seas==1]<-0
seas[seas==2]<-1

seas_block[seas_block==1]<-0
seas_block[seas_block==2]<-1

# give data
data <-list(nseason = dim(yb)[3], nblock = dim(yb)[2], nsite = dim(yb)[1], nsurvey = dim(yb)[4], 
            nout=4, y = yb, hab=hab, seas=seas, seas_block=seas_block)

# naming some parameters for loops further down in this script
nseason = dim(yb)[3]; nblock = dim(yb)[2]; nsite = dim(yb)[1]; nsurvey = dim(yb)[4]

# Initial values for state
sp_inits <- apply(yb,c(1,2,3),max)

# loop to make cases where both state 2 and 3 is observed within the same primary occasion have initial value 4
dataL <- array(0,dim=c(nsite,nblock,nseason))
for(j in 1:nsite){
  for(b in 1:nblock){
    for(i in 1:nseason){
      if (is.element(0, match(c(2,3),yb[j,b,i,], nomatch = FALSE, incomparables = FALSE))) 
        dataL[j,b,i] <- "FALSE"
      else
        dataL[j,b,i] <- "TRUE"
    }}}

for(j in 1:nsite){
  for(b in 1:nblock){
    for(i in 1:nseason){
      if(dataL[j,b,i]==TRUE){
        sp_inits[j,b,i] <- 4}
    }}}

# replace NA in initial values with the highest state
#sp_inits[is.na(sp_inits)] <- 4 

# give initial values
inits=function(){list( 
  z = sp_inits, alphaA0=runif(1,0.1,0.9), alphaB0=runif(1,0.1,0.9),
  beta0_gamA=runif(1,0.1,0.9), beta0_gamB=runif(1,0.1,0.9), beta0_gamAB=runif(1,0.1,0.9), beta0_gamBA=runif(1,0.1,0.9),
  beta0_epsA=runif(1,0.1,0.9), beta0_epsB=runif(1,0.1,0.9), beta0_epsAB=runif(1,0.1,0.9), beta0_epsBA=runif(1,0.1,0.9), 
  beta0_GamA=runif(1,0.1,0.9), beta0_GamB=runif(1,0.1,0.9), beta0_GamAB=runif(1,0.1,0.9), beta0_GamBA=runif(1,0.1,0.9),
  beta0_EpsA=runif(1,0.1,0.9), beta0_EpsB=runif(1,0.1,0.9), beta0_EpsAB=runif(1,0.1,0.9), beta0_EpsBA=runif(1,0.1,0.9),
  beta=runif(28, 0.1, 0.9), ind = rep(0,14)
)}

# Parameters monitored
params <- c("gamA","gamB","gamAB","gamBA","epsA","epsB","epsAB","epsBA","psi",
            "GamA","GamB","GamAB","GamBA","EpsA","EpsB","EpsAB","EpsBA", "pA","pB","z","x",
            "alphaA0","alphaB0",  
            "beta0_gamA", "beta0_gamAB", "beta0_gamB", "beta0_gamBA",
            "beta0_epsA", "beta0_epsAB", "beta0_epsB", "beta0_epsBA",
            "beta0_GamA", "beta0_GamAB", "beta0_GamB", "beta0_GamBA",
            "beta0_EpsA", "beta0_EpsAB", "beta0_EpsB", "beta0_EpsBA", 
            "beta", "ind", "pmdl" )

# MCMC settings
#ni <- 10   ;   nt <- 1   ;   nb <- 0 ;   nc <- 4    ;   na <- 0

# run model in jags
setwd("../")

#va_mustela_lemming_fullmod_withindicators <- jags(data, inits=inits, params, "fullmod_withindicators.txt", n.chains = nc,
#                                   n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel = T)


model <- jags.model("fullmod_withindicators.txt", 
                    data = data,
                    inits = inits)

update(model, n.iter = 10)

output <- coda.samples(model = model,
                       variable.names = params, 
                       n.iter = 10,
                       thin = 1,
                       n.adapt=5)


# Save model
setwd("./model_output")
save(va_mustela_lemming_habseas, file="va_mustela_lemming_fullmod_withindicators.rda")

# extra stuff

output <- as.mcmc(output)

print(summary(output))

ind <- function(p){
  if(p == 0) {return(t <- 0)}
  else if(p == 1) {return(t <- rbind(0, 1))}
  else if(p == 2) {return(t <- rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1)))}
  else {
    t <- rbind(cbind(ind(p - 1), rep(0, 2^(p - 1))), 
               cbind(ind(p - 1), rep(1, 2^(p - 1))))
    return(t)
  }
}

# cree toutes les combi possibles si p=5 covariables
# la premiere ligne de mat.modele correspond au modele avec intercept seulement
# la derniere ligne de mat.modele correspond au modele avec toutes les covariables
mat.modeles <- ind(14)
mat.modeles
nrow(mat.modeles)

output[1,]

grep("pmdl", names(output[1,]))
pmdl <- output[,grep("pmdl", names(output[1,]))]
pmp <- apply(pmdl,2,mean) # probabilite a posteriori des modeles

# classe les proba a posteriori des modeles de la plus grande a la plus petite
ii <- order(pmp, decreasing = T)

# affiche les modeles
res <- cbind(mat.modeles[ii,], pmp[ii])
res
#~ End of script