###################################################################################################################
##  A dynamic occupancy model for interacting species with two spatial scales model analyzing data from         ##  
##  a long-term monitoring program of small mammals on the arctic tundra                                         ##
##                 by Eivind Flittie Kleiven and Frederic Barraquand                                             ##
##                                                                                          ##
###################################################################################################################

rm(list=ls())

# Call jags(and other packages)
library(jagsUI)

## import data
setwd("./data") # set wd to where the data is stored

load("occm_var_lem_rmNA.rda")    
#load("case_study_data.RData")
yb <-occm_va_lem[,,161:171,] # change name of imported object to fit with the rest of the code

dim(yb) # check that dimensions are ok

## import hab cov
load("hab.rda")

## import seasonal cov
load("season_cov.rda")
load("season_cov_block.rda")

seas <- season_cov[,,161:171]
seas_block <- season_cov_block[,161:171]

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
  beta0_gamA=runif(1,0.1,0.9), beta0_gamB=runif(1,0.1,0.9), 
  beta0_epsA=runif(1,0.1,0.9), beta0_epsB=runif(1,0.1,0.9),  
  beta0_GamA=runif(1,0.1,0.9), beta0_GamB=runif(1,0.1,0.9), 
  beta0_EpsA=runif(1,0.1,0.9), beta0_EpsB=runif(1,0.1,0.9) 
)}

# Parameters monitored
params <- c("gamA","gamB","gamAB","gamBA","epsA","epsB","epsAB","epsBA","psi",
            "GamA","GamB","GamAB","GamBA","EpsA","EpsB","EpsAB","EpsBA", "pA","pB","z","x",
            "alphaA0","alphaB0",  
            "beta0_gamA",  "beta0_gamB",
            "beta0_epsA",  "beta0_epsB",
            "beta0_GamA",  "beta0_GamB",
            "beta0_EpsA",  "beta0_EpsB" )

# MCMC settings
ni <- 10   ;   nt <- 1   ;   nb <- 0 ;   nc <- 4    ;   na <- 0

# run model in jags
setwd("../")

va_mustela_lemming_mod1 <- jags(data, inits=inits, params, "mod1.txt", n.chains = nc,
                                   n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel = T)

# Save model
setwd("./model_output")
save(va_mustela_lemming_habseas, file="va_mustela_lemming_mod1.rda")

#~ End of script