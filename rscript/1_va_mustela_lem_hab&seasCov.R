###################################################################################################################
##  A dynamic occupancy model for interacting species with two spatial scales model analyzing data from         ##  
##  a long-term monitoring program of small mammals on the arctic tundra                                         ##
##                 by Eivind Flittie Kleiven and Frederic Barraquand                                             ##
##                                                                                          ##
###################################################################################################################

# Call jags(and other packages)
rm(list=ls())

library(jagsUI)

# set working directory
#setwd("C:/Eivind/GitProjects/MustelidsAndRodents-")
setwd("~/UiT/GitProjects/MustelidsAndRodents-")

#
# Specify model in JAGS language
sink("mod_seas_det_4stpm.txt")
cat("
    model{
    
    ###########################################
    ### Spatial Dynamic Co-Occurrence Model ###
    ###########################################
    # State 1= Unoccupied(U), State 2= rodent(A), State 3 = mustelid(B), State 4 = rodent & mustelid(AB)  
    ########################################### 
   
    ##############
    ##  Priors  ##    
    ##############

    for(i in 1:2){
      beta0_gamA[i]  ~ dnorm(0,1)
      beta0_gamAB[i] ~ dnorm(0,1)
      beta0_gamB[i]  ~ dnorm(0,1)
      beta0_gamBA[i] ~ dnorm(0,1)
    
      beta0_epsA[i]  ~ dnorm(0,1)
      beta0_epsAB[i] ~ dnorm(0,1)
      beta0_epsB[i]  ~ dnorm(0,1)
      beta0_epsBA[i] ~ dnorm(0,1)
      
      beta1_gamA[i]  ~ dnorm(0,1)
      beta1_gamAB[i] ~ dnorm(0,1)
      beta1_gamB[i]  ~ dnorm(0,1)
      beta1_gamBA[i] ~ dnorm(0,1)
    
      beta1_epsA[i]  ~ dnorm(0,1)
      beta1_epsAB[i] ~ dnorm(0,1)
      beta1_epsB[i]  ~ dnorm(0,1)
      beta1_epsBA[i] ~ dnorm(0,1)
      
      beta0_GamA[i]  ~ dnorm(0,1)
      beta0_GamAB[i] ~ dnorm(0,1)
      beta0_GamB[i]  ~ dnorm(0,1)
      beta0_GamBA[i] ~ dnorm(0,1)
    
      beta0_EpsA[i]  ~ dnorm(0,1)
      beta0_EpsAB[i] ~ dnorm(0,1)
      beta0_EpsB[i]  ~ dnorm(0,1)
      beta0_EpsBA[i] ~ dnorm(0,1)
    } # end loop
    
    # prior for the seasonal covariat on block level  
    for(b in 1:nblock){
      for(t in 1:(nseason-1)){
        seas_block[b,t] ~ dbern(0.5)
        
        for(j in 1:nsite){
          seas[j,b,t] ~ dbern(0.5)
          }
      }
    }  
      
    # interscept det prob

    alphaA0 ~ dnorm(0,1)
    alphaB0 ~ dnorm(0,1)
    
  # initial state parameters
    for(b in 1:nblock){
      for(i in 1:3){
        psi[b,i] ~ dunif(0,0.5) # site
        }
  
      # for block, which is just a function of the states of the sites within each block
   
      x[b,1] <- ifelse(sum(z[,b,1]==1) == nsite, 1,
                 ifelse(sum(z[,b,1]==2) + sum(z[,b,1]==1) == nsite, 2,
                  ifelse(sum(z[,b,1]==3) + sum(z[,b,1]==1) == nsite, 3, 4) ) )
                  
                  
      #site 
      for(j in 1:nsite){
        fsm[j, b, 1] <- 1-psi[b,1]-psi[b,2]-psi[b,3]  #-----------|U
        fsm[j, b, 2] <- psi[b,1]                      #-----------|A
        fsm[j, b, 3] <- psi[b,2]                      #-----------|B
        fsm[j, b, 4] <- psi[b,3]                      #-----------|AB

    
    # first season latent state
    # for sites   
      z[j, b, 1] ~ dcat( fsm[j, b, ( 1:nout )]) 
      } #close site loop
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  
    ###############################################
    # btpm = block transition probability matrix. #
    # All columns sum to 1.                       #
    ###############################################
     for(t in 1:(nseason-1)){
     
    # U to ...
    btpm[b,t, 1, 1] <- (1-GamA[b,t]) * (1-GamB[b,t])    #--|U
    btpm[b,t, 2, 1] <- GamA[b,t] * (1-GamB[b,t])        #--|A
    btpm[b,t, 3, 1] <- (1-GamA[b,t]) * GamB[b,t]        #--|B
    btpm[b,t, 4, 1] <- GamA[b,t] * GamB[b,t]            #--|AB
    
    # A to ...
    btpm[b,t, 1, 2] <- EpsA[b,t] * (1-GamBA[b,t])       #--|U
    btpm[b,t, 2, 2] <- (1-EpsA[b,t]) * (1-GamBA[b,t])   #--|A
    btpm[b,t, 3, 2] <- EpsA[b,t] * GamBA[b,t]           #--|B
    btpm[b,t, 4, 2] <- (1-EpsA[b,t]) * GamBA[b,t]       #--|AB
    
    # B to ...
    btpm[b,t, 1, 3] <- (1-GamAB[b,t]) * EpsB[b,t]       #--|U
    btpm[b,t, 2, 3] <- GamAB[b,t] * EpsB[b,t]           #--|A
    btpm[b,t, 3, 3] <- (1-GamAB[b,t]) * (1-EpsB[b,t])   #--|B
    btpm[b,t, 4, 3] <- GamAB[b,t] * (1-EpsB[b,t])       #--|AB
    
    # AB to ..
    btpm[b,t, 1, 4] <- EpsAB[b,t] * EpsBA[b,t]          #--|U
    btpm[b,t, 2, 4] <- (1-EpsAB[b,t]) * EpsBA[b,t]      #--|A
    btpm[b,t, 3, 4] <- EpsAB[b,t] * (1-EpsBA[b,t])      #--|B
    btpm[b,t, 4, 4] <- (1-EpsAB[b,t]) * (1-EpsBA[b,t])  #--|AB

    
   # latent block state for the rest of the seasons
    x[b, t+1] ~ dcat(btpm[b, t, (1:nout), x[b, t]])
    
    
    ####################################################
    ## stpm = site transition probability matrix.     ##
    ## These are dependent on the block level state   ##
    ## All columns sum to 1.                          ##
    ####################################################
 
 ## Logit links for col and ext probabilities
 
    logit(GamA[b,t])  <- beta0_GamA[seas_block[b,t]+1] 
    logit(GamAB[b,t]) <- beta0_GamAB[seas_block[b,t]+1] 
    logit(GamB[b,t])  <- beta0_GamB[seas_block[b,t]+1] 
    logit(GamBA[b,t]) <- beta0_GamBA[seas_block[b,t]+1] 
 
    logit(EpsA[b,t])  <- beta0_EpsA[seas_block[b,t]+1] 
    logit(EpsAB[b,t]) <- beta0_EpsAB[seas_block[b,t]+1] 
    logit(EpsB[b,t])  <- beta0_EpsB[seas_block[b,t]+1] 
    logit(EpsBA[b,t]) <- beta0_EpsBA[seas_block[b,t]+1]
    
 
     
   for(j in 1:nsite){
   
   #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # logit links

    logit(gamA[b, j, t])  <- beta0_gamA[hab[b, j]]  + beta1_gamA[seas[j, b, t]+1] 
    logit(gamAB[b, j, t]) <- beta0_gamAB[hab[b, j]] + beta1_gamAB[seas[j, b, t]+1] 
    logit(gamB[b, j, t])  <- beta0_gamB[hab[b, j]]  + beta1_gamB[seas[j, b, t]+1] 
    logit(gamBA[b, j, t]) <- beta0_gamBA[hab[b, j]] + beta1_gamBA[seas[j, b, t]+1]
 
    logit(epsA[b, j, t])  <- beta0_epsA[hab[b, j]]  + beta1_epsA[seas[j, b, t]+1] 
    logit(epsAB[b, j, t]) <- beta0_epsAB[hab[b, j]] + beta1_epsAB[seas[j, b, t]+1]
    logit(epsB[b, j, t])  <- beta0_epsB[hab[b, j]]  + beta1_epsB[seas[j, b, t]+1]
    logit(epsBA[b, j, t]) <- beta0_epsBA[hab[b, j]] + beta1_epsBA[seas[j, b, t]+1]
  

# site transition matrix   
    # blocks state (x) = U
    
    # U to ...
    stpm[b, j, t, 1, 1, 1] <- 1          #--|U
    stpm[b, j, t, 2, 1, 1] <- 0          #--|A
    stpm[b, j, t, 3, 1, 1] <- 0          #--|B
    stpm[b, j, t, 4, 1, 1] <- 0          #--|AB
    
    # A to ...
    stpm[b, j, t, 1, 2, 1] <- 1          #--|U
    stpm[b, j, t, 2, 2, 1] <- 0          #--|A
    stpm[b, j, t, 3, 2, 1] <- 0          #--|B
    stpm[b, j, t, 4, 2, 1] <- 0          #--|AB
    
    # B to ...
    stpm[b, j, t, 1, 3, 1] <- 1           #--|U
    stpm[b, j, t, 2, 3, 1] <- 0           #--|A
    stpm[b, j, t, 3, 3, 1] <- 0           #--|B
    stpm[b, j, t, 4, 3, 1] <- 0           #--|AB
    
    # AB to ..
    stpm[b, j, t, 1, 4, 1] <- 1           #--|U
    stpm[b, j, t, 2, 4, 1] <- 0           #--|A
    stpm[b, j, t, 3, 4, 1] <- 0           #--|B
    stpm[b, j, t, 4, 4, 1] <- 0           #--|AB

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 # blocks state (x) = A
  
# U to ...
    stpm[b, j, t, 1, 1, 2] <- (1-gamA[b, j, t])    #--|U
    stpm[b, j, t, 2, 1, 2] <- gamA[b, j, t]        #--|A
    stpm[b, j, t, 3, 1, 2] <- 0                    #--|B
    stpm[b, j, t, 4, 1, 2] <- 0                    #--|AB
    
    # A to ...
    stpm[b, j, t, 1, 2, 2] <- epsA[b, j, t]        #--|U
    stpm[b, j, t, 2, 2, 2] <- (1-epsA[b, j, t])    #--|A
    stpm[b, j, t, 3, 2, 2] <- 0                    #--|B
    stpm[b, j, t, 4, 2, 2] <- 0                    #--|AB
    
    # B to ...
    stpm[b, j, t, 1, 3, 2] <- (1-gamAB[b, j, t])   #--|U
    stpm[b, j, t, 2, 3, 2] <- gamAB[b, j, t]       #--|A
    stpm[b, j, t, 3, 3, 2] <- 0                    #--|B
    stpm[b, j, t, 4, 3, 2] <- 0                    #--|AB
    
    # AB to ..
    stpm[b, j, t, 1, 4, 2] <- epsAB[b, j, t]       #--|U
    stpm[b, j, t, 2, 4, 2] <- (1-epsAB[b, j, t])   #--|A
    stpm[b, j, t, 3, 4, 2] <- 0                    #--|B
    stpm[b, j, t, 4, 4, 2] <- 0                    #--|AB

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 # blocks state (x) = B
 
# U to ...
    stpm[b, j, t, 1, 1, 3] <-  (1-gamB[b, j, t])    #--|U
    stpm[b, j, t, 2, 1, 3] <- 0                     #--|A
    stpm[b, j, t, 3, 1, 3] <-  gamB[b, j, t]        #--|B
    stpm[b, j, t, 4, 1, 3] <- 0                     #--|AB
    
    # A to ...
    stpm[b, j, t, 1, 2, 3] <- (1-gamBA[b, j, t])    #--|U
    stpm[b, j, t, 2, 2, 3] <- 0                     #--|A
    stpm[b, j, t, 3, 2, 3] <- gamBA[b, j, t]        #--|B
    stpm[b, j, t, 4, 2, 3] <- 0                     #--|AB
    
    # B to ...
    stpm[b, j, t, 1, 3, 3] <-  epsB[b, j, t]        #--|U
    stpm[b, j, t, 2, 3, 3] <- 0                     #--|A
    stpm[b, j, t, 3, 3, 3] <- (1-epsB[b, j, t])     #--|B
    stpm[b, j, t, 4, 3, 3] <- 0                     #--|AB
    
    # AB to ..
    stpm[b, j, t, 1, 4, 3] <-  epsBA[b, j, t]       #--|U
    stpm[b, j, t, 2, 4, 3] <- 0                     #--|A
    stpm[b, j, t, 3, 4, 3] <-  (1-epsBA[b, j, t])   #--|B
    stpm[b, j, t, 4, 4, 3] <- 0                     #--|AB

 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 # blocks state (x) = AB
 
# U to ...
    stpm[b, j, t, 1, 1, 4] <- (1-gamA[b, j, t]) * (1-gamB[b, j, t])    #--|U
    stpm[b, j, t, 2, 1, 4] <- gamA[b, j, t]  * (1-gamB[b, j, t])       #--|A
    stpm[b, j, t, 3, 1, 4] <- (1-gamA[b, j, t]) * gamB[b, j, t]        #--|B
    stpm[b, j, t, 4, 1, 4] <- gamA[b, j, t] * gamB[b, j, t]            #--|AB
    
    # A to ...
    stpm[b, j, t, 1, 2, 4] <- epsA[b, j, t] * (1-gamBA[b, j, t])       #--|U
    stpm[b, j, t, 2, 2, 4] <- (1-epsA[b, j, t]) * (1-gamBA[b, j, t])   #--|A
    stpm[b, j, t, 3, 2, 4] <- epsA[b, j, t] * gamBA[b, j, t]           #--|B
    stpm[b, j, t, 4, 2, 4] <- (1-epsA[b, j, t]) * gamBA[b, j, t]       #--|AB
    
    # B to ...
    stpm[b, j, t, 1, 3, 4] <- (1-gamAB[b, j, t] ) * epsB[b, j, t]      #--|U
    stpm[b, j, t, 2, 3, 4] <- gamAB[b, j, t]  * epsB[b, j, t]          #--|A
    stpm[b, j, t, 3, 3, 4] <- (1-gamAB[b, j, t] ) * (1-epsB[b, j, t])  #--|B
    stpm[b, j, t, 4, 3, 4] <- gamAB[b, j, t]  * (1-epsB[b, j, t])      #--|AB
    
    # AB to ..
    stpm[b, j, t, 1, 4, 4] <- epsAB[b, j, t] * epsBA[b, j, t]          #--|U
    stpm[b, j, t, 2, 4, 4] <- (1-epsAB[b, j, t]) * epsBA[b, j, t]      #--|A
    stpm[b, j, t, 3, 4, 4] <- epsAB[b, j, t] * (1-epsBA[b, j, t])      #--|B
    stpm[b, j, t, 4, 4, 4] <- (1-epsAB[b, j, t]) * (1-epsBA[b, j, t])  #--|AB

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
      # latent site state for the rest of the seasons
        z[j, b, t+1] ~ dcat( stpm[b, j, t, ( 1:nout ) , z[ j, b, t], x[b,t+1]] + 0.01)  # +0.01 to avoide giving the dcat a prob of 0 
       
          for(day in 1:nsurvey) {      
            y[j, b, t, day] ~ dcat( dpm[( 1:nout ) , z[j, b, t]] + 0.01)  # +0.01 to avoide giving the dcat a prob of 0 
          } #end survey loop
        } # end site loop
      } # end time loop
    } #close block loop
    
    #############################################################
    ## detection matrix (OS = observed state, TS = true state) ##
    #############################################################
    
    # TS = U
    dpm[1, 1] <- 1                      #--|OS = U
    dpm[2, 1] <- 0                      #--|OS = A
    dpm[3, 1] <- 0                      #--|OS = B
    dpm[4, 1] <- 0                      #--|OS = AB
    
    # TS = A
    dpm[1, 2] <- 1-pA                #--|OS = U
    dpm[2, 2] <- pA                  #--|OS = A
    dpm[3, 2] <- 0                      #--|OS = B
    dpm[4, 2] <- 0                      #--|OS = AB
    
    # TS = B
    dpm[ 1, 3] <- 1-pB                #--|OS = U
    dpm[ 2, 3] <- 0                      #--|OS = A
    dpm[ 3, 3] <- pB                  #--|OS = B
    dpm[ 4, 3] <- 0                      #--|OS = AB
    
    # TS = AB
    dpm[ 1, 4] <- (1-pA) * (1-pB)  #--|OS = U
    dpm[ 2, 4] <- pA * (1-pB)      #--|OS = A
    dpm[ 3, 4] <- (1-pA) * pB      #--|OS = B
    dpm[ 4, 4] <- pA * pB          #--|OS = AB
    
    ## logit links for detection probs
    logit(pA) <- alphaA0 
    logit(pB) <- alphaB0 
    
    
    ## Derived parameters
    
    #diff_gamA <- gamA - gamAB
    #diff_gamB <- gamB - gamBA
    #diff_epsA <- epsA - epsAB
    #diff_epsB <- epsB - epsBA
    
    #diff_GamA <- GamA - GamAB
    #diff_GamB <- GamB - GamBA
    #diff_EpsA <- EpsA - EpsAB
    #diff_EpsB <- EpsB - EpsBA
    
    #ratio_gamA <- gamA / gamAB
    #ratio_gamB <- gamB / gamBA
    #ratio_epsA <- epsAB / epsA
    #ratio_epsB <- epsB / epsBA
    
    #ratio_GamA <- GamA / GamAB
    #ratio_GamB <- GamBA / GamB
    #ratio_EpsA <- EpsAB / EpsA
    #ratio_EpsB <- EpsB / EpsBA    
    
    }# end
    ",fill = TRUE)
sink()

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
  z = sp_inits, alphaA0=runif(1), alphaB0=runif(1),
  beta0_gamA=runif(2,0.1,0.9), beta0_gamB=runif(2,0.1,0.9), beta0_gamAB=runif(2,0.1,0.9), beta0_gamBA=runif(2,0.1,0.9),
  beta0_epsA=runif(2,0.1,0.9), beta0_epsB=runif(2,0.1,0.9), beta0_epsAB=runif(2,0.1,0.9), beta0_epsBA=runif(2,0.1,0.9), 
  beta1_gamA=runif(2,0.1,0.9), beta1_gamB=runif(2,0.1,0.9), beta1_gamAB=runif(2,0.1,0.9), beta1_gamBA=runif(2,0.1,0.9),
  beta1_epsA=runif(2,0.1,0.9), beta1_epsB=runif(2,0.1,0.9), beta1_epsAB=runif(2,0.1,0.9), beta1_epsBA=runif(2,0.1,0.9),
  beta0_GamA=runif(2,0.1,0.9), beta0_GamB=runif(2,0.1,0.9), beta0_GamAB=runif(2,0.1,0.9), beta0_GamBA=runif(2,0.1,0.9),
  beta0_EpsA=runif(2,0.1,0.9), beta0_EpsB=runif(2,0.1,0.9), beta0_EpsAB=runif(2,0.1,0.9), beta0_EpsBA=runif(2,0.1,0.9) 
)}

# Parameters monitored
params <- c("gamA","gamB","gamAB","gamBA","epsA","epsB","epsAB","epsBA","psi",
            "GamA","GamB","GamAB","GamBA","EpsA","EpsB","EpsAB","EpsBA", "pA","pB","z","x",
            "alphaA0","alphaB0", 
            "beta0_gamA", "beta0_gamAB", "beta0_gamB", "beta0_gamBA",
            "beta0_epsA", "beta0_epsAB", "beta0_epsB", "beta0_epsBA",
            "beta1_gamA", "beta1_gamAB", "beta1_gamB", "beta1_gamBA",
            "beta1_epsA", "beta1_epsAB", "beta1_epsB", "beta1_epsBA",
            "beta0_GamA", "beta0_GamAB", "beta0_GamB", "beta0_GamBA",
            "beta0_EpsA", "beta0_EpsAB", "beta0_EpsB", "beta0_EpsBA" )

# MCMC settings
ni <- 5000   ;   nt <- 5   ;   nb <- 0 ;   nc <- 4    ;   na <- 1000

# run model in jags
setwd("../")

va_mustela_lemming_habseas <- jags(data, inits=inits, params, "mod_seas_det_4stpm.txt", n.chains = nc,
                                   n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel = T)

# Save model
setwd("./model_output")
save(va_mustela_lemming_habseas, file="va_mustela_lemming_habseas_2.rda")

#~ End of script