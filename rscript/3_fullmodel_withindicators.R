###################################################################################################################
##  A dynamic occupancy model for interacting species with two spatial scales analyzing data from          ##  
##  a long-term monitoring program of small mammals on the arctic tundra                                         ##
##                 by Eivind Flittie Kleiven and Frederic Barraquand                                             ##
##                                                                                                               ##
###################################################################################################################

# Call jags(and other packages)
rm(list=ls())

# set working directory for where the model will be saved
setwd("C:/Eivind/GitProjects/MustelidsAndRodents-")
#setwd("~/UiT/GitProjects/MustelidsAndRodents-")

# Specify model in JAGS language
sink("fullmod_withindicators.txt")
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

for (j in 1:14){

# model label
ind[j]~dbern(.2)

# Create indicators for the possible variables in the model
	TempIndicator[j] <- ind[j]*pow(2, j-1) 
}

for(j in 1:28){
  # prior for regression coef (but intercept)
  beta[j]~dnorm(0,0.0001)
  }
  
    beta0_gamA  ~ dnorm(0,1)
    beta0_gamAB ~ dnorm(0,1)
    beta0_gamB  ~ dnorm(0,1)
    beta0_gamBA ~ dnorm(0,1)
    
    beta0_epsA  ~ dnorm(0,1)
    beta0_epsAB ~ dnorm(0,1)
    beta0_epsB  ~ dnorm(0,1)
    beta0_epsBA ~ dnorm(0,1)
      
    beta0_GamA  ~ dnorm(0,1)
    beta0_GamAB ~ dnorm(0,1)
    beta0_GamB  ~ dnorm(0,1)
    beta0_GamBA ~ dnorm(0,1)
    
    beta0_EpsA  ~ dnorm(0,1)
    beta0_EpsAB ~ dnorm(0,1)
    beta0_EpsB  ~ dnorm(0,1)
    beta0_EpsBA ~ dnorm(0,1)
      
      
  # interscept det prob
    alphaA0 ~ dnorm(0,1)
    alphaB0 ~ dnorm(0,1)
    

    # prior for the seasonal covariat on block level  
    for(b in 1:nblock){
      for(t in 1:(nseason-1)){
        seas_block[b,t] ~ dbern(0.5)
        
        for(j in 1:nsite){
          seas[j,b,t] ~ dbern(0.5)
          }
        }
      }  
      

    
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
 
    logit(GamA[b,t])  <- beta0_GamA  + ind[1] * beta[1] * seas_block[b,t] 
    logit(GamAB[b,t]) <- beta0_GamAB + ind[1] * beta[2] * seas_block[b,t] 
    logit(GamB[b,t])  <- beta0_GamB  + ind[2] * beta[3] * seas_block[b,t] 
    logit(GamBA[b,t]) <- beta0_GamBA + ind[2] * beta[4] * seas_block[b,t] 
 
    logit(EpsA[b,t])  <- beta0_EpsA  + ind[3] * beta[5] * seas_block[b,t] 
    logit(EpsAB[b,t]) <- beta0_EpsAB + ind[3] * beta[6] * seas_block[b,t] 
    logit(EpsB[b,t])  <- beta0_EpsB  + ind[4] * beta[7] * seas_block[b,t] 
    logit(EpsBA[b,t]) <- beta0_EpsBA + ind[4] * beta[8] * seas_block[b,t]
    

   for(j in 1:nsite){
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      # logit links

      logit(gamA[b, j, t])  <- beta0_gamA  + ind[5] * beta[9]  * hab[b, j] + ind[9]  * beta[10] * seas[j, b, t] 
      logit(gamAB[b, j, t]) <- beta0_gamAB + ind[5] * beta[11] * hab[b, j] + ind[9]  * beta[12] * seas[j, b, t] 
      logit(gamB[b, j, t])  <- beta0_gamB  + ind[6] * beta[13] * hab[b, j] + ind[10] * beta[14] * seas[j, b, t] 
      logit(gamBA[b, j, t]) <- beta0_gamBA + ind[6] * beta[15] * hab[b, j] + ind[10] * beta[16] * seas[j, b, t]
 
      logit(epsA[b, j, t])  <- beta0_epsA  + ind[7] * beta[17] * hab[b, j] + ind[11] * beta[18] * seas[j, b, t] 
      logit(epsAB[b, j, t]) <- beta0_epsAB + ind[7] * beta[19] * hab[b, j] + ind[11] * beta[20] * seas[j, b, t]
      logit(epsB[b, j, t])  <- beta0_epsB  + ind[8] * beta[21] * hab[b, j] + ind[12] * beta[22] * seas[j, b, t]
      logit(epsBA[b, j, t]) <- beta0_epsBA + ind[8] * beta[23] * hab[b, j] + ind[12] * beta[24] * seas[j, b, t]
  

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
            y[j, b, t, day] ~ dcat( dpm[b, j, t, ( 1:nout ) , z[j, b, t]] + 0.01)  # +0.01 to avoide giving the dcat a prob of 0 
          } #end survey loop

    
    #############################################################
    ## detection matrix (OS = observed state, TS = true state) ##
    #############################################################
    
    # TS = U
    dpm[b, j, t, 1, 1] <- 1                                    #--|OS = U
    dpm[b, j, t, 2, 1] <- 0                                    #--|OS = A
    dpm[b, j, t, 3, 1] <- 0                                    #--|OS = B
    dpm[b, j, t, 4, 1] <- 0                                    #--|OS = AB
    
    # TS = A
    dpm[b, j, t, 1, 2] <- 1-pA[b, j, t]                       #--|OS = U
    dpm[b, j, t, 2, 2] <- pA[b, j, t]                         #--|OS = A
    dpm[b, j, t, 3, 2] <- 0                                    #--|OS = B
    dpm[b, j, t, 4, 2] <- 0                                    #--|OS = AB
    
    # TS = B
    dpm[b, j, t, 1, 3] <- 1-pB[b, j, t]                       #--|OS = U
    dpm[b, j, t, 2, 3] <- 0                                    #--|OS = A
    dpm[b, j, t, 3, 3] <- pB[b, j, t]                         #--|OS = B
    dpm[b, j, t, 4, 3] <- 0                                    #--|OS = AB
    
    # TS = AB
    dpm[b, j, t, 1, 4] <- (1-pA[b, j, t]) * (1-pB[b, j, t])  #--|OS = U
    dpm[b, j, t, 2, 4] <- pA[b, j, t] * (1-pB[b, j, t])      #--|OS = A
    dpm[b, j, t, 3, 4] <- (1-pA[b, j, t]) * pB[b, j, t]      #--|OS = B
    dpm[b, j, t, 4, 4] <- pA[b, j, t] * pB[b, j, t]          #--|OS = AB
    
    ## logit links for detection probs
    logit(pA[b, j, t]) <- alphaA0 + ind[13] * beta[25] * seas[j, b, t] + ind[14] * beta[26] * hab[b, j]
    logit(pB[b, j, t]) <- alphaB0 + ind[13] * beta[27] * seas[j, b, t] + ind[14] * beta[28] * hab[b, j]  
    
      } # end site loop
        } # end time loop
          } #close block loop
    
    ## Derived parameters
    
    #Create a model number for each possible model
    mdl<- 1+sum(TempIndicator[]) 

    # calculate the percentage of time each model is selected (2^(nb of covariates) = 2^14 = 16384)
    for (j in 1 : 16384 ){ 
    pmdl[j] <- equals(mdl, j) 
    }
    
    
    
    }# end
    ",fill = TRUE)
sink()

