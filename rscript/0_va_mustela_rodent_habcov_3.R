###################################################################################################################
##  A dynamic occupancy model for interacting species with two spatial scales model analyzing data from         ##  
##  a long-term monitoring program of small mammals on the arctic tundra                                         ##
##                 by Eivind Flittie Kleiven and Frederic Barraquand                                             ##
##    Last updated 27.1.21                                                                                       ##
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
    
    # site parameters
    #gamA  ~ dunif(0,1)
    #gamB  ~ dunif(0,1)
    #gamAB ~ dunif(0,1)
    #gamBA ~ dunif(0,1)
    #epsA  ~ dunif(0,1)
    #epsB  ~ dunif(0,1)
    #epsAB ~ dunif(0,1)
    #epsBA ~ dunif(0,1)
    
    # block parameters
    #GamA  ~ dunif(0,1)
    #GamB  ~ dunif(0,1)
    #GamAB ~ dunif(0,1)
    #GamBA ~ dunif(0,1)
    #EpsA  ~ dunif(0,1)
    #EpsB  ~ dunif(0,1)
    #EpsAB ~ dunif(0,1)
    #EpsBA ~ dunif(0,1)
    
    # interscept det prob
    for(i in 1:2){
      beta1[i]  ~ dunif(0,1) 
      beta2[i]  ~ dunif(0,1)
      beta3[i]  ~ dunif(0,1) 
      beta4[i]  ~ dunif(0,1)
      beta5[i]  ~ dunif(0,1) 
      beta6[i]  ~ dunif(0,1)
      beta7[i]  ~ dunif(0,1) 
      beta8[i]  ~ dunif(0,1)

      beta0_gamA[i]  <- logit(beta1[i])
      beta0_gamAB[i] <- logit(beta2[i])
      beta0_gamB[i]  <- logit(beta3[i])
      beta0_gamBA[i] <- logit(beta4[i])
    
      beta0_epsA[i]  <- logit(beta5[i])
      beta0_epsAB[i] <- logit(beta6[i])
      beta0_epsB[i]  <- logit(beta7[i])
      beta0_epsBA[i] <- logit(beta8[i])
    } # end loop
    
      beta9  ~ dunif(0,1) 
      beta10 ~ dunif(0,1)
      beta11 ~ dunif(0,1) 
      beta12 ~ dunif(0,1)
      beta13 ~ dunif(0,1) 
      beta14 ~ dunif(0,1)
      beta15 ~ dunif(0,1) 
      beta16 ~ dunif(0,1)
      
      beta0_GamA  <- logit(beta9)
      beta0_GamAB <- logit(beta10)
      beta0_GamB  <- logit(beta11)
      beta0_GamBA <- logit(beta12)
    
      beta0_EpsA  <- logit(beta13)
      beta0_EpsAB <- logit(beta14)
      beta0_EpsB  <- logit(beta15)
      beta0_EpsBA <- logit(beta16)
      
    # interscept det prob
    alpha1 ~ dunif(0,1) 
    alpha2 ~ dunif(0,1)
    
    alphaA0 <- logit(alpha1)
    alphaB0 <- logit(alpha2)
    
  # initial state parameters
    for(b in 1:nblock){
    for(i in 1:3){
    psi[b,i] ~ dunif(0,0.5) # site
    }
  
      # for block, which is just a function of the states of the sites within each block
   
     x[b,1] <- ifelse(sum(z[,b,1]==1) == nsite, 1,
                ifelse(sum(z[,b,1]==2) + sum(z[,b,1]==1) == nsite, 2,
                  ifelse(sum(z[,b,1]==3) + sum(z[,b,1]==1) == nsite, 3, 4) ) )
  
  
    ###############################################
    # btpm = block transition probability matrix. #
    # All columns sum to 1.                       #
    ###############################################
    
    # U to ...
    btpm[b, 1, 1] <- (1-GamA[b]) * (1-GamB[b])    #--|U
    btpm[b, 2, 1] <- GamA[b] * (1-GamB[b])        #--|A
    btpm[b, 3, 1] <- (1-GamA[b]) * GamB[b]        #--|B
    btpm[b, 4, 1] <- GamA[b] * GamB[b]            #--|AB
    
    # A to ...
    btpm[b, 1, 2] <- EpsA[b] * (1-GamBA[b])       #--|U
    btpm[b, 2, 2] <- (1-EpsA[b]) * (1-GamBA[b])   #--|A
    btpm[b, 3, 2] <- EpsA[b] * GamBA[b]           #--|B
    btpm[b, 4, 2] <- (1-EpsA[b]) * GamBA[b]       #--|AB
    
    # B to ...
    btpm[b, 1, 3] <- (1-GamAB[b]) * EpsB[b]       #--|U
    btpm[b, 2, 3] <- GamAB[b] * EpsB[b]           #--|A
    btpm[b, 3, 3] <- (1-GamAB[b]) * (1-EpsB[b])   #--|B
    btpm[b, 4, 3] <- GamAB[b] * (1-EpsB[b])       #--|AB
    
    # AB to ..
    btpm[b, 1, 4] <- EpsAB[b] * EpsBA[b]          #--|U
    btpm[b, 2, 4] <- (1-EpsAB[b]) * EpsBA[b]      #--|A
    btpm[b, 3, 4] <- EpsAB[b] * (1-EpsBA[b])      #--|B
    btpm[b, 4, 4] <- (1-EpsAB[b]) * (1-EpsBA[b])  #--|AB



    for(t in 1:(nseason-1)){
    
   # latent block state for the rest of the seasons
    x[b, t+1] ~ dcat(btpm[b, (1:nout), x[b, t]])
    } # end time loop



    ####################################################
    ## stpm = site transition probability matrix.     ##
    ## These are dependent on the block level state   ##
    ## All columns sum to 1.                          ##
    ####################################################
 
 ## Logit links for col and ext probabilities
 
    logit(GamA[b])  <- beta0_GamA 
    logit(GamAB[b]) <- beta0_GamAB 
    logit(GamB[b])  <- beta0_GamB 
    logit(GamBA[b]) <- beta0_GamBA 
 
    logit(EpsA[b])  <- beta0_EpsA 
    logit(EpsAB[b]) <- beta0_EpsAB 
    logit(EpsB[b])  <- beta0_EpsB 
    logit(EpsBA[b]) <- beta0_EpsBA
    
 
   for(j in 1:nsite){
   
    logit(gamA[b, j])  <- beta0_gamA[hab[b, j]] 
    logit(gamAB[b, j]) <- beta0_gamAB[hab[b, j]] 
    logit(gamB[b, j])  <- beta0_gamB[hab[b, j]] 
    logit(gamBA[b, j]) <- beta0_gamBA[hab[b, j]] 
 
    logit(epsA[b, j])  <- beta0_epsA[hab[b, j]] 
    logit(epsAB[b, j]) <- beta0_epsAB[hab[b, j]] 
    logit(epsB[b, j])  <- beta0_epsB[hab[b, j]]
    logit(epsBA[b, j]) <- beta0_epsBA[hab[b, j]] 
  
   
    # blocks state (x) = U
    
    # U to ...
    stpm[b, j, 1, 1, 1] <- 1          #--|U
    stpm[b, j, 2, 1, 1] <- 0          #--|A
    stpm[b, j, 3, 1, 1] <- 0          #--|B
    stpm[b, j, 4, 1, 1] <- 0          #--|AB
    
    # A to ...
    stpm[b, j, 1, 2, 1] <- 1          #--|U
    stpm[b, j, 2, 2, 1] <- 0          #--|A
    stpm[b, j, 3, 2, 1] <- 0          #--|B
    stpm[b, j, 4, 2, 1] <- 0          #--|AB
    
    # B to ...
    stpm[b, j, 1, 3, 1] <- 1           #--|U
    stpm[b, j, 2, 3, 1] <- 0           #--|A
    stpm[b, j, 3, 3, 1] <- 0           #--|B
    stpm[b, j, 4, 3, 1] <- 0           #--|AB
    
    # AB to ..
    stpm[b, j, 1, 4, 1] <- 1           #--|U
    stpm[b, j, 2, 4, 1] <- 0           #--|A
    stpm[b, j, 3, 4, 1] <- 0           #--|B
    stpm[b, j, 4, 4, 1] <- 0           #--|AB

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 # blocks state (x) = A
  
# U to ...
    stpm[b, j, 1, 1, 2] <- (1-gamA[b, j])    #--|U
    stpm[b, j, 2, 1, 2] <- gamA[b, j]        #--|A
    stpm[b, j, 3, 1, 2] <- 0           #--|B
    stpm[b, j, 4, 1, 2] <- 0           #--|AB
    
    # A to ...
    stpm[b, j, 1, 2, 2] <- epsA[b, j]        #--|U
    stpm[b, j, 2, 2, 2] <- (1-epsA[b, j])    #--|A
    stpm[b, j, 3, 2, 2] <- 0           #--|B
    stpm[b, j, 4, 2, 2] <- 0           #--|AB
    
    # B to ...
    stpm[b, j, 1, 3, 2] <- (1-gamAB[b, j])   #--|U
    stpm[b, j, 2, 3, 2] <- gamAB[b, j]       #--|A
    stpm[b, j, 3, 3, 2] <- 0           #--|B
    stpm[b, j, 4, 3, 2] <- 0           #--|AB
    
    # AB to ..
    stpm[b, j, 1, 4, 2] <- epsAB[b, j]       #--|U
    stpm[b, j, 2, 4, 2] <- (1-epsAB[b, j])   #--|A
    stpm[b, j, 3, 4, 2] <- 0           #--|B
    stpm[b, j, 4, 4, 2] <- 0           #--|AB

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 # blocks state (x) = B
 
# U to ...
    stpm[b, j, 1, 1, 3] <-  (1-gamB[b, j])    #--|U
    stpm[b, j, 2, 1, 3] <- 0            #--|A
    stpm[b, j, 3, 1, 3] <-  gamB[b, j]        #--|B
    stpm[b, j, 4, 1, 3] <- 0            #--|AB
    
    # A to ...
    stpm[b, j, 1, 2, 3] <- (1-gamBA[b, j])    #--|U
    stpm[b, j, 2, 2, 3] <- 0            #--|A
    stpm[b, j, 3, 2, 3] <- gamBA[b, j]       #--|B
    stpm[b, j, 4, 2, 3] <- 0            #--|AB
    
    # B to ...
    stpm[b, j, 1, 3, 3] <-  epsB[b, j]        #--|U
    stpm[b, j, 2, 3, 3] <- 0            #--|A
    stpm[b, j, 3, 3, 3] <- (1-epsB[b, j])     #--|B
    stpm[b, j, 4, 3, 3] <- 0            #--|AB
    
    # AB to ..
    stpm[b, j, 1, 4, 3] <-  epsBA[b, j]       #--|U
    stpm[b, j, 2, 4, 3] <- 0            #--|A
    stpm[b, j, 3, 4, 3] <-  (1-epsBA[b, j])   #--|B
    stpm[b, j, 4, 4, 3] <- 0            #--|AB

 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 # blocks state (x) = AB
 
# U to ...
    stpm[b, j, 1, 1, 4] <- (1-gamA[b, j]) * (1-gamB[b, j])    #--|U
    stpm[b, j, 2, 1, 4] <- gamA[b, j]  * (1-gamB[b, j])       #--|A
    stpm[b, j, 3, 1, 4] <- (1-gamA[b, j]) * gamB[b, j]        #--|B
    stpm[b, j, 4, 1, 4] <- gamA[b, j] * gamB[b, j]            #--|AB
    
    # A to ...
    stpm[b, j, 1, 2, 4] <- epsA[b, j] * (1-gamBA[b, j])       #--|U
    stpm[b, j, 2, 2, 4] <- (1-epsA[b, j]) * (1-gamBA[b, j])   #--|A
    stpm[b, j, 3, 2, 4] <- epsA[b, j] * gamBA[b, j]           #--|B
    stpm[b, j, 4, 2, 4] <- (1-epsA[b, j]) * gamBA[b, j]       #--|AB
    
    # B to ...
    stpm[b, j, 1, 3, 4] <- (1-gamAB[b, j] ) * epsB[b, j]      #--|U
    stpm[b, j, 2, 3, 4] <- gamAB[b, j]  * epsB[b, j]          #--|A
    stpm[b, j, 3, 3, 4] <- (1-gamAB[b, j] ) * (1-epsB[b, j])  #--|B
    stpm[b, j, 4, 3, 4] <- gamAB[b, j]  * (1-epsB[b, j])      #--|AB
    
    # AB to ..
    stpm[b, j, 1, 4, 4] <- epsAB[b, j] * epsBA[b, j]          #--|U
    stpm[b, j, 2, 4, 4] <- (1-epsAB[b, j]) * epsBA[b, j]      #--|A
    stpm[b, j, 3, 4, 4] <- epsAB[b, j] * (1-epsBA[b, j])      #--|B
    stpm[b, j, 4, 4, 4] <- (1-epsAB[b, j]) * (1-epsBA[b, j])  #--|AB
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 # First season probabilities for each state
    #site 
    
      fsm[j, b, 1] <- 1-psi[b,1]-psi[b,2]-psi[b,3]  #-----------|U
      fsm[j, b, 2] <- psi[b,1]                      #-----------|A
      fsm[j, b, 3] <- psi[b,2]                      #-----------|B
      fsm[j, b, 4] <- psi[b,3]                      #-----------|AB

    
    # first season latent state
    # for sites   
      z[j, b, 1] ~ dcat( fsm[j, b, ( 1:nout )])    # site 
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    for(t in 1:(nseason-1)){
    
      # latent site state for the rest of the seasons
        z[j, b, t+1] ~ dcat( stpm[b, j, ( 1:nout ) , z[ j, b, t], x[b,t+1]] + 0.01)  # +0.01 to avoide giving the dcat a prob of 0 
       
          for(day in 1:nsurvey) {      
            y[j, b, t, day] ~ dcat( dpm[( 1:nout ) , z[j, b, t]] + 0.01)  # +0.01 to avoide giving the dcat a prob of 0 
          } #end survey loop
        } # end time loop
      } # end site loop
        
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

load("occm_var.rda")    
#load("case_study_data.RData")
yb <-occm_va # change name of imported object to fit with the rest of the code

dim(yb) # check that dimensions are ok

## import hab cov
load("hab.rda")

# give data
data <-list(nseason = dim(yb)[3], nblock = dim(yb)[2], nsite = dim(yb)[1], nsurvey = dim(yb)[4], 
            nout=4, y = yb, hab=hab)

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

# give initial values
inits=function(){list( 
  z = sp_inits, alpha1=runif(1), alpha2=runif(1),
  beta1=runif(2,0.1,0.9), beta2=runif(2,0.1,0.9), beta3=runif(2,0.1,0.9), beta4=runif(2,0.1,0.9),
  beta5=runif(2,0.1,0.9), beta6=runif(2,0.1,0.9), beta7=runif(2,0.1,0.9), beta8=runif(2,0.1,0.9), 
  beta9=runif(1,0.1,0.9), beta10=runif(1,0.1,0.9), beta11=runif(1,0.1,0.9), beta12=runif(1,0.1,0.9),
  beta13=runif(1,0.1,0.9), beta14=runif(1,0.1,0.9), beta15=runif(1,0.1,0.9), beta16=runif(1,0.1,0.9)
)}

# Parameters monitored
params <- c("gamA","gamB","gamAB","gamBA","epsA","epsB","epsAB","epsBA","psi",
            "GamA","GamB","GamAB","GamBA","EpsA","EpsB","EpsAB","EpsBA", "pA","pB","z","x",
            "alphaA0","alphaB0", 
            "beta0_gamA", "beta0_gamAB", "beta0_gamB", "beta0_gamBA",
            "beta0_epsA", "beta0_epsAB", "beta0_epsB", "beta0_epsBA",
            "beta0_GamA", "beta0_GamAB", "beta0_GamB", "beta0_GamBA",
            "beta0_EpsA", "beta0_EpsAB", "beta0_EpsB", "beta0_EpsBA" )

# MCMC settings
ni <- 5000   ;   nt <- 10   ;   nb <- 0 ;   nc <- 4    ;   na <- 2500

# run model in jags
setwd("../")

va_mustela_rodent_hab_ni5k <- jags(data, inits=inits, params, "mod_seas_det_4stpm.txt", n.chains = nc,
                                   n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel = T)

# Save model
setwd("./model_output")
save(va_mustela_rodent_hab_ni5k, file="va_mustela_rodent_hab_ni5k.rda")

#~ End of script