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
    beta1  ~ dunif(0,1) 
    beta2  ~ dunif(0,1)
    beta3  ~ dunif(0,1) 
    beta4  ~ dunif(0,1)
    beta5  ~ dunif(0,1) 
    beta6  ~ dunif(0,1)
    beta7  ~ dunif(0,1) 
    beta8  ~ dunif(0,1)
    beta9  ~ dunif(0,1) 
    beta10 ~ dunif(0,1)
    beta11 ~ dunif(0,1) 
    beta12 ~ dunif(0,1)
    beta13 ~ dunif(0,1) 
    beta14 ~ dunif(0,1)
    beta15 ~ dunif(0,1) 
    beta16 ~ dunif(0,1)
    
    beta0_gamA  <- logit(beta1)
    beta0_gamAB <- logit(beta2)
    beta0_gamB  <- logit(beta3)
    beta0_gamBA <- logit(beta4)
    
    beta0_epsA  <- logit(beta5)
    beta0_epsAB <- logit(beta6)
    beta0_epsB  <- logit(beta7)
    beta0_epsBA <- logit(beta8)
    
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
  
    # First season probabilities for each state
    #site    
    for(j in 1:nsite){
    fsm[j, b, 1] <- 1-psi[b,1]-psi[b,2]-psi[b,3]  #-----------|U
    fsm[j, b, 2] <- psi[b,1]                      #-----------|A
    fsm[j, b, 3] <- psi[b,2]                      #-----------|B
    fsm[j, b, 4] <- psi[b,3]                      #-----------|AB

    
    # first season latent state
    # for sites   
      z[j, b, 1] ~ dcat( fsm[j, b, ( 1:nout )])    # site 
    } # end site loop
    
    # for block, which is just a function of the states of the sites within each block
   
     x[b,1] <- ifelse(sum(z[,b,1]==1) == nsite, 1,
                ifelse(sum(z[,b,1]==2) + sum(z[,b,1]==1) == nsite, 2,
                  ifelse(sum(z[,b,1]==3) + sum(z[,b,1]==1) == nsite, 3, 4) ) )
    
    } # end block loop

    ###############################################
    # btpm = block transition probability matrix. #
    # All columns sum to 1.                       #
    ###############################################
    
    # U to ...
    btpm[ 1, 1] <- (1-GamA) * (1-GamB)    #--|U
    btpm[ 2, 1] <- GamA * (1-GamB)        #--|A
    btpm[ 3, 1] <- (1-GamA) * GamB        #--|B
    btpm[ 4, 1] <- GamA * GamB            #--|AB
    
    # A to ...
    btpm[ 1, 2] <- EpsA * (1-GamBA)       #--|U
    btpm[ 2, 2] <- (1-EpsA) * (1-GamBA)   #--|A
    btpm[ 3, 2] <- EpsA * GamBA           #--|B
    btpm[ 4, 2] <- (1-EpsA) * GamBA       #--|AB
    
    # B to ...
    btpm[ 1, 3] <- (1-GamAB) * EpsB       #--|U
    btpm[ 2, 3] <- GamAB * EpsB           #--|A
    btpm[ 3, 3] <- (1-GamAB) * (1-EpsB)   #--|B
    btpm[ 4, 3] <- GamAB * (1-EpsB)       #--|AB
    
    # AB to ..
    btpm[ 1, 4] <- EpsAB * EpsBA          #--|U
    btpm[ 2, 4] <- (1-EpsAB) * EpsBA      #--|A
    btpm[ 3, 4] <- EpsAB * (1-EpsBA)      #--|B
    btpm[ 4, 4] <- (1-EpsAB) * (1-EpsBA)  #--|AB


      
    ####################################################
    ## stpm = site transition probability matrix.     ##
    ## These are dependent on the block level state   ##
    ## All columns sum to 1.                          ##
    ####################################################
    
    # blocks state (x) = U
    
    # U to ...
    stpm[ 1, 1, 1] <- 1          #--|U
    stpm[ 2, 1, 1] <- 0          #--|A
    stpm[ 3, 1, 1] <- 0          #--|B
    stpm[ 4, 1, 1] <- 0          #--|AB
    
    # A to ...
    stpm[ 1, 2, 1] <- 1          #--|U
    stpm[ 2, 2, 1] <- 0          #--|A
    stpm[ 3, 2, 1] <- 0          #--|B
    stpm[ 4, 2, 1] <- 0          #--|AB
    
    # B to ...
    stpm[1, 3, 1] <- 1           #--|U
    stpm[2, 3, 1] <- 0           #--|A
    stpm[3, 3, 1] <- 0           #--|B
    stpm[4, 3, 1] <- 0           #--|AB
    
    # AB to ..
    stpm[1, 4, 1] <- 1           #--|U
    stpm[2, 4, 1] <- 0           #--|A
    stpm[3, 4, 1] <- 0           #--|B
    stpm[4, 4, 1] <- 0           #--|AB

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 # blocks state (x) = A
 
# U to ...
    stpm[1, 1, 2] <- (1-gamA)    #--|U
    stpm[2, 1, 2] <- gamA        #--|A
    stpm[3, 1, 2] <- 0           #--|B
    stpm[4, 1, 2] <- 0           #--|AB
    
    # A to ...
    stpm[1, 2, 2] <- epsA        #--|U
    stpm[2, 2, 2] <- (1-epsA)    #--|A
    stpm[3, 2, 2] <- 0           #--|B
    stpm[4, 2, 2] <- 0           #--|AB
    
    # B to ...
    stpm[1, 3, 2] <- (1-gamAB)   #--|U
    stpm[2, 3, 2] <- gamAB       #--|A
    stpm[3, 3, 2] <- 0           #--|B
    stpm[4, 3, 2] <- 0           #--|AB
    
    # AB to ..
    stpm[1, 4, 2] <- epsAB       #--|U
    stpm[2, 4, 2] <- (1-epsAB)   #--|A
    stpm[3, 4, 2] <- 0           #--|B
    stpm[4, 4, 2] <- 0           #--|AB

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 # blocks state (x) = B
 
# U to ...
    stpm[1, 1, 3] <-  (1-gamB)    #--|U
    stpm[2, 1, 3] <- 0            #--|A
    stpm[3, 1, 3] <-  gamB        #--|B
    stpm[4, 1, 3] <- 0            #--|AB
    
    # A to ...
    stpm[1, 2, 3] <- (1-gamBA)    #--|U
    stpm[2, 2, 3] <- 0            #--|A
    stpm[3, 2, 3] <- gamBA        #--|B
    stpm[4, 2, 3] <- 0            #--|AB
    
    # B to ...
    stpm[1, 3, 3] <-  epsB        #--|U
    stpm[2, 3, 3] <- 0            #--|A
    stpm[3, 3, 3] <- (1-epsB)     #--|B
    stpm[4, 3, 3] <- 0            #--|AB
    
    # AB to ..
    stpm[1, 4, 3] <-  epsBA       #--|U
    stpm[2, 4, 3] <- 0            #--|A
    stpm[3, 4, 3] <-  (1-epsBA)   #--|B
    stpm[4, 4, 3] <- 0            #--|AB

 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

 # blocks state (x) = AB
 
# U to ...
    stpm[1, 1, 4] <- (1-gamA) * (1-gamB)    #--|U
    stpm[2, 1, 4] <- gamA  * (1-gamB)       #--|A
    stpm[3, 1, 4] <- (1-gamA) * gamB        #--|B
    stpm[4, 1, 4] <- gamA * gamB            #--|AB
    
    # A to ...
    stpm[1, 2, 4] <- epsA * (1-gamBA)       #--|U
    stpm[2, 2, 4] <- (1-epsA) * (1-gamBA)   #--|A
    stpm[3, 2, 4] <- epsA * gamBA           #--|B
    stpm[4, 2, 4] <- (1-epsA) * gamBA       #--|AB
    
    # B to ...
    stpm[1, 3, 4] <- (1-gamAB ) * epsB      #--|U
    stpm[2, 3, 4] <- gamAB  * epsB          #--|A
    stpm[3, 3, 4] <- (1-gamAB ) * (1-epsB)  #--|B
    stpm[4, 3, 4] <- gamAB  * (1-epsB)      #--|AB
    
    # AB to ..
    stpm[1, 4, 4] <- epsAB * epsBA          #--|U
    stpm[2, 4, 4] <- (1-epsAB) * epsBA      #--|A
    stpm[3, 4, 4] <- epsAB * (1-epsBA)      #--|B
    stpm[4, 4, 4] <- (1-epsAB) * (1-epsBA)  #--|AB
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## Logit links for col and ext probabilities

    logit(gamA)  <- beta0_gamA 
    logit(gamAB) <- beta0_gamAB 
    logit(gamB)  <- beta0_gamB 
    logit(gamBA) <- beta0_gamBA 
 
    logit(epsA)  <- beta0_epsA 
    logit(epsAB) <- beta0_epsAB 
    logit(epsB)  <- beta0_epsB 
    logit(epsBA) <- beta0_epsBA 

    logit(GamA)  <- beta0_GamA 
    logit(GamAB) <- beta0_GamAB 
    logit(GamB)  <- beta0_GamB 
    logit(GamBA) <- beta0_GamBA 
 
    logit(EpsA)  <- beta0_EpsA 
    logit(EpsAB) <- beta0_EpsAB 
    logit(EpsB)  <- beta0_EpsB 
    logit(EpsBA) <- beta0_EpsBA
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
    # latent block state for the rest of the seasons
    for(t in 1:(nseason-1)){
      for(b in 1:nblock){
        x[b, t+1] ~ dcat(btpm[(1:nout), x[b, t]])

      # latent site state for the rest of the seasons
      for(j in 1:nsite){
        z[j, b, t+1] ~ dcat( stpm[( 1:nout ) , z[ j, b, t], x[b,t+1]] + 0.01)  # +0.01 to avoide giving the dcat a prob of 0 
       
          for(day in 1:nsurvey) {      
            y[j, b, t, day] ~ dcat( dpm[t, ( 1:nout ) , z[j, b, t]] + 0.01)  # +0.01 to avoide giving the dcat a prob of 0 
          } #end survey loop
        } # end site loop
      } # end block loop

    
    #############################################################
    ## detection matrix (OS = observed state, TS = true state) ##
    #############################################################
    
    # TS = U
    dpm[t, 1, 1] <- 1                      #--|OS = U
    dpm[t, 2, 1] <- 0                      #--|OS = A
    dpm[t, 3, 1] <- 0                      #--|OS = B
    dpm[t, 4, 1] <- 0                      #--|OS = AB
    
    # TS = A
    dpm[t, 1, 2] <- 1-pA[t]                #--|OS = U
    dpm[t, 2, 2] <- pA[t]                  #--|OS = A
    dpm[t, 3, 2] <- 0                      #--|OS = B
    dpm[t, 4, 2] <- 0                      #--|OS = AB
    
    # TS = B
    dpm[t, 1, 3] <- 1-pB[t]                #--|OS = U
    dpm[t, 2, 3] <- 0                      #--|OS = A
    dpm[t, 3, 3] <- pB[t]                  #--|OS = B
    dpm[t, 4, 3] <- 0                      #--|OS = AB
    
    # TS = AB
    dpm[t, 1, 4] <- (1-pA[t]) * (1-pB[t])  #--|OS = U
    dpm[t, 2, 4] <- pA[t] * (1-pB[t])      #--|OS = A
    dpm[t, 3, 4] <- (1-pA[t]) * pB[t]      #--|OS = B
    dpm[t, 4, 4] <- pA[t] * pB[t]          #--|OS = AB
    
    ## logit links for detection probs
    logit(pA[t]) <- alphaA0 
    logit(pB[t]) <- alphaB0 
    } #close time loop
    
    ## Derived parameters
    
    diff_gamA <- gamA - gamAB
    diff_gamB <- gamB - gamBA
    diff_epsA <- epsA - epsAB
    diff_epsB <- epsB - epsBA
    
    diff_GamA <- GamA - GamAB
    diff_GamB <- GamB - GamBA
    diff_EpsA <- EpsA - EpsAB
    diff_EpsB <- EpsB - EpsBA
    
    ratio_gamA <- gamA / gamAB
    ratio_gamB <- gamB / gamBA
    ratio_epsA <- epsAB / epsA
    ratio_epsB <- epsB / epsBA
    
    ratio_GamA <- GamA / GamAB
    ratio_GamB <- GamBA / GamB
    ratio_EpsA <- EpsAB / EpsA
    ratio_EpsB <- EpsB / EpsBA    
    
    }# end
    ",fill = TRUE)
sink()

## import data
setwd("./data") # set wd to where the data is stored

load("occm_var.rda")    
#load("case_study_data.RData")
yb <-occm_va # change name of imported object to fit with the rest of the code

dim(yb) # check that dimensions are ok

# give data
data <-list(nseason = dim(yb)[3], nblock = dim(yb)[2], nsite = dim(yb)[1], nsurvey = dim(yb)[4], 
            nout=4, y = yb)

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
  beta1=runif(1,0.1,0.9), beta2=runif(1,0.1,0.9), beta3=runif(1,0.1,0.9), beta4=runif(1,0.1,0.9),
  beta5=runif(1,0.1,0.9), beta6=runif(1,0.1,0.9), beta7=runif(1,0.1,0.9), beta8=runif(1,0.1,0.9), 
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
            "beta0_EpsA", "beta0_EpsAB", "beta0_EpsB", "beta0_EpsBA",
            "diff_gamA", "diff_gamB", "diff_epsA", "diff_epsB", "diff_GamA", "diff_GamB", "diff_EpsA", "diff_EpsB",
            "ratio_gamA", "ratio_gamB", "ratio_epsA", "ratio_epsB", "ratio_GamA", "ratio_GamB", "ratio_EpsA", "ratio_EpsB")

# MCMC settings
ni <- 100   ;   nt <- 1   ;   nb <- 25 ;   nc <- 4    ;   na <- 25

# run model in jags
setwd("../")

va_mustela_rodent_hab <- jags(data, inits=inits, params, "mod_seas_det_4stpm.txt", n.chains = nc,
                          n.thin = nt, n.iter = ni, n.burnin = nb, n.adapt=na, parallel = T)

# Save model
setwd("./model_output")
save(va__mustela_rodent_hab, file="va_mustela_rodent_hab.rda")

#~ End of script