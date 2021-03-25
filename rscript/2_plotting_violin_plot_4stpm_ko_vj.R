####################################################################################################################
##   Plotting the model output from the dynamic occupancy model for interacting species with two spatial scales   ##
##   analyzing the small mammal data from the Varanger peninsula                                                  ##
##                         by EFK and FB                                                                          ##
####################################################################################################################

# empty environment
rm(list=ls())

# Call jags(and other packages)
library(jagsUI)
library(ggplot2)
library(latex2exp)

# set working directory
setwd("~/UiT/Manuskript/TeoreticalModelingOfSmallRodents&Mustelids/OccupancyModel/models")

#load model
setwd("./hidden_block_ko_fromautoclass/ko_vj/model_output")
dir()

#load("va_snowbed_mustela_rodent_sdet_s1_203_ni250k.rda")
load(dir()[12])
##############################################################################

# traceplots
par(mfrow=c(2,4))
traceplot(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3, parameters = c("GamA","GamB","GamAB","GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"))
traceplot(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3, parameters = c("gamA","gamB","gamAB","gamBA", "epsA", "epsB", "epsAB", "epsBA"))

# create df in a format suitable for ggplot
dat <- data.frame(sims=unlist(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$sims.list[c(1:8,10:17)]),
                  par=c(rep("gamA", times=15000),rep("gamB", times=15000),rep("gamAB", times=15000),rep("gamBA", times=15000),
                        rep("epsA", times=15000),rep("epsB", times=15000),rep("epsAB", times=15000),rep("epsBA", times=15000),
                        rep("GamA", times=15000),rep("GamB", times=15000),rep("GamAB", times=15000),rep("GamBA", times=15000),
                        rep("EpsA", times=15000),rep("EpsB", times=15000),rep("EpsAB", times=15000),rep("EpsBA", times=15000)),
                  level=c(rep("Site parameters", times=15000*8), rep("Block parameters", times=15000*8)))


dat2 <- data.frame(sims=unlist(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$mean[c(1:8,10:17)]),
                   par=c("gamA", "gamB", "gamAB", "gamBA", "epsA", "epsB", "epsAB", "epsBA",
                         "GamA", "GamB", "GamAB", "GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"),
                   level=c(rep("Site parameters", times=8), rep("Block parameters", times=8)))

# plot colonization and extinction probabilities
ggplot(data=dat, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  geom_point(data=dat2, color="red", shape="-", size=20)+
  labs(y="", x="")+
  facet_wrap(~level, scales="free_x")+
   theme(strip.text = element_text(size=30,face="bold"),
         axis.text.x = element_text(size=40, face="bold"),
         axis.text.y = element_text(size=30, face="bold"), legend.position = "none")+
  scale_x_discrete( labels=c('gamA' = parse(text = TeX("$\\gamma_{A}$")), 'gamB' = parse(text = TeX("$\\gamma_{B}$")), 
                                        'gamAB' = parse(text = TeX("$\\gamma_{AB}$")), 'gamBA' = parse(text = TeX("$\\gamma_{BA}$")),
                                        'epsA' = parse(text = TeX("$\\epsilon_{A}$")), 'epsB' = parse(text = TeX("$\\epsilon_{B}$")),
                                        'epsAB' = parse(text = TeX("$\\epsilon_{AB}$")), 'epsBA' = parse(text = TeX("$\\epsilon_{BA}$")),
                                        'GamA' = parse(text = TeX("$\\Gamma_{A}$")), 'GamB' = parse(text = TeX("$\\Gamma_{B}$")),
                                        'GamAB' = parse(text = TeX("$\\Gamma_{AB}$")), 'GamBA' = parse(text = TeX("$\\Gamma_{BA}$")),
                                        'EpsA' = parse(text = TeX("$E_{A}$")), 'EpsB' = parse(text = TeX("$E_{B}$")), 
                                        'EpsAB' = parse(text = TeX("$E_{AB}$")), 'EpsBA' = parse(text = TeX("$E_{BA}$")) ))

# save plot
setwd("../plot")
ggsave("modperf_va_snowbed_mustela_rodent_sdet_4stpm_s1_203_ni100k_site&block_3.png", width = 60, height = 30, units="cm")

###########################################################################################################################
##    Plot estimated detection probabilities

# making df in a suitable format for ggplot
dat3 <- data.frame(sims=c(unlist(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$sims.list[18][[1]][,c(1,15)]),
                          unlist(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$sims.list[19][[1]][,c(1,15)])),
                   par=c(rep("pA_S", times=15000), rep("pA_W", times=15000), rep("pB_S", times=15000), rep("pB_W", times=15000)) )

dat4 <- data.frame(sims=c(unlist(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$mean[18][[1]][c(1,15)]),
                          unlist(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$mean[19][[1]][c(1,15)])),
                   par=c("pA_S", "pA_W", "pB_S", "pB_W"))

# specifing what order the parameters should be plotted
pos3 <-c("pA_S", "pA_W", "pB_S", "pB_W")  

#make plot
ggplot(data=dat3, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey", size=1.5)+
  geom_point(data=dat4, color="red", shape="-", size=20)+
  labs(y="", x="")+
  theme(axis.text.x = element_text(size=30, face="bold"),
        axis.text.y = element_text(size=30, face="bold"), legend.position = "none")+
  scale_x_discrete(limits=pos3, labels=c('pA_S' = parse(text = TeX("$p_{A_S}$")), 'pA_W' = parse(text = TeX("$p_{A_W}$")), 
                                         'pB_S' = parse(text = TeX("$p_{B_S}$")), 'pB_W' = parse(text = TeX("$p_{B_W}$"))))

# save plot
ggsave("modperf_va_snowbed_mustela_rodent_sdet_4stpm_s1_203_p_3.png", width = 60, height = 20, units="cm")


## plot estimated initial occupancy probability

# make df with a suitable structure for ggplot
dat5 <- data.frame(sims=c(as.vector(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$sims.list$psi[,,1]),
                          as.vector(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$sims.list$psi[,,2]),
                          as.vector(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$sims.list$psi[,,3])),
                   par=c(rep("psi1", times=120000), rep("psi2", times=120000), rep("psi3", times=120000)),
                   block=c(rep(c(rep("b1", times=15000),rep("b2",times=15000),rep("b3", times=15000),rep("b4", times=15000),rep("b5", times=15000),rep("b6", times=15000),rep("b7", times=15000),rep("b8", times=15000))
                               , times=3)))

# spesify what order to plot the simulations
pos5 <-c("psi1","psi2","psi3")  

# make plot
ggplot(data=dat5, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  geom_boxplot(data=dat5, width=0.4, color="black", alpha=0.4, fill="white")+
  geom_point(data=dat6, color="red", shape="-", size=15)+
  labs(y="", x="")+
  facet_wrap(~block, labeller = label_parsed, ncol=4)+
  theme(strip.text = element_text(size=30,face="bold"), axis.text.x = element_text(size=30, face="bold"),
        axis.text.y = element_text(size=30, face="bold"), legend.position = "none")+
  scale_x_discrete(limits=pos5, 
                   labels=c('psi1' = parse(text = TeX("$\\psi_{A}$")), 'psi2' = parse(text = TeX("$\\psi_{B}$")), 'psi3' = parse(text = TeX("$\\psi_{AB}$"))))

# save the model
ggsave("modperf_va_snowbed_mustela_rodent_sdet_4stpm_s1_203_psi_3.png", width = 60, height = 30, units="cm")

####
str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$mean[1:8])
str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$sd[1:8])
str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$q2.5[1:8])
str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$q97.5[1:8])

str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$mean[10:17])
str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$sd[10:17])
str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$q2.5[10:17])
str(va_snowbed_mustela_rodent_sdet_4stpm_ni100k_3$q97.5[10:17])
##

#~ End of Script