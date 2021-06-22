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
setwd("../model_output")

#load model
#setwd("./hidden_block_ko_fromautoclass/ko_vj/model_output")
dir()

#load("va_snowbed_mustela_rodent_sdet_s1_203_ni250k.rda")
load(dir()[2])
load(dir()[10])

##############################################################################

## plot simple model for rodents and mustelids

# traceplots
par(mfrow=c(2,4))
traceplot(va_mustela_rodent, parameters = c("GamA","GamB","GamAB","GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"))
# GamA, GamB and GamAB still have not converged properly. A longer chain might help.

traceplot(va_mustela_rodent, parameters = c("gamA","gamB","gamAB","gamBA", "epsA", "epsB", "epsAB", "epsBA"))
# all site parameters look fine


# create df in a format suitable for ggplot
dat <- data.frame(sims=unlist(va_mustela_rodent$sims.list[c(1:8,10:17)]),
                  par=c(rep("gamA", times=4000),rep("gamB", times=4000),rep("gamAB", times=4000),rep("gamBA", times=4000),
                        rep("epsA", times=4000),rep("epsB", times=4000),rep("epsAB", times=4000),rep("epsBA", times=4000),
                        rep("GamA", times=4000),rep("GamB", times=4000),rep("GamAB", times=4000),rep("GamBA", times=4000),
                        rep("EpsA", times=4000),rep("EpsB", times=4000),rep("EpsAB", times=4000),rep("EpsBA", times=4000)),
                  level=c(rep("Site parameters", times=4000*8), rep("Block parameters", times=4000*8)))


dat2 <- data.frame(sims=unlist(va_mustela_rodent$mean[c(1:8,10:17)]),
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
ggsave("modperf_va_mustela_rodent_col&ext.png", width = 60, height = 30, units="cm")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot habitat model

# traceplots
par(mfrow=c(2,4))
traceplot(va_mustela_rodent_hab_ni5k, parameters = c("GamA","GamB","GamAB","GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"))
traceplot(va_mustela_rodent_hab_ni5k, parameters = c("gamA","gamB","gamAB","gamBA", "epsA", "epsB", "epsAB", "epsBA"))

# create df in a format suitable for ggplot
dat <- data.frame(sims = c(va_mustela_rodent_hab_ni5k$sims.list["gamA"][[1]][,1,1], va_mustela_rodent_hab_ni5k$sims.list["gamA"][[1]][,1,11],
                           va_mustela_rodent_hab_ni5k$sims.list["gamB"][[1]][,1,1], va_mustela_rodent_hab_ni5k$sims.list["gamB"][[1]][,1,11],
                           va_mustela_rodent_hab_ni5k$sims.list["gamAB"][[1]][,1,1], va_mustela_rodent_hab_ni5k$sims.list["gamAB"][[1]][,1,11],
                           va_mustela_rodent_hab_ni5k$sims.list["gamBA"][[1]][,1,1], va_mustela_rodent_hab_ni5k$sims.list["gamBA"][[1]][,1,11],
                           va_mustela_rodent_hab_ni5k$sims.list["epsA"][[1]][,1,1], va_mustela_rodent_hab_ni5k$sims.list["epsA"][[1]][,1,11],
                           va_mustela_rodent_hab_ni5k$sims.list["epsB"][[1]][,1,1], va_mustela_rodent_hab_ni5k$sims.list["epsB"][[1]][,1,11],
                           va_mustela_rodent_hab_ni5k$sims.list["epsAB"][[1]][,1,1], va_mustela_rodent_hab_ni5k$sims.list["epsAB"][[1]][,1,11],
                           va_mustela_rodent_hab_ni5k$sims.list["epsBA"][[1]][,1,1], va_mustela_rodent_hab_ni5k$sims.list["epsBA"][[1]][,1,11],
                           va_mustela_rodent_hab_ni5k$sims.list["GamA"][[1]][,1], va_mustela_rodent_hab_ni5k$sims.list["GamB"][[1]][,1],
                           va_mustela_rodent_hab_ni5k$sims.list["GamAB"][[1]][,1], va_mustela_rodent_hab_ni5k$sims.list["GamBA"][[1]][,1],
                           va_mustela_rodent_hab_ni5k$sims.list["EpsA"][[1]][,1], va_mustela_rodent_hab_ni5k$sims.list["EpsB"][[1]][,1],
                           va_mustela_rodent_hab_ni5k$sims.list["EpsAB"][[1]][,1], va_mustela_rodent_hab_ni5k$sims.list["EpsBA"][[1]][,1]),
                  par=c(rep("gamA_1", times=2000),rep("gamA_2", times=2000),rep("gamB_1", times=2000),rep("gamB_2", times=2000),
                        rep("gamAB_1", times=2000),rep("gamAB_2", times=2000),rep("gamBA_1", times=2000), rep("gamBA_2", times=2000),
                        rep("epsA_1", times=2000),rep("epsA_2", times=2000), rep("epsB_1", times=2000), rep("epsB_2", times=2000),
                        rep("epsAB_1", times=2000),rep("epsAB_2", times=2000),rep("epsBA_1", times=2000),rep("epsBA_2", times=2000),
                        rep("GamA", times=2000),rep("GamB", times=2000),rep("GamAB", times=2000),rep("GamBA", times=2000),
                        rep("EpsA", times=2000),rep("EpsB", times=2000),rep("EpsAB", times=2000),rep("EpsBA", times=2000)),
                  level=c(rep("Site parameters", times=2000*8*2), rep("Block parameters", times=2000*8)))


dat2 <- data.frame(sims=c(va_mustela_rodent_hab_ni5k$mean["gamA"][[1]][1,1], va_mustela_rodent_hab_ni5k$mean["gamA"][[1]][1,11],
                          va_mustela_rodent_hab_ni5k$mean["gamB"][[1]][1,1], va_mustela_rodent_hab_ni5k$mean["gamB"][[1]][1,11],
                          va_mustela_rodent_hab_ni5k$mean["gamAB"][[1]][1,1], va_mustela_rodent_hab_ni5k$mean["gamAB"][[1]][1,11],
                          va_mustela_rodent_hab_ni5k$mean["gamBA"][[1]][1,1], va_mustela_rodent_hab_ni5k$mean["gamBA"][[1]][1,11],
                          va_mustela_rodent_hab_ni5k$mean["epsA"][[1]][1,1], va_mustela_rodent_hab_ni5k$mean["epsA"][[1]][1,11],
                          va_mustela_rodent_hab_ni5k$mean["epsB"][[1]][1,1], va_mustela_rodent_hab_ni5k$mean["epsB"][[1]][1,11],
                          va_mustela_rodent_hab_ni5k$mean["epsAB"][[1]][1,1], va_mustela_rodent_hab_ni5k$mean["epsAB"][[1]][1,11],
                          va_mustela_rodent_hab_ni5k$mean["epsBA"][[1]][1,1], va_mustela_rodent_hab_ni5k$mean["epsBA"][[1]][1,11],
                          va_mustela_rodent_hab_ni5k$mean["GamA"][[1]][1], va_mustela_rodent_hab_ni5k$mean["GamB"][[1]][1],
                          va_mustela_rodent_hab_ni5k$mean["GamAB"][[1]][1], va_mustela_rodent_hab_ni5k$mean["GamBA"][[1]][1],
                          va_mustela_rodent_hab_ni5k$mean["EpsA"][[1]][1], va_mustela_rodent_hab_ni5k$mean["EpsB"][[1]][1],
                          va_mustela_rodent_hab_ni5k$mean["EpsAB"][[1]][1], va_mustela_rodent_hab_ni5k$mean["EpsBA"][[1]][1]),
                   par=c("gamA_1", "gamA_2", "gamB_1", "gamB_2", "gamAB_1", "gamAB_2", "gamBA_1", "gamBA_2",
                         "epsA_1", "epsA_2", "epsB_1", "epsB_2", "epsAB_1", "epsAB_2", "epsBA_1", "epsBA_2",
                         "GamA", "GamB", "GamAB", "GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"),
                   level=c(rep("Site parameters", times=8*2), rep("Block parameters", times=8)))

# plot colonization and extinction probabilities
ggplot(data=dat, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  geom_point(data=dat2, color="red", shape="-", size=20)+
  labs(y="", x="")+
  facet_wrap(~level, scales="free_x")+
   theme(strip.text = element_text(size=30,face="bold"),
         axis.text.x = element_text(size=30, face="bold"),
         axis.text.y = element_text(size=30, face="bold"), legend.position = "none")+
  scale_x_discrete( labels=c('gamA_1' = parse(text = TeX("$\\gamma_{A1}$")),'gamA_2' = parse(text = TeX("$\\gamma_{A2}$")),
                             'gamB_1' = parse(text = TeX("$\\gamma_{B1}$")),'gamB_2' = parse(text = TeX("$\\gamma_{B2}$")), 
                             'gamAB_1' = parse(text = TeX("$\\gamma_{AB1}$")),'gamAB_2' = parse(text = TeX("$\\gamma_{AB2}$")),
                             'gamBA_1' = parse(text = TeX("$\\gamma_{BA1}$")),'gamBA_2' = parse(text = TeX("$\\gamma_{BA2}$")),
                             'epsA_1' = parse(text = TeX("$\\epsilon_{A1}$")),'epsA_2' = parse(text = TeX("$\\epsilon_{A2}$")),
                             'epsB_1' = parse(text = TeX("$\\epsilon_{B1}$")),'epsB_2' = parse(text = TeX("$\\epsilon_{B2}$")),
                             'epsAB_1' = parse(text = TeX("$\\epsilon_{AB1}$")),'epsAB_2' = parse(text = TeX("$\\epsilon_{AB2}$")),
                             'epsBA_1' = parse(text = TeX("$\\epsilon_{BA1}$")),'epsBA_2' = parse(text = TeX("$\\epsilon_{BA2}$")),
                             'GamA' = parse(text = TeX("$\\Gamma_{A}$")), 'GamB' = parse(text = TeX("$\\Gamma_{B}$")),
                             'GamAB' = parse(text = TeX("$\\Gamma_{AB}$")), 'GamBA' = parse(text = TeX("$\\Gamma_{BA}$")),
                             'EpsA' = parse(text = TeX("$E_{A}$")), 'EpsB' = parse(text = TeX("$E_{B}$")), 
                             'EpsAB' = parse(text = TeX("$E_{AB}$")), 'EpsBA' = parse(text = TeX("$E_{BA}$")) ))

# save plot
setwd("../plot")
ggsave("modperf_va_mustela_rodent_hab_2.png", width = 90, height = 30, units="cm")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot snow model

# traceplots
par(mfrow=c(4,4))
traceplot(va_mustela_rodent_temp_ni1k, parameters = c("GamA","GamB","GamAB","GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"))
#traceplot(va_mustela_rodent_temp_ni1k, parameters = c("gamA","gamB","gamAB","gamBA", "epsA", "epsB", "epsAB", "epsBA"))

# create df in a format suitable for ggplot
n <- 4000
dat <- data.frame(sims = c(va_mustela_rodent_temp_ni1k$sims.list["gamA"][[1]][,1,1,1], va_mustela_rodent_temp_ni1k$sims.list["gamA"][[1]][,1,1,70],
                           va_mustela_rodent_temp_ni1k$sims.list["gamB"][[1]][,1,1,1], va_mustela_rodent_temp_ni1k$sims.list["gamB"][[1]][,1,1,70],
                           va_mustela_rodent_temp_ni1k$sims.list["gamAB"][[1]][,1,1,1], va_mustela_rodent_temp_ni1k$sims.list["gamAB"][[1]][,1,1,70],
                           va_mustela_rodent_temp_ni1k$sims.list["gamBA"][[1]][,1,1,1], va_mustela_rodent_temp_ni1k$sims.list["gamBA"][[1]][,1,1,70],
                           va_mustela_rodent_temp_ni1k$sims.list["epsA"][[1]][,1,1,1], va_mustela_rodent_temp_ni1k$sims.list["epsA"][[1]][,1,1,70],
                           va_mustela_rodent_temp_ni1k$sims.list["epsB"][[1]][,1,1,1], va_mustela_rodent_temp_ni1k$sims.list["epsB"][[1]][,1,1,70],
                           va_mustela_rodent_temp_ni1k$sims.list["epsAB"][[1]][,1,1,1], va_mustela_rodent_temp_ni1k$sims.list["epsAB"][[1]][,1,1,70],
                           va_mustela_rodent_temp_ni1k$sims.list["epsBA"][[1]][,1,1,1], va_mustela_rodent_temp_ni1k$sims.list["epsBA"][[1]][,1,1,70],
                           va_mustela_rodent_temp_ni1k$sims.list["GamA"][[1]][,1], va_mustela_rodent_temp_ni1k$sims.list["GamB"][[1]][,1],
                           va_mustela_rodent_temp_ni1k$sims.list["GamAB"][[1]][,1], va_mustela_rodent_temp_ni1k$sims.list["GamBA"][[1]][,1],
                           va_mustela_rodent_temp_ni1k$sims.list["EpsA"][[1]][,1], va_mustela_rodent_temp_ni1k$sims.list["EpsB"][[1]][,1],
                           va_mustela_rodent_temp_ni1k$sims.list["EpsAB"][[1]][,1], va_mustela_rodent_temp_ni1k$sims.list["EpsBA"][[1]][,1]),
                  par=c(rep("gamA_1", times=n),rep("gamA_2", times=n),rep("gamB_1", times=n),rep("gamB_2", times=n),
                        rep("gamAB_1", times=n),rep("gamAB_2", times=n),rep("gamBA_1", times=n), rep("gamBA_2", times=n),
                        rep("epsA_1", times=n),rep("epsA_2", times=n), rep("epsB_1", times=n), rep("epsB_2", times=n),
                        rep("epsAB_1", times=n),rep("epsAB_2", times=n),rep("epsBA_1", times=n),rep("epsBA_2", times=n),
                        rep("GamA", times=n),rep("GamB", times=n),rep("GamAB", times=n),rep("GamBA", times=n),
                        rep("EpsA", times=n),rep("EpsB", times=n),rep("EpsAB", times=n),rep("EpsBA", times=n)),
                  level=c(rep("Site parameters", times=n*8*2), rep("Block parameters", times=n*8)))


dat2 <- data.frame(sims=c(va_mustela_rodent_temp_ni1k$mean["gamA"][[1]][1,1,1], va_mustela_rodent_temp_ni1k$mean["gamA"][[1]][1,1,70],
                          va_mustela_rodent_temp_ni1k$mean["gamB"][[1]][1,1,1], va_mustela_rodent_temp_ni1k$mean["gamB"][[1]][1,1,70],
                          va_mustela_rodent_temp_ni1k$mean["gamAB"][[1]][1,1,1], va_mustela_rodent_temp_ni1k$mean["gamAB"][[1]][1,1,70],
                          va_mustela_rodent_temp_ni1k$mean["gamBA"][[1]][1,1,1], va_mustela_rodent_temp_ni1k$mean["gamBA"][[1]][1,1,70],
                          va_mustela_rodent_temp_ni1k$mean["epsA"][[1]][1,1,1], va_mustela_rodent_temp_ni1k$mean["epsA"][[1]][1,1,70],
                          va_mustela_rodent_temp_ni1k$mean["epsB"][[1]][1,1,1], va_mustela_rodent_temp_ni1k$mean["epsB"][[1]][1,1,70],
                          va_mustela_rodent_temp_ni1k$mean["epsAB"][[1]][1,1,1], va_mustela_rodent_temp_ni1k$mean["epsAB"][[1]][1,1,70],
                          va_mustela_rodent_temp_ni1k$mean["epsBA"][[1]][1,1,1], va_mustela_rodent_temp_ni1k$mean["epsBA"][[1]][1,1,70],
                          va_mustela_rodent_temp_ni1k$mean["GamA"][[1]][1], va_mustela_rodent_temp_ni1k$mean["GamB"][[1]][1],
                          va_mustela_rodent_temp_ni1k$mean["GamAB"][[1]][1], va_mustela_rodent_temp_ni1k$mean["GamBA"][[1]][1],
                          va_mustela_rodent_temp_ni1k$mean["EpsA"][[1]][1], va_mustela_rodent_temp_ni1k$mean["EpsB"][[1]][1],
                          va_mustela_rodent_temp_ni1k$mean["EpsAB"][[1]][1], va_mustela_rodent_temp_ni1k$mean["EpsBA"][[1]][1]),
                   par=c("gamA_1", "gamA_2", "gamB_1", "gamB_2", "gamAB_1", "gamAB_2", "gamBA_1", "gamBA_2",
                         "epsA_1", "epsA_2", "epsB_1", "epsB_2", "epsAB_1", "epsAB_2", "epsBA_1", "epsBA_2",
                         "GamA", "GamB", "GamAB", "GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"),
                   level=c(rep("Site parameters", times=8*2), rep("Block parameters", times=8)))

# plot colonization and extinction probabilities
ggplot(data=dat, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  geom_point(data=dat2, color="red", shape="-", size=20)+
  labs(y="", x="")+
  facet_wrap(~level, scales="free_x")+
  theme(strip.text = element_text(size=30,face="bold"),
        axis.text.x = element_text(size=30, face="bold"),
        axis.text.y = element_text(size=30, face="bold"), legend.position = "none")+
  scale_x_discrete( labels=c('gamA_1' = parse(text = TeX("$\\gamma_{A1}$")),'gamA_2' = parse(text = TeX("$\\gamma_{A2}$")),
                             'gamB_1' = parse(text = TeX("$\\gamma_{B1}$")),'gamB_2' = parse(text = TeX("$\\gamma_{B2}$")), 
                             'gamAB_1' = parse(text = TeX("$\\gamma_{AB1}$")),'gamAB_2' = parse(text = TeX("$\\gamma_{AB2}$")),
                             'gamBA_1' = parse(text = TeX("$\\gamma_{BA1}$")),'gamBA_2' = parse(text = TeX("$\\gamma_{BA2}$")),
                             'epsA_1' = parse(text = TeX("$\\epsilon_{A1}$")),'epsA_2' = parse(text = TeX("$\\epsilon_{A2}$")),
                             'epsB_1' = parse(text = TeX("$\\epsilon_{B1}$")),'epsB_2' = parse(text = TeX("$\\epsilon_{B2}$")),
                             'epsAB_1' = parse(text = TeX("$\\epsilon_{AB1}$")),'epsAB_2' = parse(text = TeX("$\\epsilon_{AB2}$")),
                             'epsBA_1' = parse(text = TeX("$\\epsilon_{BA1}$")),'epsBA_2' = parse(text = TeX("$\\epsilon_{BA2}$")),
                             'GamA' = parse(text = TeX("$\\Gamma_{A}$")), 'GamB' = parse(text = TeX("$\\Gamma_{B}$")),
                             'GamAB' = parse(text = TeX("$\\Gamma_{AB}$")), 'GamBA' = parse(text = TeX("$\\Gamma_{BA}$")),
                             'EpsA' = parse(text = TeX("$E_{A}$")), 'EpsB' = parse(text = TeX("$E_{B}$")), 
                             'EpsAB' = parse(text = TeX("$E_{AB}$")), 'EpsBA' = parse(text = TeX("$E_{BA}$")) ))

# save plot
setwd("../plot")
ggsave("modperf_va_mustela_rodent_snow.png", width = 90, height = 30, units="cm")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Plot hab and snow model

# traceplots
par(mfrow=c(2,4))
traceplot(va_mustela_rodent_habtemp_ni5k, parameters = c("GamA","GamB","GamAB","GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"))
#traceplot(va_mustela_rodent_temp_ni1k, parameters = c("gamA","gamB","gamAB","gamBA", "epsA", "epsB", "epsAB", "epsBA"))

# create df in a format suitable for ggplot
n <- 2000

mod1 <- va_mustela_rodent_habtemp_ni5k

dat <- data.frame(sims = c(mod1$sims.list["gamA"][[1]][,1,1,1],  mod1$sims.list["gamA"][[1]][,1,1,70],
                           mod1$sims.list["gamA"][[1]][,1,11,1],  mod1$sims.list["gamA"][[1]][,1,11,70],
                           mod1$sims.list["gamB"][[1]][,1,1,1],  mod1$sims.list["gamB"][[1]][,1,1,70],
                           mod1$sims.list["gamB"][[1]][,1,11,1],  mod1$sims.list["gamB"][[1]][,1,11,70],
                           mod1$sims.list["gamAB"][[1]][,1,1,1], mod1$sims.list["gamAB"][[1]][,1,1,70],
                           mod1$sims.list["gamAB"][[1]][,1,11,1], mod1$sims.list["gamAB"][[1]][,1,11,70],
                           mod1$sims.list["gamBA"][[1]][,1,1,1], mod1$sims.list["gamBA"][[1]][,1,1,70],
                           mod1$sims.list["gamBA"][[1]][,1,11,1], mod1$sims.list["gamBA"][[1]][,1,11,70],
                           mod1$sims.list["epsA"][[1]][,1,1,1],  mod1$sims.list["epsA"][[1]][,1,1,70],
                           mod1$sims.list["epsA"][[1]][,1,11,1],  mod1$sims.list["epsA"][[1]][,1,11,70],
                           mod1$sims.list["epsB"][[1]][,1,1,1],  mod1$sims.list["epsB"][[1]][,1,1,70],
                           mod1$sims.list["epsB"][[1]][,1,11,1],  mod1$sims.list["epsB"][[1]][,1,11,70],
                           mod1$sims.list["epsAB"][[1]][,1,1,1], mod1$sims.list["epsAB"][[1]][,1,1,70],
                           mod1$sims.list["epsAB"][[1]][,1,11,1], mod1$sims.list["epsAB"][[1]][,1,11,70],
                           mod1$sims.list["epsBA"][[1]][,1,1,1], mod1$sims.list["epsBA"][[1]][,1,1,70],
                           mod1$sims.list["epsBA"][[1]][,1,11,1], mod1$sims.list["epsBA"][[1]][,1,11,70],
                           mod1$sims.list["GamA"][[1]][,1],      mod1$sims.list["GamB"][[1]][,1],
                           mod1$sims.list["GamAB"][[1]][,1],     mod1$sims.list["GamBA"][[1]][,1],
                           mod1$sims.list["EpsA"][[1]][,1],      mod1$sims.list["EpsB"][[1]][,1],
                           mod1$sims.list["EpsAB"][[1]][,1],     mod1$sims.list["EpsBA"][[1]][,1]),
                  
                  par=c(rep("gamA_S_H", times=n),rep("gamA_S_S", times=n), rep("gamA_W_H", times=n),rep("gamA_W_S", times=n),
                        rep("gamB_S_H", times=n),rep("gamB_S_S", times=n),rep("gamB_W_H", times=n),rep("gamB_W_S", times=n),
                        rep("gamAB_S_H", times=n),rep("gamAB_S_S", times=n),rep("gamAB_W_H", times=n),rep("gamAB_W_S", times=n),
                        rep("gamBA_S_H", times=n), rep("gamBA_S_S", times=n), rep("gamBA_W_H", times=n), rep("gamBA_W_S", times=n),
                        rep("epsA_S_H", times=n),rep("epsA_S_S", times=n), rep("epsA_W_H", times=n),rep("epsA_W_S", times=n),
                        rep("epsB_S_H", times=n), rep("epsB_S_S", times=n),rep("epsB_W_H", times=n), rep("epsB_W_S", times=n),
                        rep("epsAB_S_H", times=n),rep("epsAB_S_S", times=n),rep("epsAB_W_H", times=n),rep("epsAB_W_S", times=n),
                        rep("epsBA_S_H", times=n),rep("epsBA_S_S", times=n), rep("epsBA_W_H", times=n),rep("epsBA_W_S", times=n),
                        rep("GamA", times=n),rep("GamB", times=n),rep("GamAB", times=n),rep("GamBA", times=n),
                        rep("EpsA", times=n),rep("EpsB", times=n),rep("EpsAB", times=n),rep("EpsBA", times=n)),
                  
                  level=c(rep("Site parameters", times=n*8*4), rep("Block parameters", times=n*8)))


dat2 <- data.frame(sims=c(mod1$mean["gamA"][[1]][1,1,1],  mod1$mean["gamA"][[1]][1,1,70],
                          mod1$mean["gamA"][[1]][1,11,1],  mod1$mean["gamA"][[1]][1,11,70],
                          mod1$mean["gamB"][[1]][1,1,1],  mod1$mean["gamB"][[1]][1,1,70],
                          mod1$mean["gamB"][[1]][1,11,1],  mod1$mean["gamB"][[1]][1,11,70],
                          mod1$mean["gamAB"][[1]][1,1,1], mod1$mean["gamAB"][[1]][1,1,70],
                          mod1$mean["gamAB"][[1]][1,11,1], mod1$mean["gamAB"][[1]][1,11,70],
                          mod1$mean["gamBA"][[1]][1,1,1], mod1$mean["gamBA"][[1]][1,1,70],
                          mod1$mean["gamBA"][[1]][1,11,1], mod1$mean["gamBA"][[1]][1,11,70],
                          mod1$mean["epsA"][[1]][1,1,1],  mod1$mean["epsA"][[1]][1,1,70],
                          mod1$mean["epsA"][[1]][1,11,1],  mod1$mean["epsA"][[1]][1,11,70],
                          mod1$mean["epsB"][[1]][1,1,1],  mod1$mean["epsB"][[1]][1,1,70],
                          mod1$mean["epsB"][[1]][1,11,1],  mod1$mean["epsB"][[1]][1,11,70],
                          mod1$mean["epsAB"][[1]][1,1,1], mod1$mean["epsAB"][[1]][1,1,70],
                          mod1$mean["epsAB"][[1]][1,11,1], mod1$mean["epsAB"][[1]][1,11,70],
                          mod1$mean["epsBA"][[1]][1,1,1], mod1$mean["epsBA"][[1]][1,1,70],
                          mod1$mean["epsBA"][[1]][1,11,1], mod1$mean["epsBA"][[1]][1,11,70],
                          mod1$mean["GamA"][[1]][1],      mod1$mean["GamB"][[1]][1],
                          mod1$mean["GamAB"][[1]][1],     mod1$mean["GamBA"][[1]][1],
                          mod1$mean["EpsA"][[1]][1],      mod1$mean["EpsB"][[1]][1],
                          mod1$mean["EpsAB"][[1]][1],     mod1$mean["EpsBA"][[1]][1]),
                   par=c("gamA_S_H", "gamA_S_S","gamA_W_H", "gamA_W_S",
                         "gamB_S_H", "gamB_S_S","gamB_W_H", "gamB_W_S",
                         "gamAB_S_H", "gamAB_S_S","gamAB_W_H", "gamAB_W_S",
                         "gamBA_S_H", "gamBA_S_S","gamBA_W_H", "gamBA_W_S",
                         "epsA_S_H", "epsA_S_S","epsA_W_H", "epsA_W_S",
                         "epsB_S_H", "epsB_S_S","epsB_W_H", "epsB_W_S",
                         "epsAB_S_H", "epsAB_S_S","epsAB_W_H", "epsAB_W_S",
                         "epsBA_S_H", "epsBA_S_S","epsBA_W_H", "epsBA_W_S",
                         "GamA", "GamB", "GamAB", "GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"),
                   level=c(rep("Site parameters", times=8*4), rep("Block parameters", times=8)))

# plot colonization and extinction probabilities
ggplot(data=dat, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  geom_point(data=dat2, color="red", shape="-", size=20)+
  labs(y="", x="")+
  facet_wrap(~level, scales="free_x")+
  theme(strip.text = element_text(size=30,face="bold"),
        axis.text.x = element_text(size=30, face="bold"),
        axis.text.y = element_text(size=30, face="bold"), legend.position = "none")+
  scale_x_discrete( labels=c('gamA_S_H' = parse(text = TeX("$\\gamma_{A_{SH}}$")),'gamA_S_S' = parse(text = TeX("$\\gamma_{A_{SS}}$")),
                             'gamA_W_H' = parse(text = TeX("$\\gamma_{A_{WH}}$")),'gamA_W_S' = parse(text = TeX("$\\gamma_{A_{WS}}$")),
                             'gamB_S_H' = parse(text = TeX("$\\gamma_{B_{S_H}}$")),'gamB_S_S' = parse(text = TeX("$\\gamma_{B_{SH}}$")),
                             'gamB_W_H' = parse(text = TeX("$\\gamma_{B_{W_H}}$")),'gamB_W_S' = parse(text = TeX("$\\gamma_{B_{WH}}$")),
                             'gamAB_S_H' = parse(text = TeX("$\\gamma_{AB_{SH}}$")),'gamAB_S_S' = parse(text = TeX("$\\gamma_{AB_{SS}}$")),
                             'gamAB_W_H' = parse(text = TeX("$\\gamma_{AB_{WH}}$")),'gamAB_W_S' = parse(text = TeX("$\\gamma_{AB_{WS}}$")),
                             'gamBA_S_H' = parse(text = TeX("$\\gamma_{BA_{SH}}$")),'gamBA_S_S' = parse(text = TeX("$\\gamma_{BA_{SS}}$")),
                             'gamBA_W_H' = parse(text = TeX("$\\gamma_{BA_{WH}}$")),'gamBA_W_S' = parse(text = TeX("$\\gamma_{BA_{WS}}$")),
                             'epsA_S_H' = parse(text = TeX("$\\epsilon_{A_{SH}}$")),'epsA_S_S' = parse(text = TeX("$\\epsilon_{A_{SS}}$")),
                             'epsA_W_H' = parse(text = TeX("$\\epsilon_{A_{WH}}$")),'epsA_W_S' = parse(text = TeX("$\\epsilon_{A_{WS}}$")),
                             'epsB_S_H' = parse(text = TeX("$\\epsilon_{B_{SH}}$")),'epsB_S_S' = parse(text = TeX("$\\epsilon_{B_{SS}}$")),
                             'epsB_W_H' = parse(text = TeX("$\\epsilon_{B_{WH}}$")),'epsB_W_S' = parse(text = TeX("$\\epsilon_{B_{WS}}$")),
                             'epsAB_S_H' = parse(text = TeX("$\\epsilon_{AB_{SH}}$")),'epsAB_S_S' = parse(text = TeX("$\\epsilon_{AB_{SS}}$")),
                             'epsAB_W_H' = parse(text = TeX("$\\epsilon_{AB_{WH}}$")),'epsAB_W_S' = parse(text = TeX("$\\epsilon_{AB_{WS}}$")),
                             'epsBA_S_H' = parse(text = TeX("$\\epsilon_{BA_{S_H}}$")),'epsBA_S_S' = parse(text = TeX("$\\epsilon_{BA_{SS}}$")),
                             'epsBA_W_H' = parse(text = TeX("$\\epsilon_{BA_{W_H}}$")),'epsBA_W_S' = parse(text = TeX("$\\epsilon_{BA_{WS}}$")),
                             'GamA' = parse(text = TeX("$\\Gamma_{A}$")), 'GamB' = parse(text = TeX("$\\Gamma_{B}$")),
                             'GamAB' = parse(text = TeX("$\\Gamma_{AB}$")), 'GamBA' = parse(text = TeX("$\\Gamma_{BA}$")),
                             'EpsA' = parse(text = TeX("$E_{A}$")), 'EpsB' = parse(text = TeX("$E_{B}$")), 
                             'EpsAB' = parse(text = TeX("$E_{AB}$")), 'EpsBA' = parse(text = TeX("$E_{BA}$")) ))

# save plot
setwd("../plot")
ggsave("modperf_va_mustela_rodent_snow&hab.png", width = 90, height = 30, units="cm")





#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## plot simple model for voles and mustelids

# traceplots
par(mfrow=c(2,4))
traceplot(va_mustela_vole, parameters = c("GamA","GamB","GamAB","GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"))
# GamA, GamB and GamAB still have not converged properly. A longer chain might help.

traceplot(va_mustela_vole, parameters = c("gamA","gamB","gamAB","gamBA", "epsA", "epsB", "epsAB", "epsBA"))
# all site parameters look fine


# create df in a format suitable for ggplot
n=4000
mod1 <- va_mustela_vole

dat <- data.frame(sims=unlist(mod1$sims.list[c(1:8,10:17)]),
                  par=c(rep("gamA", times=n),rep("gamB", times=n),rep("gamAB", times=n),rep("gamBA", times=n),
                        rep("epsA", times=n),rep("epsB", times=n),rep("epsAB", times=n),rep("epsBA", times=n),
                        rep("GamA", times=n),rep("GamB", times=n),rep("GamAB", times=n),rep("GamBA", times=n),
                        rep("EpsA", times=n),rep("EpsB", times=n),rep("EpsAB", times=n),rep("EpsBA", times=n)),
                  level=c(rep("Site parameters", times=n*8), rep("Block parameters", times=n*8)))


dat2 <- data.frame(sims=unlist(mod1$mean[c(1:8,10:17)]),
                   par=c("gamA", "gamB", "gamAB", "gamBA", "epsA", "epsB", "epsAB", "epsBA",
                         "GamA", "GamB", "GamAB", "GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"),
                   level=c(rep("Site parameters", times=8), rep("Block parameters", times=8)))

# plot colonization and extinction probabilities
plotvole <- ggplot(data=dat, aes(x=par, y=sims, fill="grey"))+
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


plotvole <- ggplot(data=dat, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  geom_point(data=dat2, color="red", shape="-", size=20)+
  labs(y="", x="")+
  facet_wrap(~level, scales="free_x")+
  theme(strip.text = element_text(size=30,face="bold"),
        axis.text.x = element_blank(),
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
ggsave("modperf_va_mustela_vole_col&ext.png", width = 60, height = 30, units="cm")



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## plot simple model for lemmings and mustelids

# traceplots
par(mfrow=c(2,4))
traceplot(va_mustela_lem, parameters = c("GamA","GamB","GamAB","GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"))
# GamA, GamB and GamAB still have not converged properly. A longer chain might help.

traceplot(va_mustela_lem, parameters = c("gamA","gamB","gamAB","gamBA", "epsA", "epsB", "epsAB", "epsBA"))
# all site parameters look fine


# create df in a format suitable for ggplot
n=4000
mod2 <- va_mustela_lem

dat <- data.frame(sims=unlist(mod2$sims.list[c(1:8,10:17)]),
                  par=c(rep("gamA", times=n),rep("gamB", times=n),rep("gamAB", times=n),rep("gamBA", times=n),
                        rep("epsA", times=n),rep("epsB", times=n),rep("epsAB", times=n),rep("epsBA", times=n),
                        rep("GamA", times=n),rep("GamB", times=n),rep("GamAB", times=n),rep("GamBA", times=n),
                        rep("EpsA", times=n),rep("EpsB", times=n),rep("EpsAB", times=n),rep("EpsBA", times=n)),
                  level=c(rep("Site parameters", times=n*8), rep("Block parameters", times=n*8)))


dat2 <- data.frame(sims=unlist(mod2$mean[c(1:8,10:17)]),
                   par=c("gamA", "gamB", "gamAB", "gamBA", "epsA", "epsB", "epsAB", "epsBA",
                         "GamA", "GamB", "GamAB", "GamBA", "EpsA", "EpsB", "EpsAB", "EpsBA"),
                   level=c(rep("Site parameters", times=8), rep("Block parameters", times=8)))

# plot colonization and extinction probabilities
plotlem <- ggplot(data=dat, aes(x=par, y=sims, fill="grey"))+
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


plotlem <- ggplot(data=dat, aes(x=par, y=sims, fill="grey"))+
  geom_violin(fill="grey")+
  geom_point(data=dat2, color="red", shape="-", size=20)+
  labs(y="", x="")+
  facet_wrap(~level, scales="free_x")+
  theme(strip.text = element_blank(),
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
ggsave("modperf_va_mustela_lem_col&ext.png", width = 60, height = 30, units="cm")



cowplot::plot_grid(plotvole, plotlem, labels = c("Vole","Lemming"),label_x = 0.06, label_y=0.9, label_size = 22 , nrow=2)
ggsave("modperf_va_mustela_lem&vole_col&ext.png", width = 60, height = 30, units="cm")

#ggpubr::ggarrange(plotvole, plotlem, labels = c("Vole","Lemming"), nrow=2)



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

#~ End of Script