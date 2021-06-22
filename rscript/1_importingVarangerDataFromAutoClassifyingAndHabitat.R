# read libraries
library(lubridate)
library(dplyr)
library(tidyr)

# set working directory
#wd <- "C:/Eivind/GitProjects/MustelidsAndRodents-"
#setwd(wd)

# look at files in the directory
dir() 
setwd("./data")

# list filenames from Komag and vestre jakobselv
filenames <- c("classifications_komagdalen_2016_2020_03_02.csv",
               "classifications_komagdalen_2017_2020_03_02.csv",
               "classifications_komagdalen_2018_2020_03_02.csv",
               "classifications_komagdalen_2019_2020_03_02.csv",
               "classifications_komagdalen_2020_2020_03_02.csv",
               "classifications_vestre_jakobselv_2019_2020_03_02.csv",
               "classifications_vestre_jakosbelv_2020_2020_03_02.csv")

ctdata <- list()
for(i in 1:length(filenames)){
ctdata[[i]] <- read.csv(filenames[i])}

# check that importing was fine
head(ctdata[[1]])
str(ctdata)

# merge list objects to one df
df <- do.call(rbind, ctdata) 

# add columns with site name and date 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
datesplit <- function(filename)  # retrieve date and station from file name
{
  filename2 <- as.character(filename)
  string <- strsplit(filename2,"/")        # split filename at _
  string2 <- strsplit(filename2,"\\\\")[[1]] # split filename at \\\\
  string3 <- strsplit(filename2,"_")
  
  station <-string[[1]][10] # pick the stationID
  date <- string3[[1]][10] # pick the date
  
  # remove last "'" 
  # obtain file name isolated
  newfilename1 <- strsplit(string2[3],"'")[[1]][1]
  
  outdf <- c(station, date, newfilename1) # collect the things you want the function to output
  return(outdf)
} # end 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# run function and store as df
datesite <- as.data.frame(t(sapply(df$fileName,datesplit))) 
names(datesite) <- c("site", "date", "NewFileName")

# add more time/date info
#datesite$date <- as.Date(datesite$date)
#datesite$day <- format(datesite$date, format = "%d")
#datesite$month <- format(datesite$date, format="%m")
#datesite$year <- format(datesite$date, format="%Y")
#datesite$julian <- julian(datesite$date, origin=min(datesite$date))
#datesite$week <- week(datesite$date)
#datesite$yearweek<-format(datesite$date, format="%Y-%W")

# bind date and site df with initial df
df2 <- cbind(df,datesite) 
str(df2) # check that its ok

#0=bad quality, 1=empty, 2=bird, 3=vole, 4=least_weasel, 5=lemming, 6=shrew, 7=stoat, 

# select images classified with more than 90% certainty
df2$answer <- ifelse(df2$confidence1>0.90,df2$guess1,0)

# set answer to NA if quality is bad
table(df2$answer)

df2$answer[df2$answer==0] <- NA
na_vec <- is.na(df2$answer)

summary(df2) # check number of NA's

#remove rows with NA's
df3 <- df2[!na_vec,]

# check that correct number of rows was removed
nrow(df2)-nrow(df3)

# remove unuseful columns
df3 <- df3[,-c(1:2,4:8,10:13)]

# check that the df is fine
str(df3)
table(df3$answer)

# add column for vole, lemming, stoat and least weasel
df3$vole <- ifelse(df3$answer==3,1,0)
table(df3$vole)

df3$lemming <- ifelse(df3$answer==5,1,0)
table(df3$lemming)

df3$stoat <- ifelse(df3$answer==7,1,0)
table(df3$stoat)

df3$least_weasel <- ifelse(df3$answer==4,1,0)
table(df3$least_weasel)

df3$mustela <- df3$stoat+df3$least_weasel

#############################################################
# import metadata (which will provide info on the habitat of each site)
setwd("~/UiT/CameraTrapsForSmallMammals/Data")
metadata <- read.csv("Small mammal camera trap metadata locations.csv", header=T, sep=";")

# check that it looks fine
head(metadata)
names(metadata)[3]<-"site" # change name of site column to be identical to df2

# subset to only keep the site and habitat columns
metadata <- select(metadata, c("site","habitat"))

# add habitat column to the dataset
dat <- left_join(df3, metadata, by="site") # add on info from metadata

#check that it went fine
str(dat)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("~/UiT/GitProjects/MustelidsAndRodents-/data")
# export data frame for metadata
 write.csv(dat, "varanger_camera_answer_rmNA.csv")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# sum animals per day
vole_c <- aggregate(vole~date, data=dat,sum)
vole_c$date <- as.Date(vole_c$date)

stoat_c <- aggregate(stoat~date, data=dat,sum)
stoat_c$date <- as.Date(stoat_c$date)

lemming_c <- aggregate(lemming~date, data=dat,sum)
lemming_c$date <- as.Date(lemming_c$date)

least_weasel_c <- aggregate(least_weasel~date, data=dat,sum)
least_weasel_c$date <- as.Date(least_weasel_c$date)

plot(vole~date, data=vole_c, type="l", col="blue", ylim=c(0,1600))
lines(lemming~date, data=lemming_c, col="red")
lines(stoat~date, data=stoat_c, col="green")
lines(least_weasel~date, data=least_weasel_c, col="orange")
#legend("topleft",legend=c("vole","lemming","stoat","least weasel"),lty=1, lwd=2,
#       col=c("blue","red","green","orange"), cex=0.5)

# daily occpancy
vole_site <- aggregate(vole~date+site, data=dat, max)
vole_occ <- aggregate(vole~date, data=vole_site, mean)

stoat_site <- aggregate(stoat~date+site, data=dat,max)
stoat_occ <- aggregate(stoat~date, data=stoat_site,mean)

lemming_site <- aggregate(lemming~date+site, data=dat,max)
lemming_occ <- aggregate(lemming~date, data=lemming_site,mean)

least_weasel_site <- aggregate(least_weasel~date+site, data=dat,max)
least_weasel_occ <- aggregate(least_weasel~date, data=least_weasel_site,mean)

mustela_site <- aggregate(mustela~date+site, data=dat,max)
mustela_occ <- aggregate(mustela~date, data=mustela_site,mean)

plot(vole_occ$vole , type="l", col="blue")
lines(lemming_occ$lemming, col="red")
lines(stoat_occ$stoat, col="green")
lines(least_weasel_occ$least_weasel, col="orange")
#legend("topleft",legend=c("vole","lemming","stoat","least weasel"),lty=1, lwd=2,
#       col=c("blue","red","green","orange"), cex=0.5)

##############################################
# make occupancy tables for vole and stoat   #
##############################################

# make a column for each site
occ_vole <- spread(vole_site, site, vole)[,-1]
occ_lemming <- spread(lemming_site, site, lemming)[,-1]
occ_stoat <- spread(stoat_site, site, stoat)[,-1]
occ_least_weasel <- spread(least_weasel_site, site, least_weasel)[,-1]
occ_mustela <- spread(mustela_site, site, mustela)[,-1]

# spesifying dimentions of the multi-state occupancy table
ndays <- 7
nsite <- dim(occ_vole)[2]
nprocc <- floor(dim(occ_vole)[1]/ndays)
ntime <- ndays*nprocc

# creating empty array with correct dimentions
occ_vole2 <- array(NA, dim=c(ndays,nprocc,nsite))
occ_lemming2 <- array(NA, dim=c(ndays,nprocc,nsite))
occ_stoat2 <- array(NA, dim=c(ndays,nprocc,nsite))
occ_least_weasel2 <- array(NA, dim=c(ndays,nprocc,nsite))
occ_mustela2 <- array(NA, dim=c(ndays,nprocc,nsite))

# split the observations into primary and secondary occasions (different dimensions in the array) 
for(i in 1:nsite){
  occ_vole2[,,i] <- matrix(unlist(split(occ_vole[1:ntime,i], ceiling(seq_along(occ_vole[1:ntime,i])/ndays))), nrow = ndays, byrow = F)
  occ_lemming2[,,i] <- matrix(unlist(split(occ_lemming[1:ntime,i], ceiling(seq_along(occ_lemming[1:ntime,i])/ndays))), nrow = ndays, byrow = F)
  occ_stoat2[,,i] <- matrix(unlist(split(occ_stoat[1:ntime,i], ceiling(seq_along(occ_stoat[1:ntime,i])/ndays))), nrow = ndays, byrow = F)
  occ_least_weasel2[,,i] <- matrix(unlist(split(occ_least_weasel[1:ntime,i], ceiling(seq_along(occ_least_weasel[1:ntime,i])/ndays))), nrow = ndays, byrow = F)
  occ_mustela2[,,i] <- matrix(unlist(split(occ_mustela[1:ntime,i], ceiling(seq_along(occ_mustela[1:ntime,i])/ndays))), nrow = ndays, byrow = F)
  }

str(occ_vole2)
table(occ_vole2)
table(occ_least_weasel2)
table(occ_stoat2)
table(occ_lemming2)
table(occ_mustela2)

# combine vole and lemmings to a rodent functional group
occ_rodent <- occ_vole2+occ_lemming2 
occ_rodent[occ_rodent==2] <- 1

table(occ_rodent)
dim(occ_rodent)

#############################################################
## Make multi state occupancy tables with block structure ###
#############################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For rodents and mustelids

# set the mustelids to 2
occ_mustela2[occ_mustela2==1] <- 2

# making a multi state array for mustelids and rodents
occ_va <- occ_mustela2+occ_rodent 
table(occ_va)

# changing states from 0,1,2,3 to 1,2,3,4 (needed for the categorical dist in Jags)
occ_va[occ_va==3]<-4
occ_va[occ_va==2]<-3
occ_va[occ_va==1]<-2
occ_va[occ_va==0]<-1

# add block dim
nblocks <- 8

# making empty array to fill in the multi state occupancy matrix
# a problem here is that the sites in Komag have only 11 sites per block, 
# while the sites in VJ have 12 sites per block. Now I create empty sites in komag, which dosen't make sence.
# probably this will also cause problems when we want to include covariates. 

occm_va <- array(NA, dim=c(12,nblocks,nprocc,ndays))

occ_va <- aperm(occ_va, c(3,2,1))

occm_va[1:11,1,,] <- occ_va[1:11,,] 
occm_va[1:11,2,,] <- occ_va[12:22,,] 
occm_va[1:11,3,,] <- occ_va[23:33,,] 
occm_va[1:11,4,,] <- occ_va[34:44,,]
occm_va[,5,,] <- occ_va[45:56,,] 
occm_va[,6,,] <- occ_va[57:68,,] 
occm_va[,7,,] <- occ_va[69:80,,] 
occm_va[,8,,] <- occ_va[81:92,,] 

table(occm_va)
dim(occm_va)

# check number of NA's
summary(occm_va) 

# save 
setwd("C:/Eivind/GitProjects/MustelidsAndRodents-/data")
save(occm_va, file="occm_va_rmNA.rda")


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For voles and mustelids

# set the mustelids to 2
occ_mustela2[occ_mustela2==1] <- 2

# making a multi state array for mustelids and rodents
occ_va_vole <- occ_mustela2+occ_vole2 
table(occ_va_vole)

# changing states from 0,1,2,3 to 1,2,3,4 (needed for the categorical dist in Jags)
occ_va_vole[occ_va_vole==3]<-4
occ_va_vole[occ_va_vole==2]<-3
occ_va_vole[occ_va_vole==1]<-2
occ_va_vole[occ_va_vole==0]<-1

# add block dim
nblocks <- 8

# making empty array to fill in the multi state occupancy matrix
# a problem here is that the sites in Komag have only 11 sites per block, 
# while the sites in VJ have 12 sites per block. Now I create empty sites in komag, which dosen't make sence.
# probably this will also cause problems when we want to include covariates. 

occm_va_vole <- array(NA, dim=c(12,nblocks,nprocc,ndays))

occ_va_vole <- aperm(occ_va_vole, c(3,2,1))

occm_va_vole[1:11,1,,] <- occ_va_vole[1:11,,] 
occm_va_vole[1:11,2,,] <- occ_va_vole[12:22,,] 
occm_va_vole[1:11,3,,] <- occ_va_vole[23:33,,] 
occm_va_vole[1:11,4,,] <- occ_va_vole[34:44,,]
occm_va_vole[,5,,] <- occ_va_vole[45:56,,] 
occm_va_vole[,6,,] <- occ_va_vole[57:68,,] 
occm_va_vole[,7,,] <- occ_va_vole[69:80,,] 
occm_va_vole[,8,,] <- occ_va_vole[81:92,,] 

table(occm_va_vole)
dim(occm_va_vole)

# check number of NA's
summary(occm_va_vole) 

# save 
#setwd("C:/Eivind/GitProjects/MustelidsAndRodents-/data")
setwd("~/UiT/GitProjects/MustelidsAndRodents-/data")
save(occm_va_vole, file="occm_va_vole_rmNA.rda")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# For lemmings and mustelids

# set the mustelids to 2
occ_mustela2[occ_mustela2==1] <- 2

# making a multi state array for mustelids and rodents
occ_va_lem <- occ_mustela2+occ_lemming2 
table(occ_va_lem)

# changing states from 0,1,2,3 to 1,2,3,4 (needed for the categorical dist in Jags)
occ_va_lem[occ_va_lem==3]<-4
occ_va_lem[occ_va_lem==2]<-3
occ_va_lem[occ_va_lem==1]<-2
occ_va_lem[occ_va_lem==0]<-1

# add block dim
nblocks <- 8

# making empty array to fill in the multi state occupancy matrix
# a problem here is that the sites in Komag have only 11 sites per block, 
# while the sites in VJ have 12 sites per block. Now I create empty sites in komag, which dosen't make sence.
# probably this will also cause problems when we want to include covariates. 

occm_va_lem <- array(NA, dim=c(12,nblocks,nprocc,ndays))

occ_va_lem <- aperm(occ_va_lem, c(3,2,1))

occm_va_lem[1:11,1,,] <- occ_va_lem[1:11,,] 
occm_va_lem[1:11,2,,] <- occ_va_lem[12:22,,] 
occm_va_lem[1:11,3,,] <- occ_va_lem[23:33,,] 
occm_va_lem[1:11,4,,] <- occ_va_lem[34:44,,]
occm_va_lem[,5,,] <- occ_va_lem[45:56,,] 
occm_va_lem[,6,,] <- occ_va_lem[57:68,,] 
occm_va_lem[,7,,] <- occ_va_lem[69:80,,] 
occm_va_lem[,8,,] <- occ_va_lem[81:92,,] 

table(occm_va_lem)
dim(occm_va_lem)

# check number of NA's
summary(occm_va_lem) 

# save 
#setwd("C:/Eivind/GitProjects/MustelidsAndRodents-/data")
setwd("~/UiT/GitProjects/MustelidsAndRodents-/data")
save(occm_va_lem, file="occm_var_lem_rmNA.rda")

###########################
# make habitat covariate ##
###########################

str(dat)
# make empty arrays with correct dimentions
hab <- array(NA, dim=c(8, 12))
site_2 <- data.frame()

# extract unique site names
site <- unique(vole_site$site)

# fill in site names in the empty array with block structure
site_2[1, 1:11] <- site[1:11]
site_2[2, 1:11] <- site[12:22]
site_2[3, 1:11] <- site[23:33]
site_2[4, 1:11] <- site[34:44]
site_2[5, 1:12] <- site[45:56]
site_2[6, 1:12] <- site[57:68]
site_2[7, 1:12] <- site[69:80]
site_2[8, 1:12] <- site[81:92]

# make dummy variable describing number of sites in each block
ind <- c(11,11,11,11,12,12,12,12)

dat$habitat[dat$habitat=="hummock mire"] <- 1
dat$habitat[dat$habitat=="snowbed"] <- 2
dat$habitat <- as.integer(dat$habitat)

# extract habitat for each site and place in the empty hab array with block structure
for(b in 1:8){
  for(i in 1:ind[b]){
    hab[b,i] <- dat$habitat[dat$site==site_2[b,i]][1]
  }}

# fill in manually the habitat of that last site in each block in komag (block 1-4)
# these sites was not established before 2020, hence they are not included in the dataset
# however, a covariate from here is still need.

hab[1,12] <- 2
hab[2,12] <- 2
hab[3,12] <- 1
hab[4,12] <- 1

hab

# save habitate covariate
setwd("~/UiT/GitProjects/MustelidsAndRodents-/data")
save(hab, file="hab.rda")

#~ End of script
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~