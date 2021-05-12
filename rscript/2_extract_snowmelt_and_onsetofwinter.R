# load libraries
library(dplyr)
library(tidyr)

# set wd
wd1 <- "C:/Eivind/GitProjects/MustelidsAndRodents-/data"
setwd(wd1)

# read in temperature data
dat <- read.csv("varanger_cameradata_final.csv")

str(dat)

dat$date <- as.Date(dat$datetime)

# plot temperatures at some sites
sitenames <- sort(unique(dat$site))

par(mfrow=c(4,3))

for(i in 1:11){
  site1 <- filter(dat, site==sitenames[i])
  plot(site1$temp ~ site1$date, main=sitenames[i])
  abline(h=0, col="red")
}

par(mfrow=c(4,3))

for(i in 12:22){
  site1 <- filter(dat, site==sitenames[i])
  plot(site1$temp ~ site1$date, main=sitenames[i])
  abline(h=0, col="red")
}

par(mfrow=c(4,3))

for(i in 23:33){
  site1 <- filter(dat, site==sitenames[i])
  plot(site1$temp ~ site1$date, main=sitenames[i])
  abline(h=0, col="red")
}


par(mfrow=c(4,3))

for(i in 34:44){
  site1 <- filter(dat, site==sitenames[i])
  plot(site1$temp ~ site1$date, main=sitenames[i])
  abline(h=0, col="red")
}


par(mfrow=c(4,3))

for(i in 45:56){
  site1 <- filter(dat, site==sitenames[i])
  plot(site1$temp ~ site1$date, main=sitenames[i])
  abline(h=0, col="red")
}

for(i in 57:68){
  site1 <- filter(dat, site==sitenames[i])
  plot(site1$temp ~ site1$date, main=sitenames[i])
  abline(h=0, col="red")
}

for(i in 69:80){
  site1 <- filter(dat, site==sitenames[i])
  plot(site1$temp ~ site1$date, main=sitenames[i])
  abline(h=0, col="red")
}

for(i in 81:92){
  site1 <- filter(dat, site==sitenames[i])
  plot(site1$temp ~ site1$date, main=sitenames[i])
  abline(h=0, col="red")
}

# add julian date
#dat$julian <- julian(dat$date, origin = min(dat$date, na.rm=T))
#summary(dat$julian)

#filter out only timelapse images

unique(dat$trigger_mode)

dat2 <- filter(dat, trigger_mode == "T" |  trigger_mode == "timelapse")

str(dat2)
# check that min date is the same as in occupancy tables (2015-08-21)
min(dat2$datetime, na.rm = T)

# calculate mean daily temp
temp_day <- aggregate(temp~date+site, data=dat2, mean)

# spread into df with site as columns
temp_site <- spread(temp_day, site, temp)[, -1]

# adding empty columns for the sites without observations
temp_site$ko_ga_sn_6 <- NA
temp_site$ko_hu_sn_6 <- NA
temp_site$ko_kj_hu_6 <- NA
temp_site$ko_ry_hu_6 <- NA

# add week variable
temp_site2 <- temp_site[1:1841,]
temp_site2$week <- rep(1:(1841/7), each = 7)

# make site name vector
sitenames <- sort(c(unique(dat2$site), "ko_ga_sn_6", "ko_hu_sn_6", "ko_kj_hu_6","ko_ry_hu_6"))

# make empty list to store output when looping over sites
temp_site3 <- list()

# aggregate sd over weeks
for(i in 1:length(sitenames)){
  
  dat <- select (temp_site2, c(week , sitenames[i]))
  names(dat) <- c("week","x")

  temp_site3[[i]] <- dat %>%
                      group_by(week) %>%
                      summarize_all(sd, na.rm = TRUE)
}

# arrange in array
temp4 <- array(NA, dim=c(96, 263))

for(i in 1:96){
temp4[i,] <- temp_site3[[i]][,2]$x}

temp4[temp4<1] <- 0   # snow
temp4[temp4>1] <- 1   # no snow

hist(temp4)

par(mfrow=c(4,3))

for(i in 1:11){
plot(temp4[i,], main = sitenames[i])}

for(i in 13:24){
  plot(temp4[i,], main = sitenames[i])}

for(i in c(25:29,31:36)){
  plot(temp4[i,], main = sitenames[i])}

for(i in c(37:41,43:48)){
  plot(temp4[i,], main = sitenames[i])}

for(i in 49:60){
  plot(temp4[i,], main = sitenames[i])}

for(i in 61:72){
  plot(temp4[i,], main = sitenames[i])}

for(i in 73:84){
  plot(temp4[i,], main = sitenames[i])}

for(i in 85:96){
  plot(temp4[i,], main = sitenames[i])}
# over all it looks like there is a seasonal rythm here, so this might be ok

# fill in NA's with mean for the week
temp5 <- temp4

# loop to fill cells with mean for all sites if the site measurment is NA
for(j in 1:96){
  for(i in 1:263){
    if(is.na(temp4[j,i])==T) {
      temp5[j,i] <- round(mean(temp4[,i], na.rm=T)) 
    }
  }}

dim(temp5)

# save 
save(temp5, file="temp.rda")

#~ The End
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
