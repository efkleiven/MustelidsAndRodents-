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

# add julian date
#dat$julian <- julian(dat$date, origin = min(dat$date, na.rm=T))
#summary(dat$julian)

#filter out only timelapse images

unique(dat$trigger_mode)

dat2 <- filter(dat, trigger_mode == "T" |  trigger_mode == "timelapse")
str(dat2)

# check that min date is the same as in occupancy tables (2015-08-21)
min(dat2$datetime, na.rm = T)
max(dat2$datetime, na.rm = T)

# calculate mean daily temp
temp_day <- aggregate(temp~date+site, data=dat2, mean)
str(temp_day)

# add year column
temp_day$year<-format(temp_day$date, format="%Y")

###
sitenames[1]
years <- unique(temp_day$year)

df <- filter(temp_day, site==sitenames[1] & year==years[2])

plot(temp~date, data=df)

# you could add season(spring/autumn as a dummy variable) then filter only springs and find first last day with below 0,
# and do the same similarly for autumn, where you choose the first day with below 0.   

temp_day$seas <- NA

str(temp_day)

temp_day$seas <- ifelse(temp_day$date < as.Date("2016-03-01"), temp_day$seas <- "freez1",
                        ifelse(temp_day$date < as.Date("2016-09-01"), temp_day$seas <- "melt1",
                          ifelse(temp_day$date < as.Date("2017-03-01"), temp_day$seas <- "freez2", 
                            ifelse(temp_day$date < as.Date("2017-09-01"), temp_day$seas <-"melt2", 
                              ifelse(temp_day$date < as.Date("2018-03-01"), temp_day$seas <-"freez3",
                                ifelse(temp_day$date < as.Date("2018-09-01"), temp_day$seas <-"melt3",
                                  ifelse(temp_day$date < as.Date("2019-03-01"), temp_day$seas <-"freez4",
                                    ifelse(temp_day$date < as.Date("2019-09-01"), temp_day$seas <-"melt4",
                                      ifelse(temp_day$date < as.Date("2020-03-01"), temp_day$seas <-"freez5", temp_day$seas <-"melt5"
                                             )))))))))

head(temp_day)
tail(temp_day)
table(temp_day$seas)



#############

# find first day with below 0 temp in autumn

sitenames <- unique(temp_day$site)
years <- unique(temp_day$year)
dat15 <- c()
dat16 <- c()
dat17 <- c()
dat18 <- c()
dat19 <- c()
dat20 <- c()


for(i in 1:length(sitenames)){
dat <- filter(temp_day, year==2015 & seas=="freez1", site==sitenames[i])
dat15[i] <- as.character(min(dat$date[dat$temp<0]))
}

for(i in 1:length(sitenames)){
  dat <- filter(temp_day, year==2016 & seas=="freez2", site==sitenames[i])
  dat16[i] <- as.character(min(dat$date[dat$temp<0]))
}

for(i in 1:length(sitenames)){
  dat <- filter(temp_day, year==2017 & seas=="freez3", site==sitenames[i])
  dat17[i] <- as.character(min(dat$date[dat$temp<0]))
}

for(i in 1:length(sitenames)){
  dat <- filter(temp_day, year==2018 & seas=="freez4", site==sitenames[i])
  dat18[i] <- as.character(min(dat$date[dat$temp<0]))
}

for(i in 1:length(sitenames)){
  dat <- filter(temp_day, year==2019 & seas=="freez5", site==sitenames[i])
  dat19[i] <- as.character(min(dat$date[dat$temp<0]))
}



datfreeze <- data.frame(site=sitenames, s2015=dat15, s2016 =dat16, s2017=dat17, s2018=dat18, s2019=dat19 )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Fill in block and habitat mean for NA's

# add columns with block id and habitat
for(i in 1:length(sitenames)){
  split1 <- strsplit(datfreeze$site[i], split = "_")
  datfreeze$block[i] <- paste(split1[[1]][1],split1[[1]][2], sep="_")
  datfreeze$hab[i] <- split1[[1]][3]
}

#fill NA's with block and habitat mean
blocknames <- unique(datfreeze$block)
hab <- unique(datfreeze$hab)

datfreeze$s2015 <- as.Date(datfreeze$s2015)
datfreeze$s2016 <- as.Date(datfreeze$s2016)
datfreeze$s2017 <- as.Date(datfreeze$s2017)
datfreeze$s2018 <- as.Date(datfreeze$s2018)
datfreeze$s2019 <- as.Date(datfreeze$s2019)


# block 1
dat11 <- filter(datfreeze, block==blocknames[1] & hab==hab[1])
dat12 <- filter(datfreeze, block==blocknames[1] & hab==hab[2])

dat11$s2015[is.na(dat11$s2015)] <- mean(dat11$s2015, na.rm=T)
dat12$s2015[is.na(dat12$s2015)] <- mean(dat12$s2015, na.rm=T)

dat11$s2016[is.na(dat11$s2016)] <- mean(dat11$s2016, na.rm=T)
dat12$s2016[is.na(dat12$s2016)] <- mean(dat12$s2016, na.rm=T)

dat11$s2017[is.na(dat11$s2017)] <- mean(dat11$s2017, na.rm=T)
dat12$s2017[is.na(dat12$s2017)] <- mean(dat12$s2017, na.rm=T)

dat11$s2018[is.na(dat11$s2018)] <- mean(dat11$s2018, na.rm=T)
dat12$s2018[is.na(dat12$s2018)] <- mean(dat12$s2018, na.rm=T)

dat11$s2019[is.na(dat11$s2019)] <- mean(dat11$s2019, na.rm=T)
dat12$s2019[is.na(dat12$s2019)] <- mean(dat12$s2019, na.rm=T)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# block 2
dat21 <- filter(datfreeze, block==blocknames[2] & hab==hab[1])
dat22 <- filter(datfreeze, block==blocknames[2] & hab==hab[2])

dat21$s2015[is.na(dat21$s2015)] <- mean(dat21$s2015, na.rm=T)
dat22$s2015[is.na(dat22$s2015)] <- mean(dat22$s2015, na.rm=T)

dat21$s2016[is.na(dat21$s2016)] <- mean(dat21$s2016, na.rm=T)
dat22$s2016[is.na(dat22$s2016)] <- mean(dat22$s2016, na.rm=T)

dat21$s2017[is.na(dat21$s2017)] <- mean(dat21$s2017, na.rm=T)
dat22$s2017[is.na(dat22$s2017)] <- mean(dat22$s2017, na.rm=T)

dat21$s2018[is.na(dat21$s2018)] <- mean(dat21$s2018, na.rm=T)
dat22$s2018[is.na(dat22$s2018)] <- mean(dat22$s2018, na.rm=T)

dat21$s2019[is.na(dat21$s2019)] <- mean(dat21$s2019, na.rm=T)
dat22$s2019[is.na(dat22$s2019)] <- mean(dat22$s2019, na.rm=T)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# block 3
dat31 <- filter(datfreeze, block==blocknames[3] & hab==hab[1])
dat32 <- filter(datfreeze, block==blocknames[3] & hab==hab[2])

dat31$s2015[is.na(dat31$s2015)] <- mean(dat31$s2015, na.rm=T)
dat32$s2015[is.na(dat32$s2015)] <- mean(dat32$s2015, na.rm=T)

dat31$s2016[is.na(dat31$s2016)] <- mean(dat31$s2016, na.rm=T)
dat32$s2016[is.na(dat32$s2016)] <- mean(dat32$s2016, na.rm=T)

dat31$s2017[is.na(dat31$s2017)] <- mean(dat31$s2017, na.rm=T)
dat32$s2017[is.na(dat32$s2017)] <- mean(dat32$s2017, na.rm=T)

dat31$s2018[is.na(dat31$s2018)] <- mean(dat31$s2018, na.rm=T)
dat32$s2018[is.na(dat32$s2018)] <- mean(dat32$s2018, na.rm=T)

dat31$s2019[is.na(dat31$s2019)] <- mean(dat31$s2019, na.rm=T)
dat32$s2019[is.na(dat32$s2019)] <- mean(dat32$s2019, na.rm=T)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# block 4
dat41 <- filter(datfreeze, block==blocknames[4] & hab==hab[1])
dat42 <- filter(datfreeze, block==blocknames[4] & hab==hab[2])

dat41$s2015[is.na(dat41$s2015)] <- mean(dat41$s2015, na.rm=T)
dat42$s2015[is.na(dat42$s2015)] <- mean(dat42$s2015, na.rm=T)

dat41$s2016[is.na(dat41$s2016)] <- mean(dat41$s2016, na.rm=T)
dat42$s2016[is.na(dat42$s2016)] <- mean(dat42$s2016, na.rm=T)

dat41$s2017[is.na(dat41$s2017)] <- mean(dat41$s2017, na.rm=T)
dat42$s2017[is.na(dat42$s2017)] <- mean(dat42$s2017, na.rm=T)

dat41$s2018[is.na(dat41$s2018)] <- mean(dat41$s2018, na.rm=T)
dat42$s2018[is.na(dat42$s2018)] <- mean(dat42$s2018, na.rm=T)

dat41$s2019[is.na(dat41$s2019)] <- mean(dat41$s2019, na.rm=T)
dat42$s2019[is.na(dat42$s2019)] <- mean(dat42$s2019, na.rm=T)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# block 5, this are VJ blocks which were not deployed until summer 2018.
# hence NA's in 15-17 will be dealt with later

dat51 <- filter(datfreeze, block==blocknames[5] & hab==hab[1])
dat52 <- filter(datfreeze, block==blocknames[5] & hab==hab[2])

dat51$s2018[is.na(dat41$s2018)] <- mean(dat51$s2018, na.rm=T)
dat52$s2018[is.na(dat42$s2018)] <- mean(dat52$s2018, na.rm=T)

dat51$s2019[is.na(dat41$s2019)] <- mean(dat51$s2019, na.rm=T)
dat52$s2019[is.na(dat42$s2019)] <- mean(dat52$s2019, na.rm=T)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# block 6, this are VJ blocks which were not deployed until summer 2018.
# hence NA's in 15-17 will be dealt with later

dat61 <- filter(datfreeze, block==blocknames[6] & hab==hab[1])
dat62 <- filter(datfreeze, block==blocknames[6] & hab==hab[2])

dat61$s2018[is.na(dat61$s2018)] <- mean(dat61$s2018, na.rm=T)
dat62$s2018[is.na(dat62$s2018)] <- mean(dat62$s2018, na.rm=T)

dat61$s2019[is.na(dat61$s2019)] <- mean(dat61$s2019, na.rm=T)
dat62$s2019[is.na(dat62$s2019)] <- mean(dat62$s2019, na.rm=T)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# block 7, this are VJ blocks which were not deployed until summer 2018.
# hence NA's in 15-17 will be dealt with later

dat71 <- filter(datfreeze, block==blocknames[7] & hab==hab[1])
dat72 <- filter(datfreeze, block==blocknames[7] & hab==hab[2])

dat71$s2018[is.na(dat71$s2018)] <- mean(dat71$s2018, na.rm=T)
dat72$s2018[is.na(dat72$s2018)] <- mean(dat72$s2018, na.rm=T)

dat71$s2019[is.na(dat71$s2019)] <- mean(dat71$s2019, na.rm=T)
dat72$s2019[is.na(dat72$s2019)] <- mean(dat72$s2019, na.rm=T)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# block 8, this are VJ blocks which were not deployed until summer 2018.
# hence NA's in 15-17 will be dealt with later

dat81 <- filter(datfreeze, block==blocknames[8] & hab==hab[1])
dat82 <- filter(datfreeze, block==blocknames[8] & hab==hab[2])

dat81$s2018[is.na(dat81$s2018)] <- mean(dat81$s2018, na.rm=T)
dat82$s2018[is.na(dat82$s2018)] <- mean(dat82$s2018, na.rm=T)

dat81$s2019[is.na(dat81$s2019)] <- mean(dat81$s2019, na.rm=T)
dat82$s2019[is.na(dat82$s2019)] <- mean(dat82$s2019, na.rm=T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
freeze <- rbind(dat11,dat12,dat21,dat22,dat31,dat32,dat41,dat42,dat51,dat52,dat61,dat62,dat71,dat72,dat81,dat82 )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# find last day with below 0 temp in spring

sitenames <- unique(temp_day$site)
years <- unique(temp_day$year)

dat16 <- c()
dat17 <- c()
dat18 <- c()
dat19 <- c()
dat20 <- c()


for(i in 1:length(sitenames)){
  dat <- filter(temp_day, year==2016 & seas=="melt1", site==sitenames[i])
  dat16[i] <- as.character(max(dat$date[dat$temp<0]))
}

for(i in 1:length(sitenames)){
  dat <- filter(temp_day, year==2017 & seas=="melt2", site==sitenames[i])
  dat17[i] <- as.character(max(dat$date[dat$temp<0]))
}

for(i in 1:length(sitenames)){
  dat <- filter(temp_day, year==2018 & seas=="melt3", site==sitenames[i])
  dat18[i] <- as.character(max(dat$date[dat$temp<0]))
}

for(i in 1:length(sitenames)){
  dat <- filter(temp_day, year==2019 & seas=="melt4", site==sitenames[i])
  dat19[i] <- as.character(max(dat$date[dat$temp<0]))
}

for(i in 1:length(sitenames)){
  dat <- filter(temp_day, year==2020 & seas=="melt5", site==sitenames[i])
  dat20[i] <- as.character(max(dat$date[dat$temp<0]))
}


datmelt <- data.frame(site=sitenames, s2016 =dat16, s2017=dat17, s2018=dat18, s2019=dat19, s2020=dat20 )

# filter block by block
# add column with block id
for(i in 1:length(sitenames)){
  split1 <- strsplit(datmelt$site[i], split = "_")
  datmelt$block[i] <- paste(split1[[1]][1],split1[[1]][2], sep="_")
  datmelt$hab[i] <- split1[[1]][3]
}


#fill NA's with block and habitat mean
blocknames <- unique(datmelt$block)
hab <- unique(datmelt$hab)

datmelt$s2016 <- as.Date(datmelt$s2016)
datmelt$s2017 <- as.Date(datmelt$s2017)
datmelt$s2018 <- as.Date(datmelt$s2018)
datmelt$s2019 <- as.Date(datmelt$s2019)
datmelt$s2020 <- as.Date(datmelt$s2020)

# block 1
dat11 <- filter(datmelt, block==blocknames[1] & hab==hab[1])
dat12 <- filter(datmelt, block==blocknames[1] & hab==hab[2])

dat11$s2016[is.na(dat11$s2016)] <- mean(dat11$s2016, na.rm=T)
dat12$s2016[is.na(dat12$s2016)] <- mean(dat12$s2016, na.rm=T)

dat11$s2017[is.na(dat11$s2017)] <- mean(dat11$s2017, na.rm=T)
dat12$s2017[is.na(dat12$s2017)] <- mean(dat12$s2017, na.rm=T)

dat11$s2018[is.na(dat11$s2018)] <- mean(dat11$s2018, na.rm=T)
dat12$s2018[is.na(dat12$s2018)] <- mean(dat12$s2018, na.rm=T)

dat11$s2019[is.na(dat11$s2019)] <- mean(dat11$s2019, na.rm=T)
dat12$s2019[is.na(dat12$s2019)] <- mean(dat12$s2019, na.rm=T)

dat11$s2020[is.na(dat11$s2020)] <- mean(dat11$s2020, na.rm=T)
dat12$s2020[is.na(dat12$s2020)] <- mean(dat12$s2020, na.rm=T)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# block 2
dat21 <- filter(datmelt, block==blocknames[2] & hab==hab[1])
dat22 <- filter(datmelt, block==blocknames[2] & hab==hab[2])

dat21$s2016[is.na(dat21$s2016)] <- mean(dat21$s2016, na.rm=T)
dat22$s2016[is.na(dat22$s2016)] <- mean(dat22$s2016, na.rm=T)

dat21$s2017[is.na(dat21$s2017)] <- mean(dat21$s2017, na.rm=T)
dat22$s2017[is.na(dat22$s2017)] <- mean(dat22$s2017, na.rm=T)

dat21$s2018[is.na(dat21$s2018)] <- mean(dat21$s2018, na.rm=T)
dat22$s2018[is.na(dat22$s2018)] <- mean(dat22$s2018, na.rm=T)

dat21$s2019[is.na(dat21$s2019)] <- mean(dat21$s2019, na.rm=T)
dat22$s2019[is.na(dat22$s2019)] <- mean(dat22$s2019, na.rm=T)

dat21$s2020[is.na(dat21$s2020)] <- mean(dat21$s2020, na.rm=T)
dat22$s2020[is.na(dat22$s2020)] <- mean(dat22$s2020, na.rm=T)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# block 3
dat31 <- filter(datmelt, block==blocknames[3] & hab==hab[1])
dat32 <- filter(datmelt, block==blocknames[3] & hab==hab[2])

dat31$s2016[is.na(dat31$s2016)] <- mean(dat31$s2016, na.rm=T)
dat32$s2016[is.na(dat32$s2016)] <- mean(dat32$s2016, na.rm=T)

dat31$s2017[is.na(dat31$s2017)] <- mean(dat31$s2017, na.rm=T)
dat32$s2017[is.na(dat32$s2017)] <- mean(dat32$s2017, na.rm=T)

dat31$s2018[is.na(dat31$s2018)] <- mean(dat31$s2018, na.rm=T)
dat32$s2018[is.na(dat32$s2018)] <- mean(dat32$s2018, na.rm=T)

dat31$s2019[is.na(dat31$s2019)] <- mean(dat31$s2019, na.rm=T)
dat32$s2019[is.na(dat32$s2019)] <- mean(dat32$s2019, na.rm=T)

dat31$s2020[is.na(dat31$s2020)] <- mean(dat31$s2020, na.rm=T)
dat32$s2020[is.na(dat32$s2020)] <- mean(dat32$s2020, na.rm=T)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# block 4
dat41 <- filter(datmelt, block==blocknames[4] & hab==hab[1])
dat42 <- filter(datmelt, block==blocknames[4] & hab==hab[2])

dat41$s2016[is.na(dat41$s2016)] <- mean(dat41$s2016, na.rm=T)
dat42$s2016[is.na(dat42$s2016)] <- mean(dat42$s2016, na.rm=T)

dat41$s2017[is.na(dat41$s2017)] <- mean(dat41$s2017, na.rm=T)
dat42$s2017[is.na(dat42$s2017)] <- mean(dat42$s2017, na.rm=T)

dat41$s2018[is.na(dat41$s2018)] <- mean(dat41$s2018, na.rm=T)
dat42$s2018[is.na(dat42$s2018)] <- mean(dat42$s2018, na.rm=T)

dat41$s2019[is.na(dat41$s2019)] <- mean(dat41$s2019, na.rm=T)
dat42$s2019[is.na(dat42$s2019)] <- mean(dat42$s2019, na.rm=T)

dat41$s2020[is.na(dat41$s2020)] <- mean(dat41$s2020, na.rm=T)
dat42$s2020[is.na(dat42$s2020)] <- mean(dat42$s2020, na.rm=T)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# block 5, this is VJ, where cameras where not deployed until summer 2018,
# hence all melt dates in 2016-2018 will be NA's and these will be dealt with later

dat51 <- filter(datmelt, block==blocknames[5] & hab==hab[1])
dat52 <- filter(datmelt, block==blocknames[5] & hab==hab[2])

dat51$s2019[is.na(dat51$s2019)] <- mean(dat51$s2019, na.rm=T)
dat52$s2019[is.na(dat52$s2019)] <- mean(dat52$s2019, na.rm=T)

dat51$s2020[is.na(dat51$s2020)] <- mean(dat51$s2020, na.rm=T)
dat52$s2020[is.na(dat52$s2020)] <- mean(dat52$s2020, na.rm=T)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# block 6
dat61 <- filter(datmelt, block==blocknames[6] & hab==hab[1])
dat62 <- filter(datmelt, block==blocknames[6] & hab==hab[2])

dat61$s2019[is.na(dat61$s2019)] <- mean(dat61$s2019, na.rm=T)
dat62$s2019[is.na(dat62$s2019)] <- mean(dat62$s2019, na.rm=T)

dat61$s2020[is.na(dat61$s2020)] <- mean(dat61$s2020, na.rm=T)
dat62$s2020[is.na(dat62$s2020)] <- mean(dat62$s2020, na.rm=T)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# block 7
dat71 <- filter(datmelt, block==blocknames[5] & hab==hab[1])
dat72 <- filter(datmelt, block==blocknames[5] & hab==hab[2])

dat71$s2019[is.na(dat71$s2019)] <- mean(dat71$s2019, na.rm=T)
dat72$s2019[is.na(dat72$s2019)] <- mean(dat72$s2019, na.rm=T)

dat71$s2020[is.na(dat71$s2020)] <- mean(dat71$s2020, na.rm=T)
dat72$s2020[is.na(dat72$s2020)] <- mean(dat72$s2020, na.rm=T)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# block 8
dat81 <- filter(datmelt, block==blocknames[8] & hab==hab[1])
dat82 <- filter(datmelt, block==blocknames[8] & hab==hab[2])

dat81$s2019[is.na(dat81$s2019)] <- mean(dat81$s2019, na.rm=T)
dat82$s2019[is.na(dat82$s2019)] <- mean(dat82$s2019, na.rm=T)

dat81$s2020[is.na(dat81$s2020)] <- mean(dat81$s2020, na.rm=T)
dat82$s2020[is.na(dat82$s2020)] <- mean(dat82$s2020, na.rm=T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
melt <- rbind(dat11,dat12,dat21,dat22,dat31,dat32,dat41,dat42,dat51,dat52,dat61,dat62,dat71,dat72,dat81,dat82 )
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# use melt and freeze to make covariate with same temporal and spatial structure as multi-state occupancy table

# make time vector with dates from onset of monitoring to the end


dates <- seq(as.Date("2015-08-21"), as.Date("2020-09-05"), by="days")
season <- array(NA, dim=c(length(dates), length(sitenames)))

freeze[1,]
melt[1,]

# for sites in komag
for(j in 1:44){
  for(i in 1:length(dates)){
    season[i,j] <- ifelse(dates[i] < freeze$s2015[j], "summer", 
                    ifelse(dates[i] < melt$s2016[j], "winter",
                      ifelse(dates[i] < freeze$s2016[j], "summer",
                        ifelse(dates[i] < melt$s2017[j], "winter",
                          ifelse(dates[i] < freeze$s2017[j], "summer",
                            ifelse(dates[i] < melt$s2018[j], "winter",
                              ifelse(dates[i] < freeze$s2018[j], "summer",
                                ifelse(dates[i] < melt$s2019[j], "winter",
                                  ifelse(dates[i] < freeze$s2019[j], "summer",
                                    ifelse(dates[i] < melt$s2020[j], "winter", "summer"))))))))))
}}

# for sites in VJ
for(j in 45:92){
  for(i in 1:length(dates)){
    season[i,j] <- ifelse(dates[i] < as.Date("2018-07-01"), NA,
                     ifelse(dates[i] < freeze$s2018[j], "summer",
                       ifelse(dates[i] < melt$s2019[j], "winter",
                         ifelse(dates[i] < freeze$s2019[j], "summer",
                           ifelse(dates[i] < melt$s2020[j], "winter", "summer")))))
  }}


season[season=="winter"] <- 1
season[season=="summer"] <- 2

# make empty array with correct structure for seas covariate
nblocks <- 8
ndays <- 7
nprocc <- 263
totdays <- nprocc*ndays

season_cov <- array(NA, dim=c(12,nblocks,nprocc))


# calculate weekly means
n = 7

length(as.numeric(season[seq(1, nrow(season), n), 1]))

#block 1
for(i in 1:11){
season_cov[i,1,] <- tapply(as.numeric(season[,i]), rep(seq_along(season[,i]), each = n, length.out = length(season[,i])), median)[-264]
}

#block 2
for(i in 1:11){
  season_cov[i,2,] <- tapply(as.numeric(season[,i+11]), rep(seq_along(season[,i+11]), each = n, length.out = length(season[,i+11])), median)[-264]
}

#block 3
for(i in 1:11){
  season_cov[i,3,] <- tapply(as.numeric(season[,i+22]), rep(seq_along(season[,i+22]), each = n, length.out = length(season[,i+22])), median)[-264]
}

#block 4
for(i in 1:11){
  season_cov[i,4,] <- tapply(as.numeric(season[,i+33]), rep(seq_along(season[,i+33]), each = n, length.out = length(season[,i+33])), median)[-264]
}

#block 5
for(i in 1:12){
  season_cov[i,5,] <- tapply(as.numeric(season[,i+44]), rep(seq_along(season[,i+44]), each = n, length.out = length(season[,i+44])), FUN=median)[-264]
}

#block 6
for(i in 1:12){
  season_cov[i,6,] <- tapply(as.numeric(season[,i+56]), rep(seq_along(season[,i+56]), each = n, length.out = length(season[,i+56])), median)[-264]
}

#block 7
for(i in 1:12){
  season_cov[i,7,] <- tapply(as.numeric(season[,i+68]), rep(seq_along(season[,i+68]), each = n, length.out = length(season[,i+68])), median)[-264]
}

#block 8
for(i in 1:12){
  season_cov[i,8,] <- tapply(as.numeric(season[,i+80]), rep(seq_along(season[,i+80]), each = n, length.out = length(season[,i+80])), median)[-264]
}

summary(season_cov)
dim(season_cov)

head(season_cov)
season_cov[11,5,]

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save season_cov
# save 
save(season_cov, file="season_cov.rda")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# make block level seasonal cov
# take median block value

season_cov_block <- apply(season_cov, c(2,3), median, na.rm=T )

table(season_cov_block) # 1.5 comes when 6 sites are 1 and 6 are 2, this only happens 10 times therefore we can simply set these to 1

season_cov_block[season_cov_block==1.5] <- 1

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# save season_cov_block
save(season_cov_block, file="season_cov_block.rda")
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~ The end!
