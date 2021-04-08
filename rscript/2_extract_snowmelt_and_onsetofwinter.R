# load libraries
library(dplyr)

# set wd
wd1 <- "C:/Eivind/GitProjects/MustelidsAndRodents-/data"
setwd(wd1)

# read in temperature data
dat <- read.csv("varanger_cameradata_final.csv")

str(dat)

dat$date <- as.Date(dat$datetime)

sitenames <- sort(unique(dat$site))

par(mfrow=c(4,4))

for(i in 51:66){
  site1 <- filter(dat, site==sitenames[i])
  plot(site1$temp ~ site1$date, main=sitenames[i])
  abline(h=0)
}

site1 <- filter(dat, site==sitenames[45])
site2 <- filter(dat, site=="ko_ga_hu_2")
site2 <- filter(dat, site=="ko_ga_hu_2")


plot(site1$temp ~ site1$date)
abline(h=0)

plot(site2$temp ~ site2$date)
abline(h=0)

# an idea is to subset every unique year and site combination. 
# then snowmelt can be extracted as the first day with above for example 2 degrees
# for onset of snow it's a bit more tricky. But remember that you need these
# covariates on a weekly time scale. Hence, you could calculate weekly sd's
# then presence of snow can be described as the weeks where the weekly sd
# are below some threshold.
