# load libraries
library(dplyr)

# set wd
wd1 <- "C:/Eivind/GitProjects/MustelidsAndRodents-/data"
setwd(wd1)

# read in image annotation data
data1 <- read.csv("varanger_camera_answer.csv", head=TRUE)
str(data1)

# set metadata directory
metwd <- ("metadata")
#dirloc <-dir("metadata")
dirmet <- dir(metwd)

ctlist <- list()
listindex=0 # to track how long the list of data.frames will be
str(ctlist)

# combine meta data files

#for(k in 1: length(dirloc)){
# newdir2 <- dirloc[k]

for (i in 1:length(dirmet)){
  newdir <- dirmet[i]
  newmet <- paste0(metwd,"/",newdir)
  
  yeardir <- dir(newmet)
  
  for(j in 1:length(yeardir))
  { 
    newmet2 <- paste0(metwd,"/",newdir,"/",yeardir[j])
    sitedir <- dir(newmet2) 
    
    for(s in 1:length(sitedir)){
    listindex = listindex + 1
    ct_data <- read.table(paste0(newmet2,"/",sitedir[s]), header=TRUE)
    ctlist[[listindex]] <- ct_data
    }
  }
}

# convert list into data frame
metadf <- do.call(rbind, ctlist)

# remove columns that you don't need
metadf1 <- select(metadf, DateTimeOriginal, TriggerMode, Sequence, AmbientTemperature, NewFileName)
names(metadf1) <- c("datetime", "trigger_mode", "sequence", "temp", "NewFileName")

metadf1$datetime <- as.character(as.POSIXct(metadf1$datetime, format="%Y:%m:%d %H:%M:%S"))

str(metadf1)

# load meta data from Komag 2016-2018 which have a different format
setwd("./metadata2")
dat2 <- list()

for(i in 1:length(dir())){
dat2[[i]]<- read.table(dir()[i], header=T, sep=";")
}

metadf2 <- do.call(rbind, dat2)
metadf2$datetime <- as.character(as.POSIXct(paste(metadf2$t_date, metadf2$t_time), format="%Y-%m-%d %H:%M:%S"))
str(metadf2)

metadf3 <- select(metadf2, datetime, v_trigger_mode, v_sequence, v_temperature, v_image_name)
str(metadf3)
names(metadf3) <- c("datetime", "trigger_mode", "sequence", "temp", "NewFileName")

# combine the two different meta data files
metadf4 <- rbind(metadf1, metadf3)

# add metadata variables to haakoya camera data
cameratrap1 <- dplyr::left_join(data1, metadf4, by="NewFileName")

# format time and data
cameratrap1$datetime <- strptime(cameratrap1$datetime, format="%Y-%m-%d %H:%M:%S", tz="CET") #Central European Time

# remove unimportant variables
cameratrap2 <- select(cameratrap1, site, habitat, datetime, temp, answer, vole, lemming, stoat, least_weasel, mustela, confidence1)
tail(cameratrap2)

#plot(cameratrap1$temp)

setwd("../")
write.csv(cameratrap2, "varanger_cameradata_final.csv")

#~ End of Script