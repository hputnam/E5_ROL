#Title: Chlorophyll Data
#Project: E5
#Author: HM Putnam 
#Edited by: HM Putnam
#Date Last Modified: 20191031
#See Readme file for details

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
#if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 

Data <- read.csv("RAnalysis/Data/20191031_Plate1_CHL.csv")

# For Dinos from Jeffrey and Humphrey 1975 in 100% acetone
# chla = 11.43*A663 - 0.64*A630
# chlc2 = 27.09*A630 - 3.63*A663

#remove background at 750nm

Data$X630 <- Data$X630-Data$X750
Data$X663 <- Data$X663-Data$X750

#Calculate Chl-a and Chl-c2 in µg/ml
Data$chla.ug.ml <- 11.43*Data$X663 - 0.64*Data$X630
Data$chlc2.ug.ml <- 27.09*Data$X630 - 3.63*Data$X663

Sample.Info <- read.csv("RAnalysis/Data/20191031_Plate1_CHL_WELLS.csv")
Meta <- read.csv("RAnalysis/Data/Sample_Metadata.csv")

Data <- merge(Data, Sample.Info, by="Well")
Data <- merge(Data, Meta, by="Tube.ID")

Data.avg.a <- aggregate(chla.ug.ml ~ Sample.ID, data=Data, FUN=mean)
Data.avg.c2 <- aggregate(chlc2.ug.ml ~ Sample.ID, data=Data, FUN=mean)
Data.avg <- cbind(Data.avg.a, Data.avg.c2$chlc2.ug.ml)
Data.avg <- merge(Data.avg, Meta, by="Sample.ID")
colnames(Data.avg) <- c("Sample.ID","chla.ug.ml", "chlc2.ug.ml", "Tube.ID", "Homog.Vol.ml", "Species", "Site")            
Data.avg$chla.ug <- Data.avg$chla.ug.ml * Data.avg$Homog.Vol.ml
Data.avg$chlc2.ug <- Data.avg$chlc2.ug.ml * Data.avg$Homog.Vol.ml

Surface.Area <- read.csv("RAnalysis/Data/E5surface_area.csv") 
Data.avg <- merge(Data.avg, Surface.Area,  by="Sample.ID")

Data.avg$chla.ug.cm2 <- Data.avg$chla.ug / Data.avg$SA.cm2  
Data.avg$chlc2.ug.cm2 <- Data.avg$chlc2.ug / Data.avg$SA.cm2 
Data.avg$group <- paste0(Data.avg$Species,"_",Data.avg$Site)
boxplot(Data.avg$chla.ug.cm2 ~ Data.avg$group, las=2)

