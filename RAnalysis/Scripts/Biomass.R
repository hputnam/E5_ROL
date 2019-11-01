#Title: Biomass Data
#Project: E5
#Author: HM Putnam 
#Edited by: HM Putnam
#Date Last Modified: 20191031
#See Readme file for details

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
#if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 

Data <- read.csv("RAnalysis/Data/Biomass_Data.csv")
Meta <- read.csv("RAnalysis/Data/Sample_Metadata.csv")
Surface.Area <- read.csv("RAnalysis/Data/E5surface_area.csv") 
Data <-merge(Data, Meta, by="Tube.ID")
Data <- merge(Data, Surface.Area,  by="Sample.ID")

Data$delta.mass.g <- Data$dry.mass.g - Data$pan.mass.g
Data$mass.g.ml <- Data$delta.mass.g / Data$Vol.added.ml
Data$mass.g <- Data$mass.g.ml * Data$Homog.Vol.ml.x
Data$mass.g.cm2 <- Data$mass.g / Data$SA.cm2

Data$group <- paste0(Data$Species,"_",Data$Site)
boxplot(Data$mass.g.cm2 ~ Data$group, las=2)  