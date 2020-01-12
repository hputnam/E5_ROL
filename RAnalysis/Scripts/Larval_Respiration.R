#Title: Larval Respiration Calculations
#Author: KH Wong
#Date Last Modified: 20191029
#See Readme file for details

rm(list=ls()) #clears workspace 

## install packages if you dont already have them in your library
if ("devtools" %in% rownames(installed.packages()) == 'FALSE') install.packages('devtools') 
library(devtools)
if ("segmented" %in% rownames(installed.packages()) == 'FALSE') install.packages('segmented') 
if ("plotrix" %in% rownames(installed.packages()) == 'FALSE') install.packages('plotrix') 
if ("gridExtra" %in% rownames(installed.packages()) == 'FALSE') install.packages('gridExtra') 
if ("LoLinR" %in% rownames(installed.packages()) == 'FALSE') install_github('colin-olito/LoLinR') 
if ("lubridate" %in% rownames(installed.packages()) == 'FALSE') install.packages('lubridate') 
if ("chron" %in% rownames(installed.packages()) == 'FALSE') install.packages('chron') 
if ("plyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('plyr') 
if ("dplyr" %in% rownames(installed.packages()) == 'FALSE') install.packages('dplyr') 
if ("stringr" %in% rownames(installed.packages()) == 'FALSE') install.packages('stringr') 

#Read in required libraries
##### Include Versions of libraries
#install_github('colin-olito/LoLinR')
library("ggplot2")
library("segmented")
library("plotrix")
library("gridExtra")
library("LoLinR")
library("lubridate")
library("chron")
library('plyr')
library('dplyr')
library('stringr')
library('Rmisc')

## Data on 20101029 ##
path.p<-"RAnalysis/Data/Larval_Respirometry/20191029" #the location of all your respirometry files 

# bring in the respiration file names
file.names<-basename(list.files(path = path.p, pattern = "csv$", recursive = TRUE)) #list all csv file names in the folder and subfolders
#basename above removes the subdirectory name from the file
#add file names that include the subdirectory name
#file.names.full<- list.files(path = path.p, pattern = "csv$", recursive = TRUE) #list all csv file names in the folder and subfolders

#generate a 3 column dataframe with specific column names
Resp.R <- data.frame(matrix(NA, ncol=6))
colnames(Resp.R) <- c("Date", "Run", "Sample.ID","Chamber.ID","Intercept", "umol.L.sec")

Resp.Rb <- data.frame(matrix(NA, ncol=6))
colnames(Resp.Rb) <- c("Date", "Run", "Sample.ID","Chamber.ID","Intercept", "umol.L.sec")

#Load Sample Info
Sample.Info <- read.csv(file="RAnalysis/Data/Larval_Respirometry/Larval_Resp_Sample_Info.csv", header=T) #read sample.info data
#includes treatment, tank, chamber volume, animal size/number etc for normalization
# Sample.Info.20181102 <- Sample.Info %>%
#   filter(Date == "20181102")
Actual.T <- unique(Sample.Info$Target.T)
rename <- Sample.Info$Chamber.ID
samp <- Sample.Info$Sample.ID
run <- str_sub(file.names, 10, 13)#file.names
date <- str_sub(file.names, 1, str_length(file.names)-9)#file.names

for(i in 1:length(file.names)) { # for every file in list start at the first and run this following function
  Resp.Data <-read.table(file.path(path.p,file.names[i]), skip = 72, header=T, sep=",", na.string="NA", fill = TRUE, as.is=TRUE, fileEncoding="latin1") #reads in the data files
  Resp.Data$Temp <- Resp.Data[,31]
  Resp.Data$Time.Min <- seq.int(0.017, (nrow(Resp.Data))*0.25, by = 0.25) #set time in min
  Resp.Data.T <- Resp.Data %>% #filters data by temperature
     filter(Temp >= Actual.T[i] - 0.5 & Temp <= Actual.T[i] + 0.5)
  Resp.Data.N <- Resp.Data.T[,3:26]
  # Resp.Data.N <- Resp.Data[,3:26]
  
  for(j in 1:(ncol(Resp.Data.N))){
    model <- rankLocReg(
      xall=Resp.Data.T$Time.Min, yall=as.numeric(Resp.Data.N[, j]),
      alpha=0.4, method="pc", verbose=TRUE)
    
    pdf(paste0("RAnalysis/Output/Larval_Resp_Plots/20191029/",date[i], "_",run[i],"_",rename[j],"_regression_trunc.pdf"))
    plot(model)
    dev.off()
    
    Resp.Rb[j,1] <- as.character(date[i]) #stores the date
    Resp.Rb[j,2] <- as.character(run[i]) #stores the run number
    Resp.Rb[j,3] <- as.character(samp[j+(i-1)*ncol(Resp.Data.N)]) #stores the sample ID
    Resp.Rb[j,4] <- as.character(rename[j]) #stores the chamber ID
    Resp.Rb[j,5:6] <- model$allRegs[i,c(4,5)] #inserts slope and intercept in the dataframe
    
  }
  Resp.R <- rbind(Resp.R, Resp.Rb)
}

Resp.R <- Resp.R[-1,]

write.csv(Resp.R, paste0("RAnalysis/Output/Larval_Resp_rates_20191029.csv", sep=""))

Resp.R <- read.csv(file = "RAnalysis/Output/Larval_Resp_rates_20191029.csv")

#Merge data with sample info
Data <- merge(Resp.R, Sample.Info, by="Sample.ID")

#Account for chamber volume to convert from umol L-1 s-1 to umol s-1. This standardizes across water volumes (different because of coral size) and removes per Liter
Data$umol.sec <- Data$umol.L.sec * Data$Chamber.Vol.L

Data.BK <- subset(Data, Type == "Blank") #subset to blank data only
plot(Data.BK$umol.sec , ylab="umol O2 sec-1")

resp.blnk <- aggregate(umol.sec ~ Type*Target.T, data=Data.BK, mean)
Data$Blank <- resp.blnk$umol.sec[match(Data$Target.T,resp.blnk$Target.T)]
# 
# blank.means <- summarySE(Data.BK, measurevar="umol.sec", groupvars=c("Run.y"))
# colnames(blank.means)[which(names(blank.means) == "umol.sec")] <- "Blank.umol.sec"
# 
# Data.final <- merge(Data, blank.means, by = "Run.y")

#Account for blank rate Subtract Blank by the temperature blank
Data$umol.sec.corr <- Data$umol.sec-Data$Blank

#normalize to Org.number and min (and larval size eventually)

# sperm.conc <- read.csv("Data/CASA_Output2.csv")
# sperm.conc$Merge.ID <- paste0(as.character(sperm.conc$Date), "_", as.character(sperm.conc$Coral.ID))
# Data.Sample$Merge.ID <- paste0(as.character(Data.Sample$Date.x), "_", as.character(Data.Sample$Coral.ID))

# Data.Sperm <- merge(Data.Sample,sperm.conc, by="Merge.ID" )

Data.final <- Data %>% filter(Type == "Sample")

Data.final$umol.larva.sec <- Data.final$umol.sec.corr/Data.final$Org.Number
Data.final$umol.larva.min <- Data.final$umol.larva.sec/60
Data.final$nmol.larva.min <- Data.final$umol.larva.min*1000

Data.final$Target.T <- as.factor(Data.final$Target.T)

write.csv(Data.final, paste0("RAnalysis/Output/Larval_Resp_rates_20191029_calc.csv", sep=""))

#Plotting
Box.Larval.resp <- ggplot(Data.final, aes(x=Target.T, y=nmol.larva.min))+
  geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) 

ggsave(file="RAnalysis/Output/Larval.Resp.20191029.pdf", Box.Larval.resp, width = 11, height = 11, units = c("in"))

t.test(nmol.larva.min ~ Target.T, data = Data.final)
