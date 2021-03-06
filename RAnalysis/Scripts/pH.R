library(lubridate)

pH <- read.csv("RAnalysis/Data/20200102_SN20697597_Site1_pH.csv", header=F, sep=",", na.strings="NA", skip=700) #load  data
pH <- pH[,1:4]
colnames(pH) <- c("Date.Time", "Temp", "mV", "pH")
pH$Date.Time <- parse_date_time(pH$Date.Time, "Y!-m!*-d! H!:M!:S!") 

pdf("RAnalysis/Output/20200102_pH_Temp.pdf")
par(mfrow=c(2,1))
par(mar = c(5,5,2,5))
with(pH, plot(Date.Time, pH, type="l", col="blue", 
             ylim=c(7.9,8.5), ylab =NA))
axis(2, col="blue",col.axis="blue")
mtext("pH",side=2,col="blue",line=4) 
par(new = T)
with(pH, plot(Date.Time, Temp,type="l", axes=F,col="black", xlab=NA, ylab=NA, ylim=c(27,32)))
axis(side = 4)
mtext(side = 4, line = 3, 'Temperature°C')
#dev.off()

pH <- read.csv("RAnalysis/Data/20200103_SN20697596_Site2_pH.csv", header=F, sep=",", na.strings="NA", skip=700) #load  data
pH <- pH[,1:4]
colnames(pH) <- c("Date.Time", "Temp", "mV", "pH")
pH$Date.Time <- parse_date_time(pH$Date.Time, "Y!-m!*-d! H!:M!:S!") 

#pdf("RAnalysis/Output/20200103_Site2_pH_Temp.pdf")
par(mar = c(5,5,2,5))
with(pH, plot(Date.Time, pH, type="l", col="blue", 
              ylim=c(7.9,8.5), ylab =NA))
axis(2, col="blue",col.axis="blue")
mtext("pH",side=2,col="blue",line=4) 
par(new = T)
with(pH, plot(Date.Time, Temp,type="l", axes=F,col="black", xlab=NA, ylab=NA, ylim=c(27,32)))
axis(side = 4)
mtext(side = 4, line = 3, 'Temperature°C')
dev.off()