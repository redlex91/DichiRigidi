# carico i dati
setwd("~/Documents/LabFisComp/DischiRigidi")
data<-read.table("speed.dat")
vx<-data[,1]
vy<-data[,2]

setwd("~/Documents/LabFisComp/DischiRigidi/R")
pdf(file="speed_dist.pdf")

# faccio i grafici
par(mfrow=c(2,1))
histvx<-hist(vx, breaks=50, freq=0, main="Distribuzione delle velocità lungo x", xlab="Velocità delle particelle lungo x")
par(new=TRUE)
lines(histvx$breaks, dnorm(histvx$breaks,mean(vx), sd(vx)), "l", col="red", lwd=3, main="", xlab="", ylab="")
histvy<-hist(vy, breaks=50, freq=0, main="Distribuzione delle velocità lungo y", xlab="Velocità delle particelle lungo y")
par(new=TRUE)
lines(histvy$breaks, dnorm(histvy$breaks,mean(vy), sd(vy)), "l", col="red", lwd=3, main="", xlab="", ylab="")

dev.off()

