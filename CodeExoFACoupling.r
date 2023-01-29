##########################################################################
##########################################################################
## Colocalization/co-appearance between focal adhesions and exocytosis events
## Code from "Spatial organization of lysosomal exocytosis relies on membrane tension gradients", please cite Lachuer et al.
## Written by Hugo LACHUER
## R version 4.1.2 (2021-11-01) -- "Bird Hippie"
##########################################################################
##########################################################################


##########################################################################
##1) Librairies

library(raster)
library(spatstat)
library(ggplot2)

##########################################################################
##2) Import functions and dataset

#function
source("Functions.R")
colfuncGray <- colorRampPalette(c("white", "black"))

#Dataset
setwd()
Paxillin <- list()	#list of FA images (e.g. paxillin staining) (raster transformed in matrix)
Bin <- list()	#Mask of the cell (raster transformed in matrix)
Result <- list()	#Exocytosis coordinates (X and Y coordinates) (each element of the list is one cell)


##########################################################################
##3) Colocalization analysis

Nsim <- 100
pvalue <- Ratio <- vector(,length=length(Paxillin))

for(i in 1:length(Paxillin)){
	
	FA <- Paxillin[[i]]
	A1 <- Result[[i]]
	Nexo <- nrow(A1)

	plot(raster(FA), col=colfuncGray(100))
	points(x=A1[,1]/ncol(FA), y=abs(nrow(FA)-A1[,2])/nrow(FA), col="red", pch=19)
	
	Obs <- MeasureMap(Map=FA, X=A1[,1], Y=A1[,2])
	Sim <- vector(, length=Nsim)
	
	for(j in 1:Nsim){
		A2 <- CSRCell(n=nrow(A1), Temporal=FALSE, Shape="Free", Mask=Bin[[i]])
		Sim[j] <- MeasureMap(Map=FA, X=A2[,1], Y=A2[,2])
	}
	
	pvalue[i] <- (length(which(Sim>Obs))+1)/(Nsim+1)
	Ratio[i] <- Obs/mean(Sim)
}

##########################################################################
##4) p-value histogram

hist(pvalue, xlim=c(0,1), breaks=seq(0,1,by=0.1), main = "Exact test colocalization testing", ylab="Count", xlab="p-value", col="firebrick", lwd=4, cex.lab=1.5, cex.main=2, cex.axis = 1.2)
ks.test(pvalue,"punif",0,1)

##########################################################################
##5) Colocalization analysis

V1 <- rep("Ctrl",length=length(Ratio))
V2 <- Ratio
DF <- as.data.frame(cbind(V1,V2))
DF[,2] <- V2
p.value <- t.test(x=V2,mu=1, alternative="two.sided")$p.value
p <- ggplot(DF, aes(x=DF[,1], y=DF[,2], fill="Gray")) + geom_boxplot(outlier.shape = NA, fill='#A4A4A4', color="black") + geom_jitter(color="black", size=2, alpha=1, width = 0.25)
p <- p + theme(axis.text.x = element_text(size=18))+theme(axis.title.x = element_text(size=18))+theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))
p <- p + theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
p <- p + theme(plot.title = element_text(size=22)) + theme(legend.position = "none")+geom_hline(yintercept=1, linetype="dashed", color = "red", size=2)
p <- p + ggtitle(paste0("p=",round(p.value,4))) + xlab("") + ylab("Colocalization index") + theme(plot.title = element_text(size=22))
p