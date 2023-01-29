##########################################################################
##########################################################################
## FLIM analysis two samples
## Code from "Spatial organization of lysosomal exocytosis relies on membrane tension gradients", please cite Lachuer et al. PNAS (2023)
## Written by Hugo LACHUER
## R version 4.1.2 (2021-11-01) -- "Bird Hippie"
##########################################################################
##########################################################################

##########################################################################
##1) Librairies

library(raster)
library(viridis)
library(ggplot2)
library(ape)
library(spatstat)
library(sfsmisc)

##########################################################################
##2) Import functions and dataset

source("Functions.R")
colfuncGreen <- colorRampPalette(c("black", "Green"))
colfuncRed <- colorRampPalette(c("black", "Red"))

setwd("")
Ctrl <- list()	#list of FLIM images (raster transformed in matrix) (NA when no lifetime/outisde the cell) (fluoresence lifetime in ps) (Condition 1)
Osmo <- list()	#list of FLIM images (raster transformed in matrix) (NA when no lifetime/outisde the cell) (fluoresence lifetime in ps) (Condition 2)


PixelSize <- 0.087	#In µm

##########################################################################
##3) Average properties

#Computation
Cond1Average <- Cond1SD <- vector()
for(i in 1:length(Cond1)){
	Cond1Average[i] <- mean(Cond1[[i]], na.rm=TRUE)
	Cond1SD[i] <- sd(Cond1[[i]], na.rm=TRUE)
}
Cond2Average <- Cond2SD <- vector()
for(i in 1:length(Cond2)){
	Cond2Average[i] <- mean(Cond2[[i]], na.rm=TRUE)
	Cond2SD[i] <- sd(Cond2[[i]], na.rm=TRUE)
}

#Average
V1 <- c(Cond1Average, Cond2Average)
V2 <- c(rep("Cond1", length(Cond1Average)), rep("Cond2", length(Cond2Average)))
DF <- cbind(V1, V2)
DF <- as.data.frame(DF)
DF[,1] <- V1
DF[,2] <- factor(DF[,2], levels = c("Cond1", "Cond2"))
p.value <- t.test(Cond1Average, Cond2Average)$p.value
p <- ggplot(DF, aes(x=DF[,2], y=DF[,1], fill="Gray")) + geom_boxplot(outlier.shape = NA, fill='#A4A4A4', color="black") + geom_jitter(color="black", size=2, alpha=1, width = 0.25)
p <- p + theme(axis.text.x = element_text(size=18))+theme(axis.title.x = element_text(size=18))+theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))
p <- p + theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
p <- p + theme(plot.title = element_text(size=22)) + theme(legend.position = "none")
p <- p + ggtitle(paste0("Lifetime (p=", round(p.value,4),")")) + xlab("Condition") + ylab("Average lifetime (ps)") + theme(plot.title = element_text(size=22))
p



#SD
V1 <- c(Cond1SD, Cond2SD)
V2 <- c(rep("Cond1", length(Cond1Average)), rep("Cond2", length(Cond2Average)))
DF <- cbind(V1, V2)
DF <- as.data.frame(DF)
DF[,1] <- V1
DF[,2] <- factor(DF[,2], levels = c("Cond1", "Cond2"))
p.value <- t.test(Cond1SD, Cond2SD)$p.value
p <- ggplot(DF, aes(x=DF[,2], y=DF[,1], fill="Gray")) + geom_boxplot(outlier.shape = NA, fill='#A4A4A4', color="black") + geom_jitter(color="black", size=2, alpha=1, width = 0.25)
p <- p + theme(axis.text.x = element_text(size=18))+theme(axis.title.x = element_text(size=18))+theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))
p <- p + theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
p <- p + theme(plot.title = element_text(size=22)) + theme(legend.position = "none")
p <- p + ggtitle(paste0("Standard deviation (p=", round(p.value,4),")")) + xlab("Condition") + ylab("Standard-deviation of the lifetime (ps)") + theme(plot.title = element_text(size=22))
p



#Variation coefficient
V1 <- c(Cond1SD/Cond1Average, Cond2SD/Cond2Average)
V2 <- c(rep("Cond1", length(Cond1Average)), rep("Cond2", length(Cond2Average)))
DF <- cbind(V1, V2)
DF <- as.data.frame(DF)
DF[,1] <- V1
DF[,2] <- factor(DF[,2], levels = c("Cond1", "Cond2"))
p.value <- t.test(Cond1SD/Cond1Average, Cond2SD/Cond2Average)$p.value
p <- ggplot(DF, aes(x=DF[,2], y=DF[,1], fill="Gray")) + geom_boxplot(outlier.shape = NA, fill='#A4A4A4', color="black") + geom_jitter(color="black", size=2, alpha=1, width = 0.25)
p <- p + theme(axis.text.x = element_text(size=18))+theme(axis.title.x = element_text(size=18))+theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))
p <- p + theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
p <- p + theme(plot.title = element_text(size=22)) + theme(legend.position = "none")
p <- p + ggtitle(paste0("Variation coefficient (p=", round(p.value,4),")")) + xlab("Condition") + ylab("Variation coefficient") + theme(plot.title = element_text(size=22))
p

##########################################################################
##4) Fluoresence lifetime distribution


N <- nrow(Cond1[[1]])*ncol(Cond1[[1]])
X <- 2000:8000

#Cond1
V1 <- unlist(Cond1)
V2 <- as.character(sort(rep(1:length(Cond1), N)))
DF <- as.data.frame(cbind(V2,V1))
DF[,2] <- V1
Index <- which(is.na(V1))
DF <- DF[-Index,]
A1 <- AverageDensity(Data=DF, X=X, Border=FALSE, Circular=FALSE)
AverageDensityCond1 <- A1[[1]]
SEMDensityCond1 <- A1[[2]]
Indices <- which(is.na(SEMDensityCond1))
if(length(Indices)==0){Indices <- length(SEMDensityCond1)+1}
plot(x=X, y=AverageDensityCond1, type="l", main= "", ylab="Density", xlab="Lifetime (ps)", col="red", lwd=2, cex.lab=1.5, cex.main=2, cex.axis = 1.2)
xshade <- c(X[-Indices],rev(X[-Indices]))
yshade <- c(AverageDensityCond1[-Indices]-SEMDensityCond1[-Indices],rev(AverageDensityCond1[-Indices]+SEMDensityCond1[-Indices]))
polygon(xshade,yshade, border = NA, col=adjustcolor("red", 0.4))
#Single cell
plot(density(Cond1[[1]], na.rm=TRUE), col=inferno(length(Cond1))[1], lwd=2, xlim=c(2000,6000), ylim=c(0,0.003))
for(i in 2:length(Cond1)){
	lines(density(Cond1[[i]], na.rm=TRUE), col=inferno(length(Cond1))[i], lwd=2)
}

#Cond2
V1 <- unlist(Cond2)
V2 <- as.character(sort(rep(1:length(Cond2), N)))
DF <- as.data.frame(cbind(V2,V1))
DF[,2] <- V1
Index <- which(is.na(V1))
DF <- DF[-Index,]
A1 <- AverageDensity(Data=DF, X=X, Border=FALSE, Circular=FALSE)
AverageDensityCond2 <- A1[[1]]
SEMDensityCond2 <- A1[[2]]
Indices <- which(is.na(SEMDensityCond2))
if(length(Indices)==0){Indices <- length(SEMDensityCond2)+1}
plot(x=X, y=AverageDensityCond2, type="l", main= "", ylab="Density", xlab="Lifetime (ps)", col="red", lwd=2, cex.lab=1.5, cex.main=2, cex.axis = 1.2)
xshade <- c(X[-Indices],rev(X[-Indices]))
yshade <- c(AverageDensityCond2[-Indices]-SEMDensityCond2[-Indices],rev(AverageDensityCond2[-Indices]+SEMDensityCond2[-Indices]))
polygon(xshade,yshade, border = NA, col=adjustcolor("red", 0.4))
#Single cell
plot(density(Cond2[[1]], na.rm=TRUE), col=inferno(length(Cond2))[1], lwd=2, xlim=c(2000,6000), ylim=c(0,0.003))
for(i in 2:length(Cond2)){
	lines(density(Cond2[[i]], na.rm=TRUE), col=inferno(length(Cond2))[i], lwd=2)
}


#Pool
Ymax <- max(c(AverageDensityCond1, AverageDensityCond2), na.rm=TRUE)
plot(x=X, y=AverageDensityCond1, type="l", main= "", ylab="Density", xlab="Lifetime (ps)", col="red", lwd=2, cex.lab=1.5, cex.main=2, cex.axis = 1.2, ylim=c(0,Ymax))
lines(x=X, y=AverageDensityCond2, col="blue", lwd=2)
Indices <- which(is.na(SEMDensityCond1))
if(length(Indices)==0){Indices <- length(SEMDensityCond1)+1}
xshade1 <- c(X[-Indices],rev(X[-Indices]))
yshade1 <- c(AverageDensityCond1[-Indices]-SEMDensityCond1[-Indices],rev(AverageDensityCond1[-Indices]+SEMDensityCond1[-Indices]))
polygon(xshade1,yshade1, border = NA, col=adjustcolor("red", 0.4))
Indices <- which(is.na(SEMDensityCond2))
if(length(Indices)==0){Indices <- length(SEMDensityCond2)+1}
xshade2 <- c(X[-Indices],rev(X[-Indices]))
yshade2 <- c(AverageDensityCond2[-Indices]-SEMDensityCond2[-Indices],rev(AverageDensityCond2[-Indices]+SEMDensityCond2[-Indices]))
polygon(xshade2,yshade2, border = NA, col=adjustcolor("blue", 0.4))
legend(x=6000,y=6*10^-4, lwd=2, legend=c("Cond1", "Cond2"), col=c("red", "blue"))


##########################################################################
##5) Moran's index


#Cond1
MoranCond1 <- MoranIndex(Cond1, Bin=5, Sampling=TRUE, Nmax=1000)	
MoranCond1Matrix <- matrix(,ncol=4,nrow=length(MoranCond1))
colnames(MoranCond1Matrix) <- c("I observed", "I expected", "sd", "pvalue")
for(i in 1:length(MoranCond1)){
	MoranCond1Matrix[i,] <- as.vector(as.numeric(MoranCond1[[i]]))
}

#Cond2
MoranCond2 <- MoranIndex(Cond2, Bin=5, Sampling=TRUE, Nmax=1000)	
MoranCond2Matrix <- matrix(,ncol=4,nrow=length(MoranCond2))
colnames(MoranCond2Matrix) <- c("I observed", "I expected", "sd", "pvalue")
for(i in 1:length(MoranCond2)){
	MoranCond2Matrix[i,] <- as.vector(as.numeric(MoranCond2[[i]]))
}

#Plot
V1 <- c(MoranCond1[,1], MoranCond2[,1])
V2 <- c(rep("Cond1", nrow(MoranCond1)), rep("Cond2", nrow(MoranCond2)))
DF <- cbind(V1, V2)
DF <- as.data.frame(DF)
DF[,1] <- V1
DF[,2] <- factor(DF[,2], levels = c("Cond1", "Cond2"))
p.value <- t.test(MoranCond1[,1], MoranCond2[,1])$p.value
p <- ggplot(DF, aes(x=DF[,2], y=DF[,1], fill="Gray")) + geom_boxplot(outlier.shape = NA, fill='#A4A4A4', color="black") + geom_jitter(color="black", size=2, alpha=1, width = 0.25)
p <- p + theme(axis.text.x = element_text(size=18))+theme(axis.title.x = element_text(size=18))+theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))
p <- p + theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
p <- p + theme(plot.title = element_text(size=22)) + theme(legend.position = "none")
p <- p + ggtitle(paste0("Moran's Index (p=", round(p.value,4),")")) + xlab("Condition") + ylab("Moran's index") + theme(plot.title = element_text(size=22))
p