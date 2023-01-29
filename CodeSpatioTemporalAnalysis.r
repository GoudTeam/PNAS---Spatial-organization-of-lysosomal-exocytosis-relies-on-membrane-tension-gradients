##########################################################################
##########################################################################
## Spatial analysis of exocytosis
## Code from "Spatial organization of lysosomal exocytosis relies on membrane tension gradients", please cite Lachuer et al. PNAS (2023)
## Written by Hugo LACHUER
## R version 4.1.2 (2021-11-01) -- "Bird Hippie"
##########################################################################
##########################################################################

##########################################################################
##1) Librairies

library(ggpubr)
library(spatstat)
library(raster)
library(splancs)
library(plot3D)
library(imager)
library(evmix)
library(ggplot2)

##########################################################################
##2) Import functions and dataset

#Import functions
source("Functions.R")

#Import dataset
setwd("")

load("Ctrl - Cell masks.RData")	#BinList
BinList <- Mask
load("Ctrl - Exocytosis coordinates.RData")	#Exo
Exo <- Exo2
load("Ctrl - Acquisition time.RData")	#TIMER
TIMER <- Timer2
load("Ctrl - Number of frames.RData")	#NFrame
load("Ctrl - Exocytosis rate.RData")	#Frequency

PixelSize <- 0.160	#In µm

##########################################################################
##3) Correlation Ripley's K AUC and exocytosis frequency

AUC <- vector()

for(i in 1:length(Exo)){

	#Create events matrix
	Result <- Exo[[i]]
	Result[,1] <- round(Result[,1])
	Result[,2] <- round(Result[,2])
	
	#Mask
	Bin2 <- BinList[[i]]
	
	#Create spatstat object
	w <- owin(xrange=c(1,ncol(Bin2)*0.160), yrange=c(0,nrow(Bin2)*0.160), mask=(Bin2==1))	#Il faut une matrice booléenne
	p <- ppp(x=Result[,1]*0.160, y=Result[,2]*0.160, window=w)
		
	#Compute Ripley's AUC at 8µm
	A3 <- Kest(p, correction="best", rmax=8)
	A4 <- 0
	Deltar <- A3$r[2] - A3$r[1]	#Deltar is a constant
	for(j in 1:length(A3$r)){
		A4 <- A4 + A3$trans[j]*Deltar
	}
	AUC[i] <- A4
}


##Correlation
LR <- lm(Frequency ~ AUC)
summary(LR)
plot(x=AUC,y=Frequency, pch=19, col="blue",xlab="AUC",ylab="Frequency (event/µm2/s-1)")
LinearFit <- function(x)	{return(as.numeric(LR$coefficients[2])*x+as.numeric(LR$coefficients[1]))}
curve(LinearFit, 0.8*min(AUC,na.rm=TRUE), 1.2*max(AUC,na.rm=TRUE), add=TRUE, lwd=3, col="red", lty="dashed")


#Exocytosis rate
V1 <- rep("Ctrl",length=length(Exo))
V2 <- Frequency
DF <- as.data.frame(cbind(V1,V2))
DF[,2] <- V2
p <- ggplot(DF, aes(x=DF[,1], y=DF[,2], fill="Gray")) + geom_boxplot(outlier.shape = NA, fill='#A4A4A4', color="black") + geom_jitter(color="black", size=2, alpha=1, width = 0.25)
p <- p + theme(axis.text.x = element_text(size=18))+theme(axis.title.x = element_text(size=18))+theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))
p <- p + theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
p <- p + theme(plot.title = element_text(size=22)) + theme(legend.position = "none")
p <- p + ggtitle("") + xlab("") + ylab("event/µm²/s") + theme(plot.title = element_text(size=22))
p


##########################################################################
##4) NND

pvalue <- Clustering <- NND_MC <- B <- vector()

pb <- winProgressBar(title = "progress bar", min = 1,max = length(Exo), width = 300)	#progression bar
for(i in length(Exo)){

	#Create events matrix
	Result <- Exo[[i]]
	Result[,1] <- round(Result[,1])
	Result[,2] <- round(Result[,2])
	
	#Mask
	Bin2 <- BinList[[i]]
	
	Out <- NNDtest(X=Result[,1], Y=Result[,2], N_CSR_Sim=500, Shape="Free", R, Xlength, Ylength, Mask=Bin2)
	pvalue[i] <- Out[[1]]
	Clustering[i] <- Out[[2]]
	NND_MC[i] <- mean(Out[[3]])
	B[i] <- Out[[4]]
	
	setWinProgressBar(pb, i, title=paste0("Computation:", round(i/length(Exo)*100,2),"%"))
}
close(pb)

##p-value histogram
hist(pvalue, xlim=c(0,1), breaks=seq(0,1,by=0.1), main = "p-value histogram of CSR test (Monte-Carlo)", ylab="Count", xlab="p-value", col="cornflowerblue", lwd=4, cex.lab=1.5, cex.main=2, cex.axis = 1.2)
mtext(paste0("Clustering coefficient :", round(sum(Clustering)/length(pvalue),3), "        Dispersion coefficient :", round(1-sum(Clustering)/length(pvalue),3), "        p-valueKS =", round(ks.test(pvalue,"punif",0,1)$p,3)), cex=1.3)


##Absolute NND plot
V1 <- rep("Ctrl",length=length(B))
DF <- as.data.frame(cbind(V1,B))
DF[,2] <- B*PixelSize
colnames(DF) <- c("Condition","NND")
p <- ggplot(DF, aes(x=DF[,1], y=DF[,2], fill="Gray")) + geom_boxplot(outlier.shape = NA, fill='#A4A4A4', color="black") + geom_jitter(color="black", size=2, alpha=1, width = 0.25)
p <- p + theme(axis.text.x = element_text(size=18))+theme(axis.title.x = element_text(size=18))+theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))
p <- p + theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
p <- p + theme(plot.title = element_text(size=22)) + theme(legend.position = "none")
p <- p + ggtitle("") + xlab("") + ylab("NND (µm)") + theme(plot.title = element_text(size=22))
p


##Normalized NND plot
V1 <- rep("Ctrl",length=length(B))
DF <- as.data.frame(cbind(V1,B/NND_MC))
DF[,2] <- B/NND_MC
colnames(DF) <- c("Condition","NND")
p.value <- t.test(x=B, y=NND_MC, alternative="two.sided", paired=FALSE)$p.value
p <- ggplot(DF, aes(x=DF[,1], y=DF[,2], fill="Gray")) + geom_boxplot(outlier.shape = NA, fill='#A4A4A4', color="black") + geom_jitter(color="black", size=2, alpha=1, width = 0.25)
p <- p + theme(axis.text.x = element_text(size=18))+theme(axis.title.x = element_text(size=18))+theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))
p <- p + theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
p <- p + theme(plot.title = element_text(size=22)) + theme(legend.position = "none")
p <- p + ggtitle(paste0("p.value=",round(p.value,4))) + xlab("") + ylab("Relative NND") + theme(plot.title = element_text(size=22))
p <- p + geom_hline(yintercept=1, linetype="dashed", color = "red", size=2)
p


##########################################################################
##5) Ripley

N_sim <- 100

#Observation windows list
wmList <- ExoList <- list()
for(i in length(Exo)){
	
	#Mask
	Bin2 <- BinList[[i]]
	
	#Spatstat object
	w <- owin(xrange=c(1,ncol(Bin2)*PixelSize), yrange=c(0,nrow(Bin2)*PixelSize), mask=(Bin2==1))
	wmList[[i]] <- w
	
	Result  <- Exo[[i]]
	ExoList[[i]] <- Result[,1:2]*PixelSize
}

#Monte-Carlo CSR simulations
CSRList <- list()
for(i in 1:N_sim){
	CSRSubList <- list()
	for(j in length(Exo)){
		n_exo <- nrow(Exo[[j]])
		A1 <- CSRCell(n=n_exo, Temporal=FALSE, Shape="Free", Mask=wmList[[j]])
		CSRSubList[[j]] <- A1*PixelSize
	}
	CSRList[[i]] <- CSRSubList
}


#Average Ripley
A1 <- AverageRipley(DataList=ExoList, wmList=wmList, Inhom=FALSE, sigma=0)
R <- A1[[1]]
MuRipley <- A1[[2]]
SEMRipley <- A1[[3]]

#Error bands of Ripley
ErrorDensity <- matrix(,ncol=N_sim, nrow=length(R))
pb <- winProgressBar(title = "progress bar", min = 1,max = N_sim, width = 300)
for(i in 1:N_sim){
	A1 <- AverageRipley(DataList=CSRList[[i]], wmList=wmList, Inhom=FALSE, sigma=0)
	ErrorDensity[,i] <- A1[[2]]
	setWinProgressBar(pb, i, title=paste0("Computation:", round(i/N_sim*100,2),"%"))
}
close(pb)
QuantilesRipley <- matrix(,ncol=2,nrow=length(R))
colnames(QuantilesRipley) <- c("2.5% quantile", "97.5% quantile")
for(i in 1:length(R)){
	QuantilesRipley[i,1] <-  quantile(ErrorDensity[i,], prob = 0.025)
	QuantilesRipley[i,2] <-  quantile(ErrorDensity[i,], prob = 0.975)
}

#Plot Ripley
plot(x=R, y=MuRipley - pi*R^2, type="l", main= "", ylab=expression(paste("K(r)-", pi, "r²")), xlab="Distance r (µm)", col="red", lwd=2, cex.lab=1.5, cex.main=2, cex.axis = 1.2)
lines(x=R, y=rep(0,length(R)), col="blue", lty=3, lwd=3)
xshade1 <- c(R,rev(R))
yshade1 <- c((MuRipley- pi*R^2)-SEMRipley,rev((MuRipley- pi*R^2)+SEMRipley))
yshade2 <- c(QuantilesRipley[,1]- pi*R^2,rev(QuantilesRipley[,2]- pi*R^2))
polygon(xshade1,yshade1, border = NA, col=adjustcolor("red", 0.4))
polygon(xshade1,yshade2, border = NA, col=adjustcolor("blue", 0.4))


##########################################################################
##6) tRipley

N_sim <- 100

##Data preparation
TmaxVector <- DeltaTVector <- vector()
TemporalData <- list()
for(i in length(Exo)){
	A1 <- Exo[[i]]
	A1 <- A1[, which(colnames(A1)=="time" | colnames(A1)=="Slice")]
	TemporalData[[i]] <- (A1-1)*TIMER[i]/(NFrame[i]-1)
	TmaxVector[i] <- TIMER[i]
	DeltaTVector[i] <- TIMER[i]/(NFrame[i]-1)
}

#Monte-Carlo CSR simulations
CSRList <- list()
for(i in 1:N_sim){
	CSRSubList <- list()
	for(j in length(Exo)){
		n_exo <- nrow(Exo[[j]])
		A1 <- sample(seq(0, TmaxVector[j], DeltaTVector[j]), n_exo, replace=TRUE)
		CSRSubList[[j]] <- A1
	}
	CSRList[[i]] <- CSRSubList
}

#Average Ripley
A1 <- AverageTemporalRipley(DataList=TemporalData, DeltaTVector=DeltaTVector, TmaxVector=TmaxVector, Centered=TRUE, Standardized=FALSE)
R <- A1[[1]]
MuRipley <- A1[[2]]
SEMRipley <- A1[[3]]

#Error bands of Ripley
ErrorDensity <- matrix(,ncol=N_sim, nrow=length(R))
for(i in 1:N_sim){
	A1 <- AverageTemporalRipley(DataList=CSRList[[i]], DeltaTVector=DeltaTVector, TmaxVector=TmaxVector, Centered=TRUE, Standardized=FALSE)
	ErrorDensity[,i] <- A1[[2]]
}
QuantilesRipley <- matrix(,ncol=2,nrow=length(R))
colnames(QuantilesRipley) <- c("2.5% quantile", "97.5% quantile")
for(i in 1:length(R)){
	QuantilesRipley[i,1] <-  quantile(ErrorDensity[i,], prob = 0.025)
	QuantilesRipley[i,2] <-  quantile(ErrorDensity[i,], prob = 0.975)
}

#Plot Ripley
plot(x=R, y=MuRipley, type="l", main= "", ylab="K(t)-2t", xlab="Time t (s)", col="red", lwd=2, cex.lab=1.5, cex.main=2, cex.axis = 1.2)
lines(x=R, y=rep(0,length(R)), col="blue", lty=3, lwd=3)
xshade1 <- c(R,rev(R))
yshade1 <- c(MuRipley-SEMRipley,rev(MuRipley+SEMRipley))
yshade2 <- c(QuantilesRipley[,1],rev(QuantilesRipley[,2]))
polygon(xshade1,yshade1, border = NA, col=adjustcolor("red", 0.4))
polygon(xshade1,yshade2, border = NA, col=adjustcolor("blue", 0.4))

##########################################################################
##7) xytRipley

pvalueT <- matrix(,ncol=2,nrow=length(Exo))

pb <- winProgressBar(title = "progress bar", min = 1,max = length(Exo), width = 300)
for(i in length(Exo)){

	#Events matrix
	Result <- Exo[[i]]
	Result[,1] <- round(Result[,1])
	Result[,2] <- round(Result[,2])
	
	#Mask
	Bin2 <- BinList[[i]]
	Surface <- sum(Bin2)*PixelSize^2
	TypicalLenght <- sqrt(Surface/pi)	#TypicalSize of cell if it was disk-shaped
	
	#Splancs objetcs
	SpacePoints <- as.points(Result[,1]*PixelSize, Result[,2]*PixelSize)
	SpacePolygon <- MooreNeighborTracing(Bin2)*PixelSize
	TotalTime <- TIMER[i]
	DeltaT <- TotalTime/(NFrame[i]-1)
	TimePoints <- Result[,3]*DeltaT
	KST <- stkhat(pts=SpacePoints, times=TimePoints, poly=SpacePolygon, tlimits=c(DeltaT, TotalTime), s=seq(0,TypicalLenght,0.1), tm=seq(0,TotalTime/4,TotalTime/200))


	#plot3D
	x <- KST$s
	y <- KST$t
	z <- KST$kst 
	persp3D(x, y, z, theta=-45, phi=15, xlab="Space", ylab="Time", zlab="Spatio-temporal Ripley's K function")

	#Plot Kobs - pi.r².t
	CSRst <- function(x,t){KCSR <- 2*pi*x^2*t}
	z2 <- outer(x, y, CSRst)
	z3 <- z - z2
	persp3D(x, y, z3, theta = -35, phi = 25, expand = 0.5,  xlab="Space", ylab="Time", zlab="Kobs - 2πr²t")
	
	#Plot Kst - Ks*Kt
	Ks_Kt <- KST$ks%*%t(KST$kt)
	persp3D(x, y, KST$kst - Ks_Kt, theta = -35, phi = 25, expand = 0.5,  xlab="Space", ylab="Time", zlab="Kst - ks.kt")
	
	
	Simulated <- stmctest(pts=SpacePoints, times=TimePoints, poly=SpacePolygon, tlimits=c(DeltaT, TotalTime), s=seq(0,TypicalLenght,0.1), tt=seq(0,TotalTime/4,TotalTime/200), nsim=1000, quiet=TRUE, returnSims=FALSE)
	Nsup <- min(length(which(Simulated$t>=Simulated$t0)), length(which(Simulated$t<=Simulated$t0)))
	pvalue <- min(2*(Nsup+1)/(length(Simulated$t)+1),1)
	hist(Simulated$t, col="lightblue", lwd=3, xlab="U statistic", ylab="Count", main="Simulation histogram", xlim=c(1.2*min(unlist(Simulated)), 1.2*max(unlist(Simulated))))
	abline(v = Simulated$t0, col="red", lwd=3, lty=2)
	
	pvalueT[i,1] <- pvalue
	pvalueT[i,2] <- Simulated$t0
	setWinProgressBar(pb, i, title=paste0("Computation:", round(i/length(Exo)*100,2),"%"))
}
close(pb)

#p-value histogram
pvalue <- ks.test(pvalueT[,1], "punif", 0, 1)$p.value
hist(pvalueT[,1], col="lightblue", lwd=3, xlab="p-value", ylab="Frequency", main= paste0("Simulation histogram (pKS=", round(pvalue,4),")"),breaks=seq(0,1.02,0.1), freq=TRUE)

#Average statistic
V1 <- rep("Ctrl",length=length(Exo))
DF <- as.data.frame(cbind(V1,pvalueT[,2]))
DF[,2] <- pvalueT[,2]
p.value <- t.test(x=pvalueT[,2], mu=0, alternative="two.sided", paired=FALSE)$p.value
p <- ggplot(DF, aes(x=DF[,1], y=DF[,2], fill="Gray")) + geom_boxplot(outlier.shape = NA, fill='#A4A4A4', color="black") + geom_jitter(color="black", size=2, alpha=1, width = 0.25)
p <- p + theme(axis.text.x = element_text(size=18))+theme(axis.title.x = element_text(size=18))+theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))
p <- p + theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
p <- p + theme(plot.title = element_text(size=22)) + theme(legend.position = "none")
p <- p + ggtitle(paste0("p.value=",round(p.value,4))) + xlab("") + ylab("D") + theme(plot.title = element_text(size=22))
p <- p + geom_hline(yintercept=0, linetype="dashed", color = "red", size=2)+ ylim(-5*10^8, 5*10^8)
p


# Median curve
ExoList <- list()
dT <- T <- dX <- vector()
for(i in length(Exo)){

	#Events matrix
	Result <- Exo[[i]]
	Result[,1] <- round(Result[,1])
	Result[,2] <- round(Result[,2])
	ExoList[[i]] <- Result

	#Steps and boundaries
	dT[i] <- TIMER[i]/(NFrame[i]-1)
	T[i] <- TIMER[i]
	dX[i] <- PixelSize
}
A1 <- MedianSpatioTemporalRipley(ExoList=ExoList, dT=dT, T=T, RasterList=BinList, dX=dX)
x <- A1[[2]]
y <- A1[[3]]
z <- A1[[1]]
persp3D(x, y, z, theta = -35, phi = 25, expand = 0.5,  xlab="Space", ylab="Time", zlab="Kst-KsKt")


##########################################################################
##8) Cell border distribution

Nsim <- 100
pvalue <- Clustering <- NND_MC <- B <- vector()
colfuncGray <- colorRampPalette(c("white", "black"))
colfuncRed <- colorRampPalette(c("white", "blue", "red"))

AverageDCSR <- AverageD <- xmax <- NormalizeFactor <- rep(0,length(Exo))
DistributionD <- as.data.frame(matrix(,ncol=2,nrow=0))
CSRList <- list()
for(i in 1:Nsim){
	CSRList[[i]] <- DistributionD
}
pb <- winProgressBar(title = "progress bar", min = 1,max = length(Exo), width = 300)
for(i in length(Exo)){

	#Events matrix
	Result <- Exo[[i]]
	Result[,1] <- round(Result[,1])
	Result[,2] <- round(Result[,2])
	
	#Mask
	Bin2 <- BinList[[i]]
	plot(raster(Bin2), col=colfuncGray(2))
	points(x=Result[,1]/ncol(A2), y=abs(nrow(Bin2)-Result[,2])/nrow(A2), col="red", pch=19)
	
	
	#Distance Transform Map
	DistanceMap <- DistanceTransform(Bin2)
	DistanceMap <- DistanceMap*PixelSize
	plot(raster(DistanceMap), col=colfuncRed(100))
	points(x=Result[,1]/ncol(A2), y=abs(nrow(Bin2)-Result[,2])/nrow(A2), col="black", pch=19)
	NormalizeFactor[i] <- sum(DistanceMap)/sum(Bin2)
		
	#Remove events out of the cell
	for(j in nrow(Result):1){
		if(DistanceMap[Result[j,2], Result[j,1]] == 1){
			Result <- Result[-j,]
		}
	}
		
	#Distribution of distances
	DistanceVector <- vector(, length=nrow(Result))
	for(j in 1:nrow(Result)){
		DistanceVector[j] <- DistanceMap[Result[j,2], Result[j,1]]
	}
	AverageD[i] <- mean(DistanceVector)
	xmax[i] <- max(DistanceMap)
	B1 <- as.data.frame(cbind(rep(as.character(i),length(DistanceVector)),DistanceVector))
	B1[,2] <- DistanceVector
	DistributionD <- rbind(DistributionD,B1)
	
	#Monte-Carlo CSR simulations
	for(j in 1:Nsim){
		Coord <- CSRCell(n=nrow(Result), Temporal=FALSE, Shape="Free", Mask=Bin2)
		DistanceVectorCSR <- vector(, length=nrow(Coord))
		for(l in 1:nrow(Coord)){
			DistanceVectorCSR[l] <- DistanceMap[Coord[l,2], Coord[l,1]]
		}
		B1 <- as.data.frame(cbind(rep(as.character(i),length(DistanceVector)),DistanceVectorCSR))
		B1[,2] <- DistanceVectorCSR
		CSRList[[j]]<- rbind(CSRList[[j]],B1)
		AverageDCSR[i] <- AverageDCSR[i]+mean(DistanceVectorCSR)
	}
	AverageDCSR[i] <- AverageDCSR[i]/Nsim
	setWinProgressBar(pb, i, title=paste0("Computation:", round(i/length(Exo)*100,2),"%"))
}
close(pb)

#Average distribution
X <- seq(0,median(xmax),median(xmax)/100)
A1 <- AverageDensity(Data=DistributionD, X=X, Border=TRUE, Circular=FALSE, xmax=xmax)	#warnings because some points are at d=0 --> removed
Mu <- A1[[1]]
SEM <- A1[[2]]

#Error bands
ErrorDensity <- matrix(,ncol=N_sim, nrow=length(X))
for(i in 1:N_sim){
	A1 <- AverageDensity(Data=CSRList[[i]], X=X, Border=TRUE, Circular=FALSE, xmax=xmax)
	ErrorDensity[,i] <- A1[[1]]
}
Quantiles <- matrix(,ncol=2,nrow=length(X))
MuCSR <- vector(,length=length(X))
colnames(Quantiles) <- c("2.5% quantile", "97.5% quantile")
for(i in 1:length(X)){
	Quantiles[i,1] <-  quantile(ErrorDensity[i,], prob = 0.025)
	Quantiles[i,2] <-  quantile(ErrorDensity[i,], prob = 0.975)
	MuCSR[i] <- mean(ErrorDensity[i,])
}

#Plot
plot(x=X, y=Mu, type="l", main= "", ylab="Density", xlab="Distance from cell border", col="red", lwd=2, cex.lab=1.5, cex.main=2, cex.axis = 1.2, ylim=c(0,0.30))
lines(x=X, y=-2/(max(X)^2-max(X))*(X-max(X)), col="blue", lty=3, lwd=3)
lines(x=X, y=MuCSR, col="blue", lwd=3)
xshade1 <- c(X,rev(X))
yshade1 <- c(Mu-SEM,rev(Mu+SEM))
yshade2 <- c(Quantiles[,1],rev(Quantiles[,2]))
polygon(xshade1,yshade1, border = NA, col=adjustcolor("red", 0.4))
polygon(xshade1,yshade2, border = NA, col=adjustcolor("blue", 0.4))


#Average cell border distance
A1 <- unique(DistributionD[,1])
MuD <- vector()
for(i in 1:length(A1)){
	D <- DistributionD[which(DistributionD[,1]==A1[i]),2]
	MuD[i] <- mean(D)
}
V1 <- rep("Ctrl",length=length(Exo))
DF <- as.data.frame(cbind(V1,MuD))
DF[,2] <- MuD
p <- ggplot(DF, aes(x=DF[,1], y=DF[,2], fill="Gray")) + geom_boxplot(outlier.shape = NA, fill='#A4A4A4', color="black") + geom_jitter(color="black", size=2, alpha=1, width = 0.25)
p <- p + theme(axis.text.x = element_text(size=18))+theme(axis.title.x = element_text(size=18))+theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))
p <- p + theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
p <- p + theme(plot.title = element_text(size=22)) + theme(legend.position = "none")
p <- p + ggtitle("") + xlab("") + ylab("Distance to cell border (µm)") + theme(plot.title = element_text(size=22))
p

#Average normalized cell border distance
V1 <- rep("Ctrl",length=length(Exo))
DF <- as.data.frame(cbind(V1,MuD/AverageDCSR))
DF[,2] <- MuD/AverageDCSR
p.value <- t.test(x=MuD/AverageDCSR, mu=1, alternative="two.sided", paired=FALSE)$p.value
p <- ggplot(DF, aes(x=DF[,1], y=DF[,2], fill="Gray")) + geom_boxplot(outlier.shape = NA, fill='#A4A4A4', color="black") + geom_jitter(color="black", size=2, alpha=1, width = 0.25)
p <- p + theme(axis.text.x = element_text(size=18))+theme(axis.title.x = element_text(size=18))+theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))
p <- p + theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
p <- p + theme(plot.title = element_text(size=22)) + theme(legend.position = "none")
p <- p + ggtitle(paste0("p.value=",round(p.value,4))) + xlab("") + ylab("Normalized distance to cell border") + theme(plot.title = element_text(size=22))
p <- p + geom_hline(yintercept=1, linetype="dashed", color = "red", size=2)
p

##########################################################################
##9) Fourier


#Prepare Data
DeltaTVector <- vector()
SignalList <- list()
for(i in length(Exo)){
	A1 <- Exo[[i]]
	A1 <- A1[,3]
	A2 <- vector(length=NFrame[i])
	for(j in 1:NFrame[i]){
		A2[j] <- length(which(A1==j))
	}
	SignalList[[i]]	<- A2
	DeltaTVector[i] <- TIMER[i]/(NFrame[i]-1)
}


A1 <- AverageFourier(TemporalDataList=SignalList, DeltaTVector=DeltaTVector, Shuffling=FALSE, Normalization=TRUE)
f <- A1[[1]]
MuFourier <- A1[[2]]
SEMFourier <- A1[[3]]
Frequencies <- A1[[4]]

#Plot
plot(x=f*10^3, y=MuFourier, type="l", main= "", ylab="Modulus", xlab="Frequency (mHz)", col="red", lwd=1, cex.lab=1.5, cex.main=2, cex.axis = 1.2)
xshade1 <- c(f,rev(f))
yshade1 <- c(MuFourier-SEMFourier,rev(MuFourier+SEMFourier))
polygon(xshade1,yshade1, border = NA, col=adjustcolor("red", 0.4))

##########################################################################
##10) Anisotropy/Polarization


Npart <- 30
vec.breaks <- seq(from = 0, to = 2*pi, by = pi*0.2)
graduation <- c("0", "0.2*pi", "0.4*pi", "0.6*pi", "0.8*pi", "pi", "1.2*pi", "1.4*pi", "1.6*pi", "1.8*pi", "2*pi")
graduation2 <- c(0, 0.2*pi, 0.4*pi, 0.6*pi, 0.8*pi, pi, 1.2*pi, 1.4*pi, 1.6*pi, 1.8*pi, 2*pi)
vec.expr <- parse(text = graduation)

CSRPolarIndex <- PolarIndex <- Polarpvalue <- vector()

for(i in length(Exo)){

	#Events matrix
	Result <- Exo[[i]]
	Result[,1] <- round(Result[,1])
	Result[,2] <- round(Result[,2])
	
	#Mask
	Bin2 <- BinList[[i]]
	Center <- CenterMass(Bin2)
	
	#Generate radial coordinates of cell pixels
	Indices <- which(Bin2==1, arr.ind=TRUE)
	RadialCell <- vector()
	for(j in 1:nrow(Indices)){
		theta <- atan2(Indices[j,1]-Center[2], Indices[j,2]-Center[1]) + pi 	#Angle between 0 and 2pi
		RadialCell[j] <- theta
	}
	
	#Generate radial coordinates of exocytosis events
	RadialExo <- vector()
	for(j in 1:nrow(Result)){
		theta <- atan2(Result[j,2]-Center[2], Result[j,1]-Center[1]) + pi
		RadialExo[j] <- theta
	}

	#Angular distribution
	CellDistribution <- vector(,length=Npart)
	Angles <- seq(from=0,to=2*pi, length.out=Npart+1)
	for(j in 1:Npart){
		Surface <- length(which((RadialCell > Angles[j]) & (RadialCell < Angles[j+1])))
		Nexo <- length(which((RadialExo > Angles[j]) & (RadialExo < Angles[j+1])))
		CellDistribution[j] <- Nexo/Surface
	}
	
	#Quantification
	Angles2 <- Angles + 2*pi/(2*Npart)
	Angles2 <- Angles2[-(Npart+1)]
	PolarIndex[i] <- ComputPolarIndex(value=CellDistribution, angle=Angles2)
	A2 <- TestPolarExocytosis(Bin=Bin2, Coordinates=Result[,1:2], NPart=NPart, Nsim=300)
	Polarpvalue[i] <- A2[1]
	CSRPolarIndex[i] <- A2[2]
	
	#Angular plot
	DF <- as.data.frame(matrix(,ncol=2,nrow=Npart))
	DF[,1] <- CellDistribution
	DF[,2] <- Angles2
	colnames(DF) <- c("Surface density", "Angle")
	p <- ggplot(DF, aes(x=DF[,2], y=DF[,1], fill=DF[,1])) + geom_bar(stat="identity") +theme_light() + scale_fill_gradient(low="white", high="red", limits=c(0.9*min(DF[,1]),max(DF[,1]))) + theme(axis.title.y=element_text(angle=0))
	p <- p + theme(axis.text.x = element_text(angle=45, vjust = 1, hjust=1)) + coord_polar() + geom_hline(yintercept = seq(0, max(DF[,1]), by = max(DF[,1])/4), color = "black", size = 0.5)
	p <- p + theme(axis.text.x = element_text(size=18))+theme(axis.title.x = element_text(size=18))+theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))
	p <- p + theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
	p <- p + theme(plot.title = element_text(size=22)) + theme(legend.position = "none") + scale_x_continuous(breaks = graduation2, labels = vec.expr )
	p <- p + ggtitle(paste0("Polarization Index = ",round(PolarIndex[i],4)," pvalue=",round(Polarpvalue[i],4))) + xlab("") + ylab("") + theme(plot.title = element_text(size=22))
	p
}

#p-value histogram
hist(Polarpvalue, xlim=c(0,1), breaks=seq(0,1,by=0.1), main = "Exact test polarization testing", ylab="Count", xlab="p-value", col="firebrick", lwd=4, cex.lab=1.5, cex.main=2, cex.axis = 1.2)
ks.test(Polarpvalue,"punif",0,1)

#Average statistic
V1 <- rep("Ctrl",length=length(Exo))
DF <- as.data.frame(cbind(V1,PolarIndex))
DF[,2] <- PolarIndex
p <- ggplot(DF, aes(x=DF[,1], y=DF[,2], fill="Gray")) + geom_boxplot(outlier.shape = NA, fill='#A4A4A4', color="black") + geom_jitter(color="black", size=2, alpha=1, width = 0.25)
p <- p + theme(axis.text.x = element_text(size=18))+theme(axis.title.x = element_text(size=18))+theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))
p <- p + theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
p <- p + theme(plot.title = element_text(size=22)) + theme(legend.position = "none")
p <- p + ggtitle("") + xlab("") + ylab("Polarization index") + theme(plot.title = element_text(size=22))
p

#Average statistic (normalized)
V1 <- rep("Ctrl",length=length(Exo))
DF <- as.data.frame(cbind(V1,PolarIndex/CSRPolarIndex))
DF[,2] <- PolarIndex/CSRPolarIndex
p.value <- t.test(x=DF[,2], mu=1, alternative="two.sided", paired=FALSE)$p.value
p <- ggplot(DF, aes(x=DF[,1], y=DF[,2], fill="Gray")) + geom_boxplot(outlier.shape = NA, fill='#A4A4A4', color="black") + geom_jitter(color="black", size=2, alpha=1, width = 0.25)
p <- p + theme(axis.text.x = element_text(size=18))+theme(axis.title.x = element_text(size=18))+theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))
p <- p + theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
p <- p + theme(plot.title = element_text(size=22)) + theme(legend.position = "none")
p <- p + ggtitle(paste0("pvalue=",round(p.value,4))) + xlab("") + ylab("Relative polarization index") + theme(plot.title = element_text(size=22))
p <- p + geom_hline(yintercept=1, linetype="dashed", color = "red", size=2)
p

length(which(DF[,2]>1))/nrow(DF)	#Number of cells more polarized than expected
length(which(DF[,2]>2 & Polarpvalue<0.05))/length(which(Polarpvalue<0.05))	#Number of cells significantly more polarized than expected