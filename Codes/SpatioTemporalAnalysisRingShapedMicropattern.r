##########################################################################
##########################################################################
## Spatio-temporal analysis Ring shaped micropattern
## Code from "Spatial organization of lysosomal exocytosis relies on membrane tension gradients", please cite Lachuer et al. PNAS (2023)
## Written by Hugo LACHUER
## R version 4.1.2 (2021-11-01) -- "Bird Hippie"
##########################################################################
##########################################################################

##########################################################################
##1) Librairies

library(circular)
library(DescTools)
library(spatstat)
library(ggplot2)
library(evmix)
library(raster)

##########################################################################
##2) Import functions and dataset


source("Functions.R")
setwd()

Cell <- list()	#Exocytosis coordinates (X,Y and T coordinates) (each element of the list is one cell)
Timer <- vector()	#Vector of movie durations in secondes
Geometry <- matrix(, ncol=7, nrow=length(Cell))	#Geometry of the pattern
colnames(Geometry) <- c("X_Center", "Y_Center", "CellRadius", "PatternRadius",  "DeltaRadius", "PatternThickness", "NormalizedThickness")
#X_center/Y_center correspond to the center coordinates of the cell
#CellRadius corresponds to the cell radius
#PatternRadius corresponds to the pattern radius
#DeltaRadius corresponds to PatternRadius - CellRadius
#PatternThickness corresponds to the thickness of the adhesive band of the ring-shaped micropattern
#Normalized Thickness corresponds to (PatternThickness-DeltaRadius)/CellRadius

#Unit disk (observation windows)
w <- owin(c(-1,1), c(-1,1), mask=matrix(TRUE, 1000,1000))
X <- raster.x(w)
Y <- raster.y(w)
wm <- owin(w$xrange, w$yrange, mask=(X^2 + Y^2 <= 1))
plot(wm)

PixelSize <- 0.160	#In µm (vector if not constant)

##########################################################################
##3) Normalize coordinates

Exo <- Cell

for(i in 1: length(Cell)){
	A1 <- Cell[[i]]
	A1[,1:2] <- A1[,1:2]*PixelSize	#Convert in µm
	
	#Change axis
	X <- Geometry[i,1]
	Y <- Geometry[i,2]
	R <- Geometry[i,3]
	A1[,1] <- A1[,1]-X
	A1[,2] <- (A1[,2]-Y)*(-1)

	#Compute polar coordinates
	Data <- as.data.frame(matrix(0, ncol=6, nrow=nrow(A1)))
	colnames(Data) <- c("time", "Theta", "A", "A normalized", "x", "y")
	Data[,1] <- A1[,3]
	Data[,2] <- atan2(A1[,2], A1[,1]) + pi #Angle between 0 and 2pi
	Data[,3] <- sqrt(A1[,1]^2 + A1[,2]^2)
	Data[,4] <- Data[,3]/R
	Data[,5] <- A1[,1]/R
	Data[,6] <- A1[,2]/R
	
	#Remove events when A>1
	k <- 0
	for(j in nrow(Data):1){
		if(Data[j,4] > 1){
			Data <- Data[-j,]
			k <- k + 1
		}
	}
	print(paste0(k, " points removed in cell ", i))

	Exo[[i]] <- Data
}

##########################################################################
##4) Single cell plot

for(i in 1:length(Exo)){
	A1 <- Exo[[i]]
	A1 <- A1[,1:2]
	pattern_radius <- Geometry[i,7]
	
	plot(x=A1[,1], y=A1[,2], xlim=c(-1.2, 1.2), ylim=c(-1.2, 1.2), pch=16, col="red", lwd=3, asp=1)
	DrawCircle(x = 0, y = 0, r.out = 1, r.in = 0, theta.1 = 0, theta.2 = 2*pi, border = par("fg"), lty = par("lty"), lwd = 3, nv = 1000, plot = TRUE)
	DrawCircle(x = 0, y = 0, r.out = 1, r.in = (1-pattern_radius), theta.1 = 0, theta.2 = 2*pi, border = NA, lty = par("lty"), lwd = 3, nv = 1000,col= rgb(0.06,0.306,0.545,alpha=0.40), plot = TRUE)

	#2D KDE
	p <- ppp(x=A1[,1], y=A1[,2], window=wm)
	ds <- density(p, kernel="gaussian", sigma=bw.diggle, edge=TRUE, diggle=TRUE)
	plot(ds ,ylab="Y coordinate", xlab="X coordinate", cex.lab=1.5, cex.main=2, cex.axis = 1.2 )
	points(p, pch=20)
	
	R <- A1[,1]^2 + A1[,2]^2
	Theta <- atan2(A1[,2], A1[,1])+pi
	
	LR <- lm(R ~ Theta)
	print(paste0("R² : ", summary(LR)$r.squared	))
	plot(x=Theta, y=R, type="p", col="blue", ylim=c(0,1), xlim=c(0, 2*pi), cex.lab=1.5, cex.main=2, lwd=3, pch=16, xlab="Angle", ylab="Normalized modulus")
	LinearFit <- function(x)	{return(as.numeric(LR$coefficients[2])*x+as.numeric(LR$coefficients[1]))}
	curve(LinearFit, 0, 2*pi, add=TRUE, lwd=3, col="red", lty="dashed")
}


##########################################################################
##5) Modification of the ExoList


NExo <- vector(,length=length(Exo))
for(i in 1:length(Exo)){
	A1 <- Exo[[i]]
	NExo[i] <- nrow(A1)
}
NExo


#Re-arrangement
for(i in 1:length(Exo)){
	A1 <- Exo[[i]]
	DeltaT <- TimerRing[i,1]/(TimerRing[i,2]-1)
	A1[,1] <- (A1[,1]-1)*DeltaT
	Exo[[i]] <- A1[,c(5,6,1)]
}

##########################################################################
##6) NND analysis

pvaluevector <- NNDRing <- vector()
Clustering <- Dispersing <- 0
for(i in 1:length(Exo)){
	A1 <- Exo[[i]]
	Output <- NNDtest(X=A1[,1], Y=A1[,2], N_CSR_Sim=100, Shape="Disk", R=1, Xlength=0, Ylength=0)
	d1 <- density(Output[[3]], bw="SJ", kernel = "gaussian")
	plot(d1, main= paste0("NND distribution (p-value=", round(Output[[1]],3), ")"), ylab="Density", xlab="NND", col="blue", lwd=4, cex.lab=1.5, cex.main=2, cex.axis = 1.2)
    abline(v = Output[[4]], col="red", lwd=3, lty=2)
	pvaluevector[i] <- Output[[1]]
	ifelse(Output[[2]]==TRUE, Clustering <- Clustering + 1, Dispersing <- Dispersing + 1)
	NNDRing[i] <- Output[[4]]
}

hist(pvaluevector, xlim=c(0,1), breaks=seq(0,1,by=0.1), main = "p-value histogram of CSR test (Monte-Carlo)", ylab="Count", xlab="p-value", col="cornflowerblue", lwd=4, cex.lab=1.5, cex.main=2, cex.axis = 1.2)
mtext(paste0("Clustering coefficient :", round(Clustering/length(pvaluevector),3), "        Dispersion coefficient :", round(Dispersing/length(pvaluevector),3), "        p-valueKS =", round(ks.test(pvaluevector,"punif",0,1)$p,3)), cex=1.3)


##########################################################################
##7) 2D KDE

#Ring
r <- 1 - mean(GeometryRing[,7])
Output <- Average2DKDE(CoordList=Exo, wm=wm)
plot(raster(Output[[1]]), col=inferno(256))	#Moyenne
DrawCircle(x = 0.5, y = 0.5, r.out = r/2, r.in = 0, border="white",  theta.1 = 0, theta.2 = 2*pi, lty = 2, lwd = 3, nv = 1000, plot = TRUE)	#Pattern
plot(raster(Output[[2]]), col=inferno(256))	#SEM
DrawCircle(x = 0.5, y = 0.5, r.out = r/2, r.in = 0, border="white",  theta.1 = 0, theta.2 = 2*pi, lty = 2, lwd = 3, nv = 1000, plot = TRUE)	#Pattern
plot(raster(Output[[3]]), col=inferno(256))	#Median (évite les outliers)
DrawCircle(x = 0.5, y = 0.5, r.out = r/2, r.in = 0, border="white",  theta.1 = 0, theta.2 = 2*pi, lty = 2, lwd = 3, nv = 1000, plot = TRUE)	#Pattern


##########################################################################
##8) 1D KDE

Xmodulus <- seq(0,1,0.01)
Xangle <- seq(0,2*pi,0.01)
N_sim <- 100
vec.breaks <- seq(from = 0, to = 2*pi, by = pi*0.2)
graduation <- c("0", "0.2*pi", "0.4*pi", "0.6*pi", "0.8*pi", "pi", "1.2*pi", "1.4*pi", "1.6*pi", "1.8*pi", "2*pi")
graduation2 <- c(0, 0.2*pi, 0.4*pi, 0.6*pi, 0.8*pi, pi, 1.2*pi, 1.4*pi, 1.6*pi, 1.8*pi, 2*pi)
vec.expr <- parse(text = graduation)

#Data preparation
MatrixRing <- as.data.frame(matrix(,ncol=3,nrow=0))
colnames(MatrixRing) <- c("r", "theta", "name")
for(i in 1:length(Exo)){
	A1 <- Exo[[i]]
	A2 <- matrix(,ncol=3,nrow=nrow(A1))
	A2[,1] <- sqrt(A1[,1]^2 + A1[,2]^2)
	A2[,2] <- atan2(A1[,2],A1[,1])+pi
	A2[,3] <- rep(paste0("Cell",i), length=nrow(A1))
	MatrixRing <- rbind(MatrixRing, A2)
}
MatrixRing[,1] <- as.numeric(as.vector(MatrixRing[,1]))
MatrixRing[,2] <- as.numeric(as.vector(MatrixRing[,2]))

#Monte-Carlo CSR simulation
CSRList <- list()
for(i in 1:N_sim){
	MatrixA2 <- as.data.frame(matrix(,ncol=3,nrow=0))
	for(j in 1:length(Exo)){
		n_exo <- nrow(Exo[[j]])
		A1 <- CSRCell(n_exo, Temporal=FALSE, Shape="Disk", R=1, Xlength=0, Ylength=0)
		A2 <- matrix(,ncol=3,nrow=nrow(A1))
		A2[,1] <- sqrt(A1[,1]^2 + A1[,2]^2)
		A2[,2] <- atan2(A1[,2],A1[,1])+pi
		A2[,3] <- rep(paste0("Cell",j), length=nrow(A1))
		MatrixA2 <- rbind(MatrixA2, A2)
	}
	MatrixA2[,1] <- as.numeric(as.vector(MatrixA2[,1]))
	MatrixA2[,2] <- as.numeric(as.vector(MatrixA2[,2]))
	CSRList[[i]] <- MatrixA2
}


#Average modulus density
A1 <- AverageDensity(Data=MatrixRing[,c(3,1)], X=Xmodulus, Border=TRUE, Circular=FALSE, xmax=1)
AverageDensityModulusRing <- A1[[1]]
SEMDensityModulusRing <- A1[[2]]

#Error bands of modulus ECDF
ErrorDensity <- matrix(,ncol=N_sim, nrow=length(Xmodulus))
for(i in 1:N_sim){
	A1 <- AverageDensity(Data=CSRList[[i]][,c(3,1)], X=Xmodulus, Border=TRUE, Circular=FALSE, xmax=1)
	ErrorDensity[,i] <- A1[[1]]
}
QuantilesDensity_ModulusRing <- matrix(,ncol=2,nrow=length(Xmodulus))
colnames(QuantilesDensity_ModulusRing) <- c("2.5% quantile", "97.5% quantile")
for(i in 1:length(Xmodulus)){
	QuantilesDensity_ModulusRing[i,1] <-  quantile(ErrorDensity[i,], prob = 0.025)
	QuantilesDensity_ModulusRing[i,2] <-  quantile(ErrorDensity[i,], prob = 0.975)
}

#Plot modulus density
plot(x=Xmodulus, y=AverageDensityModulusRing, type="l", main= "", ylab="Density", xlab="Normalized modulus", col="red", lwd=2, cex.lab=1.5, cex.main=2, cex.axis = 1.2, ylim=c(0,2))
lines(x=Xmodulus, y=2*Xmodulus, col="blue", lty=3, lwd=3)
xshade <- c(Xmodulus,rev(Xmodulus))
yshade1 <- c(AverageDensityModulusRing-SEMDensityModulusRing,rev(AverageDensityModulusRing+SEMDensityModulusRing))
yshade2 <- c(QuantilesDensity_ModulusRing[,1],rev(QuantilesDensity_ModulusRing[,2]))
polygon(xshade,yshade1, border = NA, col=adjustcolor("red", 0.4))
polygon(xshade,yshade2, border = NA, col=adjustcolor("blue", 0.4))

r <- 1-mean(GeometryRing[,7])
SEMr <- sd(1-GeometryRing[,7])/sqrt(length(GeometryRing[,7]))
rect(r, 0, 1, 1000, col= adjustcolor("darkgray",alpha=0.40), border = NA)
errorbar(r, 0.5, r-SEMr, r+SEMr, 0.1)


#Average modulus ECDF
A1 <- AverageECDF(Data=MatrixRing[,c(3,1)], X=Xmodulus)
AverageECDFModulusRing <- A1[[1]]
SEMECDFModulusRing <- A1[[2]]

#Error bands of modulus ECDF
ErrorECDF <- matrix(,ncol=N_sim, nrow=length(Xmodulus))
for(i in 1:N_sim){
	A1 <- AverageECDF(Data=CSRList[[i]][,c(3,1)], X=Xmodulus)
	ErrorECDF[,i] <- A1[[1]]
}
QuantilesECDF_ModulusRing <- matrix(,ncol=2,nrow=length(Xmodulus))
colnames(QuantilesECDF_ModulusRing) <- c("2.5% quantile", "97.5% quantile")
for(i in 1:length(Xmodulus)){
	QuantilesECDF_ModulusRing[i,1] <-  quantile(ErrorECDF[i,], prob = 0.025)
	QuantilesECDF_ModulusRing[i,2] <-  quantile(ErrorECDF[i,], prob = 0.975)
}

#Plot modulus ECDF
plot(x=Xmodulus, y=AverageECDFModulusRing, type="l", main= "", ylab="ECDF", xlab="Normalized modulus", col="red", lwd=2, cex.lab=1.5, cex.main=2, cex.axis = 1.2)
lines(x=Xmodulus, y=Xmodulus^2, col="blue", lty=3, lwd=3)
xshade <- c(Xmodulus,rev(Xmodulus))
yshade1 <- c(AverageECDFModulusRing-SEMECDFModulusRing,rev(AverageECDFModulusRing+SEMECDFModulusRing))
yshade2 <- c(QuantilesECDF_ModulusRing[,1],rev(QuantilesECDF_ModulusRing[,2]))
polygon(xshade,yshade1, border = NA, col=adjustcolor("red", 0.4))
polygon(xshade,yshade2, border = NA, col=adjustcolor("blue", 0.4))

r <- 1-mean(GeometryRing[,7])
SEMr <- sd(1-GeometryRing[,7])/sqrt(length(GeometryRing[,7]))
rect(r, 0, 1, 1000, col= adjustcolor("darkgray",alpha=0.40), border = NA)
errorbar(r, 0.2, r-SEMr, r+SEMr, 0.1)


#Average angle density
A1 <- AverageDensity(Data=MatrixRing[,c(3,2)], X=Xangle, Border=FALSE, Circular=TRUE)
AverageDensityAngleRing <- A1[[1]]
SEMDensityAngleRing <- A1[[2]]

#Error bands of modulus ECDF
ErrorDensity <- matrix(,ncol=N_sim, nrow=length(Xangle))
for(i in 1:N_sim){
	A1 <- AverageDensity(Data=CSRList[[i]][,c(3,2)], X=Xangle, Border=FALSE, Circular=TRUE)
	ErrorDensity[,i] <- A1[[1]]
}
QuantilesDensity_AngleRing <- matrix(,ncol=2,nrow=length(Xangle))
colnames(QuantilesDensity_AngleRing) <- c("2.5% quantile", "97.5% quantile")
for(i in 1:length(X)){
	QuantilesDensity_AngleRing[i,1] <-  quantile(ErrorDensity[i,], prob = 0.025)
	QuantilesDensity_AngleRing[i,2] <-  quantile(ErrorDensity[i,], prob = 0.975)
}

#Plot modulus density
plot(x=Xangle, y=AverageDensityAngleRing, type="l", main= "", ylab="Density", xlab="Normalized modulus", col="forestgreen", lwd=2, cex.lab=1.5, cex.main=2, cex.axis = 1.2, ylim=c(0,0.25), axes=FALSE)
axis(side=1, at=graduation2, labels=vec.expr, cex.axis = 1.2, lwd=4)
axis(side=2, cex.axis = 1.2, lwd=4)
lines(x=Xangle, y=rep(1/(2*pi), length=length(Xangle)), col="blue", lty=3, lwd=3)
xshade <- c(Xangle,rev(Xangle))
yshade1 <- c(AverageDensityAngleRing-SEMDensityAngleRing,rev(AverageDensityAngleRing+SEMDensityAngleRing))
yshade2 <- c(QuantilesDensity_AngleRing[,1],rev(QuantilesDensity_AngleRing[,2]))
polygon(xshade,yshade1, border = NA, col=adjustcolor("forestgreen", 0.4))
polygon(xshade,yshade2, border = NA, col=adjustcolor("blue", 0.4))

#Average angle ECDF
A1 <- AverageECDF(Data=MatrixRing[,c(3,2)], X=Xangle)
AverageECDFAngleRing <- A1[[1]]
SEMECDFAngleRing <- A1[[2]]

#Error bands of angle ECDF
ErrorECDF <- matrix(,ncol=N_sim, nrow=length(Xangle))
for(i in 1:N_sim){
	A1 <- AverageECDF(Data=CSRList[[i]][,c(3,2)], X=Xangle)
	ErrorECDF[,i] <- A1[[1]]
}
QuantilesECDF_AngleRing <- matrix(,ncol=2,nrow=length(Xangle))
colnames(QuantilesECDF_AngleRing) <- c("2.5% quantile", "97.5% quantile")
for(i in 1:length(Xangle)){
	QuantilesECDF_AngleRing[i,1] <-  quantile(ErrorECDF[i,], prob = 0.025)
	QuantilesECDF_AngleRing[i,2] <-  quantile(ErrorECDF[i,], prob = 0.975)
}

#Plot angle ECDF
plot(x=Xangle, y=AverageECDFAngleRing, type="l", main= "", ylab="ECDF", xlab="Normalized angle", col="forestgreen", lwd=2, cex.lab=1.5, cex.main=2, cex.axis = 1.2, axes=FALSE)
axis(side=1, at=graduation2, labels=vec.expr, cex.axis = 1.2, lwd=4)
axis(side=2, cex.axis = 1.2, lwd=4)
lines(x=Xangle, y=Xangle/(2*pi), col="blue", lty=3, lwd=3)
xshade <- c(Xangle,rev(Xangle))
yshade1 <- c(AverageECDFAngleRing-SEMECDFAngleRing,rev(AverageECDFAngleRing+SEMECDFAngleRing))
yshade2 <- c(QuantilesECDF_AngleRing[,1],rev(QuantilesECDF_AngleRing[,2]))
polygon(xshade,yshade1, border = NA, col=adjustcolor("forestgreen", 0.4))
polygon(xshade,yshade2, border = NA, col=adjustcolor("blue", 0.4))


##########################################################################
##9) Polarization

pvalueRing <- ResultingLengthRing <- ModulusRing <- InertiaRing <- vector()
for(i in 1:length(Exo)){
	A1 <- Exo[[i]]
	A2 <- atan2(A1[,2],A1[,1])+pi
	ResultingLengthRing[i] <- sqrt(sum(cos(A2))^2+sum(sin(A2))^2)/length(A2)
	A2 <- circular(A2)
	ModulusRing[i] <- mean(sqrt(A1[,1]^2+A1[,2]^2))
	Barycenter <- c(mean(A1[,1]), mean(A1[,2]))
	InertiaRing[i] <- mean( (A1[,1] - Barycenter[1])^2 + (A1[,2] - Barycenter[2])^2)
	pvalueRing[i] <- rayleigh.test(A2, mu=NULL)$p.value
}

hist(pvalueRing, xlim=c(0,1), breaks=seq(0,1,by=0.1), main = "Rayleigh test's p-values distribution (VAMP7Ring)", ylab="Count", xlab="p-value", col="firebrick", lwd=4, cex.lab=1.5, cex.main=2, cex.axis = 1.2)
ks.test(pvalueRing,"punif",0,1)

##########################################################################
##10) Single cell spatial analysis

Nripley <- 100

for(i in 1:length(Exo)){
	A1 <- Exo[[i]]
	p <- ppp(x=A1[,1], y=A1[,2], window=wm)
	Kfunction <- envelope(p, Kest, nsim=Nripley, nrank=1, transform = expression(. -pi*r^2), correction="best", savefuns=TRUE)
	plot(Kfunction, xlab="Distance (µm)", ylab=expression(paste("K(r)-", pi, "r²")), main=paste0("Ripley's K function (Cell", i,")"), lwd=4, cex.lab=1.5, cex.main=2, cex.axis = 1.2, legend=FALSE)
}


##########################################################################
##11) Average spatial analysis (normalized distances)

N_sim <- 100

#Windows list
wmList <- list()
for(i in 1:length(Exo)){
	wmList[[i]] <- wm
}

#Monte-Carlo CSR simulations
CSRList <- list()
for(i in 1:N_sim){
	CSRSubList <- list()
	for(j in 1:length(Exo)){
		n_exo <- nrow(Exo[[j]])
		A1 <- CSRCell(n_exo, Temporal=FALSE, Shape="Disk", R=1, Xlength=0, Ylength=0)
		CSRSubList[[j]] <- A1
	}
	CSRList[[i]] <- CSRSubList
}

#Average Ripley
A1 <- AverageRipley(DataList=Exo, wmList=wmList, Inhom=FALSE, sigma=0)
RRing <- A1[[1]]
AverageRipleyRing <- A1[[2]]
SEMRipleyRing <- A1[[3]]

#Error bands of Ripley
ErrorDensity <- matrix(,ncol=N_sim, nrow=length(RRing))
for(i in 1:N_sim){
	A1 <- AverageRipley(DataList=CSRList[[i]], wmList=wmList, Inhom=FALSE, sigma=0)
	ErrorDensity[,i] <- A1[[2]]
}
QuantilesRipleyRing <- matrix(,ncol=2,nrow=length(RRing))
colnames(QuantilesRipleyRing) <- c("2.5% quantile", "97.5% quantile")
for(i in 1:length(RRing)){
	QuantilesRipleyRing[i,1] <-  quantile(ErrorDensity[i,], prob = 0.025)
	QuantilesRipleyRing[i,2] <-  quantile(ErrorDensity[i,], prob = 0.975)
}

#Plot Ripley
plot(x=RRing, y=AverageRipleyRing - pi*RRing^2, type="l", main= "", ylab=expression(paste("K(r)-", pi, "r²")), xlab="Normalized distance", col="red", lwd=2, cex.lab=1.5, cex.main=2, cex.axis = 1.2)
lines(x=RRing, y=rep(0,length(RRing)), col="blue", lty=3, lwd=3)
xshade1 <- c(RRing,rev(RRing))
yshade1 <- c((AverageRipleyRing- pi*RRing^2)-SEMRipleyRing,rev((AverageRipleyRing- pi*RRing^2)+SEMRipleyRing))
yshade2 <- c(QuantilesRipleyRing[,1]- pi*RRing^2,rev(QuantilesRipleyRing[,2]- pi*RRing^2))
polygon(xshade1,yshade1, border = NA, col=adjustcolor("red", 0.4))
polygon(xshade1,yshade2, border = NA, col=adjustcolor("blue", 0.4))

##########################################################################
##12) Average spatial analysis (absolute distances)

N_sim <- 100

#Windows list
wmList <- list()
for(i in 1:length(Exo)){
	R <- GeometryRing[i,3]
	w <- owin(c(-R,R), c(-R,R), mask=matrix(TRUE, 1000,1000))
	X <- raster.x(w)
	Y <- raster.y(w)
	wm <- owin(w$xrange, w$yrange, mask=(X^2 + Y^2 <= R^2))
	wmList[[i]] <- wm
}

#Change coordinates
Exo2 <- list()
for(i in 1:length(Exo)){
	R <- GeometryRing[i,3]
	A1 <- Exo[[i]]
	A1[,1:2] <- A1[,1:2]*R
	Exo2[[i]] <- A1
}

#Monte-Carlo CSR simulations
CSRList <- list()
for(i in 1:N_sim){
	CSRSubList <- list()
	for(j in 1:length(Exo2)){
		n_exo <- nrow(Exo2[[j]])
		A1 <- CSRCell(n_exo, Temporal=FALSE, Shape="Disk", R=GeometryRing[j,3], Xlength=0, Ylength=0)
		CSRSubList[[j]] <- A1
	}
	CSRList[[i]] <- CSRSubList
}


#Average Ripley
A1 <- AverageRipley(DataList=Exo2, wmList=wmList, Inhom=FALSE, sigma=0)
RRing <- A1[[1]]
AverageRipleyRing <- A1[[2]]
SEMRipleyRing <- A1[[3]]


#Error bands of Ripley
ErrorDensity <- matrix(,ncol=N_sim, nrow=length(RRing))
for(i in 1:N_sim){
	A1 <- AverageRipley(DataList=CSRList[[i]], wmList=wmList, Inhom=FALSE, sigma=0)
	ErrorDensity[,i] <- A1[[2]]
}
QuantilesRipleyRing <- matrix(,ncol=2,nrow=length(RRing))
colnames(QuantilesRipleyRing) <- c("2.5% quantile", "97.5% quantile")
for(i in 1:length(RRing)){
	QuantilesRipleyRing[i,1] <-  quantile(ErrorDensity[i,], prob = 0.025)
	QuantilesRipleyRing[i,2] <-  quantile(ErrorDensity[i,], prob = 0.975)
}


#Plot Ripley
plot(x=RRing, y=AverageRipleyRing - pi*RRing^2, type="l", main= "", ylab=expression(paste("K(r)-", pi, "r²")), xlab="Distance (µm)", col="red", lwd=2, cex.lab=1.5, cex.main=2, cex.axis = 1.2, ylim=c(-20,200))
lines(x=RRing, y=rep(0,length(RRing)), col="blue", lty=3, lwd=3)
xshade1 <- c(RRing,rev(RRing))
yshade1 <- c((AverageRipleyRing- pi*RRing^2)-SEMRipleyRing,rev((AverageRipleyRing- pi*RRing^2)+SEMRipleyRing))
yshade2 <- c(QuantilesRipleyRing[,1]- pi*RRing^2,rev(QuantilesRipleyRing[,2]- pi*RRing^2))
polygon(xshade1,yshade1, border = NA, col=adjustcolor("red", 0.4))
polygon(xshade1,yshade2, border = NA, col=adjustcolor("blue", 0.4))
