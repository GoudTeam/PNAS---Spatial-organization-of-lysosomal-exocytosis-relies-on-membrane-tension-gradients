##########################################################################
##########################################################################
## Functions
## Code from "Spatial organization of lysosomal exocytosis relies on membrane tension gradients", please cite Lachuer et al. PNAS (2023)
## Written by Hugo LACHUER
## R version 4.1.2 (2021-11-01) -- "Bird Hippie"
##########################################################################
##########################################################################

#Error bar
errorbar <- function(x, y, xmax, xmin, h, col="black", lwd=2) {
	segments(xmin, y, x1 = xmax, y1 = y, lwd=lwd, col=col)
	segments(xmin, y-h/2, x1 = xmin, y1 = y+h/2, lwd=lwd, col=col)
	segments(xmax, y-h/2, x1 = xmax, y1 = y+h/2, lwd=lwd, col=col)
}

#SEM Bar
SEMbar <- function(x, y, SEM, h, Color, lwd) {
	segments(x0 = x, y0 = y-SEM, x1 = x, y1 = y+SEM, lwd=lwd, col=Color)
	segments(x0 = x-h/2, y0 = y-SEM, x1 = x+h/2, y1 = y-SEM, lwd=lwd, col=Color)
	segments(x0 = x-h/2, y0 = y+SEM, x1 = x+h/2, y1 = y+SEM, lwd=lwd, col=Color)
}

#Generate CSR cell
CSRCell <- function(n, Temporal, DeltaT, T, Shape, R, Xlength, Ylength, Mask){
	random_cell <- matrix(, ncol=2, nrow=n)
	
	if(Shape=="Disk"){
		for(i in 1:n) {
			x <- R
			y <- R
			while((x^2 + y^2)>=R^2){
				x <- runif(1, min = -R, max = R)
				y <- runif(1, min = -R, max = R)
			}
			random_cell[i,] <- c(x,y)
		}
	}
	if(Shape=="Rectangle"){
		random_cell[,1] <- runif(n, min=0, max=Xlength)
		random_cell[,2] <- runif(n, min=0, max=Ylength)
	}
	if(Shape=="Free"){
		if(is.owin(Mask)){
			Mask <- Mask[[10]]*1
		}
		for(i in 1:n){
			First <- TRUE
			x <- y <- 1
			while(Mask[y,x]==0 | First == TRUE){
				First <- FALSE
				x <- sample(1:ncol(Mask),1)
				y <- sample(1:nrow(Mask),1)
			}
			random_cell[i,] <- c(x,y)
		}
	}
	
	if(Temporal){
		TimeLine <- seq(DeltaT, T, DeltaT)
		Time <- sample(TimeLine, n, replace=TRUE)
		random_cell <- cbind(random_cell, Time)
	}
	
	return(random_cell)
}

#Measure average value of the map (matrix) at positions X and Y
MeasureMap <- function(Map, X, Y){
	A <-  0
	for(i in 1:length(X)){
		A <- A + Map[Y[i],X[i]]
	}
	A <- A/length(X)
	return(A)
}

##Compute all NND
NND <- function(X,Y){
	Output <- vector(, length=length(X))
	for(i in 1:length(X)){
		Index <- 1:length(X)
		Index <- Index[-i]
		Output[i] <- min(sqrt((X[i]-X[Index])^2+(Y[i]-Y[Index])^2))
	}
	return(Output)
}

##NND test
NNDtest <- function(X, Y, N_CSR_Sim, Shape, R, Xlength, Ylength, Mask){
	
	NND_Observed <- NND(X=X,Y=Y)
    B <- mean(NND_Observed)

    number_exo <- length(X)
	NND_MC <- vector(, length = N_CSR_Sim)
    for(k in 1:N_CSR_Sim) {
		random_cell <- CSRCell(n=number_exo, Temporal=FALSE, Shape=Shape, R=R, Xlength=Xlength, Ylength=Ylength, Mask=Mask)
		NND_simulated <- NND(X=random_cell[,1], Y=random_cell[,2])
		NND_MC[k] <- mean(NND_simulated)
    }

    NLess <- length(which(NND_MC < B))
    NGreater <- length(which(NND_MC > B))
	pvalueLess <- (NLess+1)/(length(NND_MC)+1)
    pvalueGreater <- (NGreater+1)/(length(NND_MC)+1)	#Empirical p-value
    pvalue <- min(pvalueLess,pvalueGreater)*2	#two-sided
	if(pvalue>1){pvalue<-1}
    ifelse(B>mean(NND_MC), Clustering <- FALSE, Clustering <- TRUE)
    d1 <- density(NND_MC, bw="SJ", kernel = "gaussian")
	
	Output <- list(pvalue, Clustering, NND_MC, B)
	return(Output)
}

#Average 2D KDE
Average2DKDE <- function(CoordList, wm, sigma=0){
	
	HeatMap <- array(,dim=c(nrow(wm),ncol(wm),length(CoordList)))
	for(i in 1:length(CoordList)){
		A1 <- CoordList[[i]]
		if(sigma>0){
			A1 <- A1 + cbind(rnorm(nrow(A1), mean=0, sd= sigma), rnorm(nrow(A1), mean=0, sd= sigma))	#sigma add a small noise to avoid duplicates
		}
			
		#2D KDE
		p <- ppp(x=A1[,1], y=A1[,2], window=wm)
		ds <- density(p, kernel="gaussian", sigma=bw.diggle, edge=TRUE, diggle=TRUE)
		HeatMap[,,i] <- ds[[1]]/sum(ds[[1]], na.rm=TRUE)	#Normalized total to 1
	}
	
	AverageHeatMap <- SEMHeatMap <- MedianHeatMap <- matrix(,ncol=ncol(wm), nrow=nrow(wm))
	for(i in 1:nrow(wm)){
		for(j in 1:ncol(wm)){
			AverageHeatMap[i,j] <- mean(HeatMap[i,j,], na.rm=TRUE)
			SEMHeatMap[i,j] <- sd(HeatMap[i,j,], na.rm=TRUE)/sqrt(length(CoordList))
			MedianHeatMap[i,j] <- median(HeatMap[i,j,], na.rm=TRUE)
		}
	}
	rm(HeatMap)
	
	Output <- list()
	Output[[1]] <- AverageHeatMap
	Output[[2]] <- SEMHeatMap
	Output[[3]] <- MedianHeatMap
	return(Output)
}

#ECDF
ECDF <- function(x,Data){
	y <- length(which(Data<=x))/length(Data)
	return(y)
}

#Average ECDF
AverageECDF <- function(Data, X){
	Data[,1] <- as.factor(Data[,1])
	A1 <- matrix(,ncol=length(X), nrow=length(levels(Data[,1])))
	for(i in 1:nrow(A1)){
		V1 <- Data[which(Data[,1] == levels(Data[,1])[i]),2]
		for(j in 1:length(X)){
			A1[i,j] <- ECDF(x=X[j], Data=V1)
		}
	}
	
	Average <- SEM <- vector(, length=length(X))
	for(i in 1:length(X)){
		Average[i] <- mean(A1[,i])
		SEM[i] <- sd(A1[,i])/sqrt(nrow(A1))
	}
	A2 <- list()
	A2[[1]] <- Average
	A2[[2]] <- SEM
	return(A2)
}

#Average densities
AverageDensity <- function(Data, X, Border, Circular, xmax){
	
	if(Border){
		if(length(xmax)>1){
			Name <- unique(Data[,1])
			names(xmax) <- Name
		}
	}

	Data[,1] <- as.factor(Data[,1])
	A1 <- matrix(,ncol=length(X), nrow=length(levels(Data[,1])))
	for(i in 1:nrow(A1)){
		V1 <- Data[which(Data[,1] == levels(Data[,1])[i]),2]
		
		if(Border){
			if(length(xmax)>1){
				A1[i,] <- dbckden(X,V1, lambda = 0.02, kernel = "gaussian", bcmethod = "beta1", proper = TRUE, nn = "jf96", offset = NULL, xmax = xmax[which(names(xmax)==levels(Data[,1])[i])], log = FALSE)
			} else {
				A1[i,] <- dbckden(X,V1, lambda = 0.02, kernel = "gaussian", bcmethod = "beta1", proper = TRUE, nn = "jf96", offset = NULL, xmax = xmax, log = FALSE)
			}
		}
		
		if(Circular){
			V1 <- as.circular(V1, type="angles", units="radian", template="none", modulo="asis", zero=0, rotation="counter")
			bw <- bw.nrd.circular(V1)
			d1 <- density.circular(V1, bw=bw, from=circular(0), to=circular(2*pi), kernel = "wrappednormal")
			A1[i,] <- approx(d1$x, d1$y, xout = X)$y
		}
		
		if(Circular==FALSE & Border==FALSE){
			d1 <- density(x=V1, bw="SJ", kernel="gaussian")
			if(integrate.xy(d1$x, d1$y)>1.1){	#Sometimes SJ bandwidth failed
				d1 <- density(x=V1, kernel="gaussian")
			}
			A1[i,] <- approx(d1$x, d1$y, xout = X)$y
		}
	}
	Average <- SEM <- vector(, length=length(X))
	for(i in 1:length(X)){
		Average[i] <- mean(A1[,i], na.rm=TRUE)
		SEM[i] <- sd(A1[,i], na.rm=TRUE)/sqrt(length(which(is.na(A1[,i])==FALSE)))
	}
	A2 <- list()
	A2[[1]] <- Average
	A2[[2]] <- SEM
	return(A2)
}

#Average Ripley
AverageRipley <- function(DataList, wmList, Inhom, sigma){

	KfunctionList <- list()
	Index <- which(sapply(DataList, is.null))
	if(length(Index)>0){
		DataList <- DataList[-Index]
		wmList <- wmList[-Index]
	}
	
	
	for(i in 1:length(DataList)){
		A1 <- DataList[[i]]
		p <- ppp(x=A1[,1]+rnorm(nrow(A1),0,sd=sigma), y=A1[,2]+rnorm(nrow(A1),0,sd=sigma), window=wmList[[i]])	#sigma add a small noise to avoid duplicates
		
		if(Inhom){
			A2 <- Kinhom(p, correction="best", sigma=bw.diggle)
			KfunctionList[[i]] <- cbind(A2$r, A2$trans)
		}else{
			A2 <- Kest(p, correction="best")
			KfunctionList[[i]] <- cbind(A2$r, A2$trans)
		}
	}
	

	Rmax <- vector(, length=length(KfunctionList))
	for(i in 1:length(Rmax)){
		Rmax[i] <- max(KfunctionList[[i]][,1])
	}
	Rindex <- seq(0, median(Rmax), median(Rmax)/1000)
	

	RipleyMatrix <- matrix(,ncol=length(KfunctionList), nrow=length(Rindex))
	for(i in 1:length(KfunctionList)){
		RipleyMatrix[,i] <- approx(x=KfunctionList[[i]][,1], y=KfunctionList[[i]][,2], xout=Rindex, method="linear")$y
	}
	
	RipleyAverage  <- RipleySEM <- vector(, length=length(Rindex))
	for(i in 1:length(Rindex)){
		RipleyAverage[i] <- mean(RipleyMatrix[i,], na.rm=TRUE)
		RipleySEM[i] <- sd(RipleyMatrix[i,],na.rm=TRUE)/sqrt(length(which(is.na(RipleyMatrix[i,])==FALSE)))
	}
	
	Output <- list()
	Output[[1]] <- Rindex
	Output[[2]] <- RipleyAverage
	Output[[3]] <- RipleySEM
	return(Output)
}


#AverageCurve
AverageCurve <- function(DataList){

	Index <- which(sapply(DataList, is.null))
	if(length(Index)>0){
		DataList <- DataList[-Index]
	}

	Rmax <- vector(, length=length(DataList))
	for(i in 1:length(Rmax)){
		Rmax[i] <- max(DataList[[i]][,1])
	}
	Rmin <- vector(, length=length(DataList))
	for(i in 1:length(Rmin)){
		Rmin[i] <- min(DataList[[i]][,1])
	}
	Rindex <- seq(median(Rmin), median(Rmax), (median(Rmax)-median(Rmin))/1000)
	
	DataMatrix <- matrix(,ncol=length(DataList), nrow=length(Rindex))
	for(i in 1:length(DataList)){
		DataMatrix[,i] <- approx(x=DataList[[i]][,1], y=DataList[[i]][,2], xout=Rindex, method="linear")$y
	}
	
	DataAverage  <- DataSEM <- vector(, length=length(Rindex))
	for(i in 1:length(Rindex)){
		DataAverage[i] <- mean(DataMatrix[i,], na.rm=TRUE)
		DataSEM[i] <- sd(DataMatrix[i,],na.rm=TRUE)/sqrt(length(which(is.na(DataMatrix[i,])==FALSE)))
	}
	
	#Output
	Output <- list()
	Output[[1]] <- Rindex
	Output[[2]] <- DataAverage
	Output[[3]] <- DataSEM
	return(Output)
}


#1D Ripley's K function (Ripley's borders correction)
Ripley1D <- function(Points, dL, L, Centered, Standardized){
	Ntime <- seq(2*dL,L/4, by=dL*2)
	Q <- (length(Points)-1)/L
	KVector <- vector(, length=length(Ntime)+1)
	Nb <- 2
	KVector[1] <- 0
	for(i in Ntime){
		K <- 0
		for(j in 1:length(Points)){
			Index <- 1:length(Points)
			Index <- Index[-j]
			for(k in Index){
				if(abs(Points[j]-Points[k]) <= i){
					d <- abs(Points[j]-Points[k])
					d1 <- min(Points[j], L-Points[j])
					d2 <- min(Points[k], L-Points[k])
					ifelse(d>d1, k1 <- 2, k1 <- 1)
					ifelse(d>d2, k2 <- 2, k2 <- 1)
					K <- K + 1/2*(k1+k2)
				}
			}
		}
		KVector[Nb] <- K
		Nb <- Nb + 1
	}
	KVector <- KVector/(length(Points)*Q)	#averaging and normalization
	
	if(Centered == TRUE){
		Nb <- 1
		for(i in c(0,Ntime)){
			KVector[Nb] <- KVector[Nb] - 2*i
			Nb <- Nb+1
		}
	}
	
	if(Standardized == TRUE){
		Nb <- 2
		for(i in Ntime[-1]){
			KVector[Nb] <- KVector[Nb]/sqrt(RipleyVar(length(Points), i, L))
			Nb <- Nb+1
		}
	}
	Output <- cbind(c(0,Ntime),KVector)
	return(Output)
}

#Average 1D Ripley's K function
AverageTemporalRipley <- function(DataList, DeltaTVector, TmaxVector, Centered, Standardized){

	KfunctionList <- list()
	for(i in 1:length(DataList)){
		A1 <- DataList[[i]]
		KfunctionList[[i]] <- Ripley1D(Points=DataList[[i]], dL=DeltaTVector[i], L=TmaxVector[i], Centered=Centered, Standardized=Standardized)
	}
	
	Rmax <- vector(, length=length(KfunctionList))
	for(i in 1:length(Rmax)){
		Rmax[i] <- max(KfunctionList[[i]][,1])
	}
	Rindex <- seq(0, median(Rmax), median(Rmax)/1000)
	
	RipleyMatrix <- matrix(,ncol=length(KfunctionList), nrow=length(Rindex))
	for(i in 1:length(KfunctionList)){
		RipleyMatrix[,i] <- approx(x=KfunctionList[[i]][,1], y=KfunctionList[[i]][,2], xout=Rindex, method="linear")$y
	}
	
	RipleyAverage  <- RipleySEM <- vector(, length=length(Rindex))
	for(i in 1:length(Rindex)){
		RipleyAverage[i] <- mean(RipleyMatrix[i,], na.rm=TRUE)
		RipleySEM[i] <- sd(RipleyMatrix[i,],na.rm=TRUE)/sqrt(length(which(is.na(RipleyMatrix[i,])==FALSE)))
	}
	
	#Output
	Output <- list()
	Output[[1]] <- Rindex
	Output[[2]] <- RipleyAverage
	Output[[3]] <- RipleySEM
	return(Output)
}

#Median spatio-temporal Ripley's K function
MedianSpatioTemporalRipley <- function(ExoList, dT, T, RasterList, dX){

	#remove cells with less than 10 events
	for(i in length(ExoList):1){
		if(nrow(ExoList[[i]])<10){
			ExoList[[i]] <- RasterList[[i]] <- NULL
			dT <- dT[-i]
			T <- T[i]
			dX <- dX[i]
		}
	}

	#Boundaries
	Tmax <- median(T)
	TypicalLenghts <- vector()
	for(i in 1:length(ExoList)){
		Mask <- RasterList[[i]]
		Surface <- sum(Mask)*dX[i]^2
		TypicalLenghts[i] <- sqrt(Surface/pi)
	}
	TypicalLenght <- median(TypicalLenghts)
	
	KST <- array(dim=c(length(seq(0,TypicalLenght,0.1)),length(seq(0,Tmax/4,Tmax/200)),length(ExoList)))
	for(i in 1:length(ExoList)){
		Mask <- RasterList[[i]]
		Events <- ExoList[[i]]
		Tmax <- T[i]
		dT2 <- dT[i]
		dX2 <- dX[i]
		

		SpacePoints <- as.points(Events[,1]*dX2, Events[,2]*dX2)
		SpacePolygon <- MooreNeighborTracing(Mask)*dX2
		TimePoints <- Events[,3]*dT2
	
		A1 <- stkhat(pts=SpacePoints, times=TimePoints, poly=SpacePolygon, tlimits=c(0, Tmax), s=seq(0,TypicalLenght,0.1), tm=seq(0,Tmax/4,Tmax/200))
		Ks_Kt <- A1$ks%*%t(A1$kt)
		KST[,,i] <- A1$kst - Ks_Kt	#test independency
	}
	
	MuKST <- matrix(,nrow=dim(KST)[1],ncol=dim(KST)[2])
	for(i in 1:dim(KST)[1]){
		for(j in 1:dim(KST)[2]){
			MuKST[i,j] <- median(KST[i,j,])
		}
	}
	
	Output <- list()
	Output[[1]] <- MuKST
	Output[[2]] <- seq(0,TypicalLenght,0.1)
	Output[[3]] <- seq(0,Tmax/4,Tmax/200)
	return(Output)

}

#Fourier Analysis
FourierAnalysis <- function(Temporal_vector, Fs, Shuffling){

	if(Shuffling == TRUE){
		Temporal_vector <- sample(Temporal_vector)
	}

	N <- length(Temporal_vector)
	Fn <- Fs/2	#Nyquist frquency, note that the FFT resolution is Fs/N
	critic_index <- floor(Fn*N/Fs)
	Temporal_vector <- Temporal_vector - mean(Temporal_vector)
	FFTComplexe <- fft(Temporal_vector)
	
	MatrixResults <- matrix(, ncol=2, nrow=N)
	colnames(MatrixResults) <- c("Fréquence (Hz)", "Module")
	MatrixResults[,1] <- (0:(N-1))*Fs/N
	MatrixResults[,2] <- Mod(FFTComplexe)
	
	Spectrum <- plot(y=MatrixResults[1:critic_index,2], x=MatrixResults[1:critic_index,1]*1000,type="l", main = "Modulus of Fast Fourier Transform ", ylab = "Modulus", xlab ="Frequency (mHz)", cex.lab=1.5, cex.main=2, lwd=4, col = "orange")
	
	
	k <- which.max(MatrixResults[1:critic_index,2])
	Frequency <- paste0("fréquence propre :", MatrixResults[k,1]*1000, " mHz")
	
	return(list(Spectrum, MatrixResults[k,1], MatrixResults))
}

#Average Fourier
AverageFourier <- function(TemporalDataList, DeltaTVector, Shuffling, Normalization){

	FourierList <- list()
	Frequencies <- vector()
	for(i in 1:length(TemporalDataList)){
		A1 <- FourierAnalysis(Temporal_vector = TemporalDataList[[i]], Fs=1/DeltaTVector[[i]], Shuffling=Shuffling)
		Frequencies[i] <- A1[[2]]
		if(Normalization==TRUE){
			A1[[3]] <- A1[[3]]/max(A1[[3]][,2])
		}
		FourierList[[i]] <- A1[[3]]
	}


	fmax <- vector()
	for(i in 1:length(FourierList)){
		fmax[i] <- max(FourierList[[i]][,1])
	}
	findex <- seq(0, median(fmax), median(fmax)/1000)
	

	FourierMatrix <- matrix(,ncol=length(FourierList), nrow=length(findex))
	for(i in 1:length(FourierList)){
		FourierMatrix[,i] <- approx(x=FourierList[[i]][,1], y=FourierList[[i]][,2], xout=findex, method="linear")$y
	}
	
	FourierAverage  <- FourierSEM <- vector(, length=length(findex))
	for(i in 1:length(findex)){
		FourierAverage[i] <- mean(FourierMatrix[i,], na.rm=TRUE)
		FourierSEM[i] <- sd(FourierMatrix[i,],na.rm=TRUE)/sqrt(length(which(is.na(FourierMatrix[i,])==FALSE)))
	}
	
	#Output
	Output <- list()
	Output[[1]] <- findex
	Output[[2]] <- FourierAverage
	Output[[3]] <- FourierSEM
	Output[[4]] <- Frequencies
	return(Output)
}


#FloodFill algorithm
FloodFill <- function(I,Row,Col){

	I <- rbind(rep(0,length=ncol(I)), I, rep(0,length=ncol(I)))
	I <- cbind(rep(0,length=nrow(I)), I, rep(0,length=nrow(I)))
	Row <- Row + 1
	Col <- Col + 1

	if(I[Row,Col] == 1){
		P <- matrix(,ncol=2)
		P[1,] <- c(Row,Col)
		I[Row,Col] <- 2
	} else {
		A1 <- which(I==1, arr.ind=TRUE)
		Row <- A1[1,1]
		Col <- A1[1,2]
		P <- matrix(,ncol=2)
		P[1,] <- c(Row,Col)
		I[Row,Col] <- 2
	}
	
	Pfinal <- P
	
	while(nrow(P) > 0){
	
		P1 <- P
		P <- matrix(,nrow=0,ncol=2)
		
		#4-connexe
		for(i in 1:nrow(P1)){
			if(I[P1[i,1]+1,P1[i,2]] == 1){
				I[P1[i,1]+1,P1[i,2]] <- 2
				P <- rbind(P,c(P1[i,1]+1,P1[i,2]))
			}
			if(I[P1[i,1]-1,P1[i,2]] == 1){
				I[P1[i,1]-1,P1[i,2]] <- 2
				P <- rbind(P,c(P1[i,1]-1,P1[i,2]))
			}
			if(I[P1[i,1],P1[i,2]+1] == 1){
				I[P1[i,1],P1[i,2]+1] <- 2
				P <- rbind(P,c(P1[i,1],P1[i,2]+1))
			}
			if(I[P1[i,1],P1[i,2]-1] == 1){
				I[P1[i,1],P1[i,2]-1] <- 2
				P <- rbind(P,c(P1[i,1],P1[i,2]-1))
			}
		}
		if(nrow(P) > 0)	{Pfinal <- rbind(Pfinal,P)}
	}

	I2 <- matrix(0, ncol=ncol(I), nrow=nrow(I))
	for(i in 1:nrow(Pfinal)){
		I2[Pfinal[i,1], Pfinal[i,2]] <- 1
	}
	
	I2 <- I2[-nrow(I2),]
	I2 <- I2[-1,]
	I2 <- I2[,-ncol(I2)]
	I2 <- I2[,-1]
	
	return(I2)
}

#Remove holes in a mask
RemoveHoles <- function(I){
	
	I <- rbind(rep(0,length=ncol(I)), I, rep(0,length=ncol(I)))
	I <- cbind(rep(0,length=nrow(I)), I, rep(0,length=nrow(I)))

	I <- rbind(rep(1,length=ncol(I)), I, rep(1,length=ncol(I)))
	I <- cbind(rep(1,length=nrow(I)), I, rep(1,length=nrow(I)))
	
	Row <- 2
	Col <- 2
	
	#FloodFill
	P <- matrix(,ncol=2)
	P[1,] <- c(Row,Col)
	I[Row,Col] <- 2
	Pfinal <- P
	
	while(nrow(P) > 0){
	
		P1 <- P
		P <- matrix(,nrow=0,ncol=2)
		
		#4-connexe
		for(i in 1:nrow(P1)){
			if(I[P1[i,1]+1,P1[i,2]] == 0){
				I[P1[i,1]+1,P1[i,2]] <- 2
				P <- rbind(P,c(P1[i,1]+1,P1[i,2]))
			}
			if(I[P1[i,1]-1,P1[i,2]] == 0){
				I[P1[i,1]-1,P1[i,2]] <- 2
				P <- rbind(P,c(P1[i,1]-1,P1[i,2]))
			}
			if(I[P1[i,1],P1[i,2]+1] == 0){
				I[P1[i,1],P1[i,2]+1] <- 2
				P <- rbind(P,c(P1[i,1],P1[i,2]+1))
			}
			if(I[P1[i,1],P1[i,2]-1] == 0){
				I[P1[i,1],P1[i,2]-1] <- 2
				P <- rbind(P,c(P1[i,1],P1[i,2]-1))
			}
		}
		if(nrow(P) > 0)	{Pfinal <- rbind(Pfinal,P)}
	}
	
	I2 <- matrix(1, ncol=ncol(I), nrow=nrow(I))
	for(i in 1:nrow(Pfinal)){
		I2[Pfinal[i,1], Pfinal[i,2]] <- 0
	}

	I2 <- I2[-c(nrow(I2)-1, nrow(I2)),]
	I2 <- I2[-c(1,2),]
	I2 <- I2[,-c(ncol(I2)-1, ncol(I2))]
	I2 <- I2[,-c(1,2)]
	
	return(I2)
}

#Binning
Binning <- function(Image, Bin){
	Image2 <- matrix(,ncol=floor(ncol(Image)/Bin), nrow=floor(nrow(Image)/Bin))
	ki <- 1
	for(i in seq(1, floor(ncol(Image)/Bin)*Bin, Bin)){
		kj <- 1
		for(j in seq(1, floor(nrow(Image)/Bin)*Bin, Bin)){
			Image2[kj,ki] <- mean(Image[j:(j+Bin-1),i:(i+Bin-1)], na.rm=TRUE)
			kj <- kj + 1
		}
		ki <- ki + 1
	}
	return(Image2)
}

#Center of Mass
CenterMass <- function(Image){
	Xmatrix <- matrix(rep(1:ncol(Image), nrow(Image)), ncol=ncol(Image), nrow=nrow(Image), byrow=TRUE)
	Ymatrix <- matrix(rep(1:nrow(Image), ncol(Image)), ncol=ncol(Image), nrow=nrow(Image), byrow=FALSE)
	Xmass <- sum(Xmatrix * Image)/sum(Image)
	Ymass <- sum(Ymatrix * Image)/sum(Image)
	return(c(Xmass, Ymass))
}


#Moore's neighbors
MooreNeighbors <- function(x,y){
	Moore <- matrix(, ncol=2, nrow=8)
	colnames(Moore) <- c("x", "y")
	Moore[1:3,1] <- (x-1):(x+1)
	Moore[1:3,2] <- rep(y-1, 3)
	Moore[4,] <- c(x+1,y)
	Moore[5:7,1] <- (x+1):(x-1)
	Moore[5:7,2] <- rep(y+1, 3)
	Moore[8,] <- c(x-1,y)
	return(Moore)
}

#Moore Nieghbor Boundary tracing algorithm
MooreNeighborTracing <- function(Bin){

	IndexMatrix <- matrix(c(1,2,3,8,0,4,7,6,5), ncol=3, byrow=TRUE)

	Bin <- rbind(rep(0,length=ncol(Bin)), Bin, rep(0,length=ncol(Bin)))
	Bin <- cbind(rep(0,length=nrow(Bin)), Bin, rep(0,length=nrow(Bin)))

	Bin <- Bin/max(Bin)

	#Found the starting point (right to left, bottom to top)
	icol <- 2
	while(sum(Bin[,icol])==0){
		icol <- icol+1
	}
	irow <- nrow(Bin)
	while(Bin[irow, icol]==0){
		irow <- irow-1
	}
	
	#First point
	Boundaries <- matrix(,ncol=2,nrow=1)
	colnames(Boundaries) <- c("x", "y")
	Boundaries[1,] <- c(icol,irow)
	Pointer <- c(Boundaries[1,1], Boundaries[1,2]+1)	#Current location
	Diss <- Pointer - Boundaries[1,]
	Index <- IndexMatrix[Diss[2]+2,Diss[1]+2]
	Neighboors <- MooreNeighbors(Boundaries[1,1],Boundaries[1,2])
	Neighboors <- rbind(Neighboors,Neighboors)
	Neighboors <- Neighboors[Index:(Index+7),]
	k <- 2
	while(Bin[Neighboors[k,2], Neighboors[k,1]]==0){
		k <- k + 1
	}
	Boundaries<- rbind(Boundaries,c(Neighboors[k,1], Neighboors[k,2]))
	Pointer <- c(Neighboors[k-1,1], Neighboors[k-1,2])
	
	#Iteration
	j <- 2
	while((Boundaries[j,1] != Boundaries[1,1]) | (Boundaries[j,2] != Boundaries[1,2])){
		Diss <- Pointer - Boundaries[j,]
		Index <- IndexMatrix[Diss[2]+2,Diss[1]+2]
		Neighboors <- MooreNeighbors(Boundaries[j,1],Boundaries[j,2])
		Neighboors <- rbind(Neighboors,Neighboors)
		Neighboors <- Neighboors[Index:(Index+7),]
		k <- 2
		while(Bin[Neighboors[k,2], Neighboors[k,1]]==0){
			k <- k + 1
		}
		Boundaries<- rbind(Boundaries,c(Neighboors[k,1], Neighboors[k,2]))
		Pointer <- c(Neighboors[k-1,1], Neighboors[k-1,2])
		j <- j + 1
	}
	
	#Rescale
	Boundaries <- Boundaries - 1
	
	return(Boundaries)
}

DistanceTransform <- function(Image, MaskSize=3){
	Bin <- IsoData(Image, N=1000)
	
	#Bin2 <- Closing(Bin)	#remove small holes
	
	#Remove holes
	Bin2 <- FloodFill(Bin, Row=floor(nrow(Bin)/2), Col=floor(ncol(Bin)/2))
	Bin2 <- RemoveHoles(I=Bin2)
	
	Bin3 <- distance_transform(as.cimg(Bin2), value=0, metric=2)
	Bin3 <- as.matrix(Bin3)
	return(Bin3)
}


#Compute Moran Indices
MoranIndex <- function(DataList, Bin, Sampling, Nmax){

	Moran <- list()

	for(i in 1:length(DataList)){
		Map <- DataList[[i]]
		
		if(Bin>1){
			Map <- Binning(Image=Map, Bin=Bin)
		}
		

		DF <- as.data.frame(matrix(,ncol=3,nrow=ncol(Map)*nrow(Map)))
		colnames(DF) <- c("x", "y", "value")

		DF[,1] <- rep(1:ncol(Map), each=nrow(Map))
		DF[,2] <- rep(1:nrow(Map),ncol(Map))
		DF[,3] <- as.vector(Map)
		Indices <- which(is.na(DF[,3]))
		DF <- DF[-Indices,]
		
		if(Sampling==TRUE){
			if(nrow(DF)>Nmax){
				DF <- DF[sample(1:nrow(DF), Nmax),]	#random sampling
			}
		}
		
		
		Data <- DF[,3]
		WeightMatrix <- matrix(0,ncol=Nmax, nrow=Nmax)
		for(j in 1:Nmax){
			for(k in j:Nmax){
				WeightMatrix[j,k] <- 1/sqrt((DF[j,1] - DF[k,1])^2 + (DF[j,2] - DF[k,2])^2)
			}
		}
		
		diag(WeightMatrix) <- 0	#0 diagonal
		WeightMatrix[lower.tri(WeightMatrix)] = t(WeightMatrix)[lower.tri(WeightMatrix)]	#Symetry
		WeightMatrix <- WeightMatrix/sum(WeightMatrix, na.rm=TRUE)	#Normalized
		
		Moran[[i]] <- Moran.I(Data, WeightMatrix, na.rm=TRUE)
		
		print(paste0("Progression: ",round(i/length(DataList)*100,2),"%"))
	}
	rm(WeightMatrix)
	rm(Data)
	return(Moran)
}

#Compute polarization index
ComputPolarIndex <- function(value,angle){
	P <- value/sum(value, na.rm=TRUE)
	C <- sum(P*cos(angle), na.rm=TRUE)
	S <- sum(P*sin(angle), na.rm=TRUE)
	Vobs <- sqrt(C^2+S^2)	#Polarization index
	return(Vobs)
}

#Polarization of exocytosis
TestPolarExocytosis <- function(Bin, Coordinates, NPart, Nsim){

	#if only one event
	if(is.vector(Coordinates)){
		Coordinates <- matrix(Coordinates,ncol=2,nrow=1)
	}
	
	#Required angles
	Angles <- seq(from=0,to=2*pi, length.out=Npart+1)
	Angles2 <- Angles + 2*pi/(2*Npart)
	Angles2 <- Angles2[-(Npart+1)]
		
	#Generate radial coordinates of the Bin
	Center <- CenterMass(Bin)
	Indices <- which(Bin==1, arr.ind=TRUE)
	RadialCell <- vector()
	for(j in 1:nrow(Indices)){
		theta <- atan2(Indices[j,1]-Center[2], Indices[j,2]-Center[1]) + pi 	#Angle between 0 and 2pi
		RadialCell[j] <- theta
	}
	
	#Compute surfaces associated to each angles interval
	Surface <- vector(,length=Npart)
	for(j in 1:Npart){
		Surface[j] <- length(which((RadialCell > Angles[j]) & (RadialCell < Angles[j+1])))
	}
		
	#Generate radial coordinates of exocytosis events
	RadialExo <- vector()
	for(j in 1:nrow(Coordinates)){
		theta <- atan2(Coordinates[j,2]-Center[2], Coordinates[j,1]-Center[1]) + pi
		RadialExo[j] <- theta
	}
		
	#Compute observed polarization index
	CellDistribution <- vector(,length=Npart)
	for(j in 1:Npart){
		Nexo <- length(which((RadialExo > Angles[j]) & (RadialExo < Angles[j+1])))
		CellDistribution[j] <- Nexo/Surface[j]
	}
	Vobs <- ComputPolarIndex(value=CellDistribution, angle=Angles2)
	print(paste0("Observed polarization index: ", round(Vobs,4)))
		
	#Compute simulations
	Vsim <- vector()
	for(i in 1:Nsim){
		
		ResultCSR <- CSRCell(n=nrow(Coordinates), Temporal=FALSE, Shape="Free", Mask=Bin2)
			
		#Generate radial coordinates of exocytosis events
		RadialExo <- vector()
		for(j in 1:nrow(ResultCSR)){
			theta <- atan2(ResultCSR[j,2]-Center[2], ResultCSR[j,1]-Center[1]) + pi
			RadialExo[j] <- theta
		}
		#Compute observed polarization index
		CellDistribution <- vector(,length=Npart)
		for(j in 1:Npart){
			Nexo <- length(which((RadialExo > Angles[j]) & (RadialExo < Angles[j+1])))
			CellDistribution[j] <- Nexo/Surface[j]
		}
		Vsim[i] <- ComputPolarIndex(value=CellDistribution, angle=Angles2)
	}
	print(paste0("Average simulated polarization index: ", round(mean(Vsim),4)))
		
	#compute the pvalue
	Nless <- length(which(Vsim<Vobs))
	Nmore <- length(which(Vsim>Vobs))
	pvalue <- min(c(Nless/Nsim,Nmore/Nsim))*2
	return(c(pvalue,mean(Vsim)))
}