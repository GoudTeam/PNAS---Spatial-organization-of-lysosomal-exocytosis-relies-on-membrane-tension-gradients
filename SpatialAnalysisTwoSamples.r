##########################################################################
##########################################################################
## Spatial analysis two paired samples
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

##########################################################################
##2) Import functions and dataset

source("Functions.R")

setwd()

Timer <- vector()	#list of movie durations in secondes
Area <- vector()	#list of cell area in µm²
ExoBefore <- list()	#list of exocytosis events coordinates (X,Y) (each element of the list is one cell) in condition 1
ExoAfter <- list()	#list of exocytosis events coordinates (X,Y) (each element of the list is one cell) in condition 2
MaskBefore <- list()	#Mask of the cell (raster transformed in matrix)
MaskAfter <- list()	#Mask of the cell (raster transformed in matrix)
#Before/After notation is used for paired data

PixelSize <- 0.160	#in µm

##########################################################################
##3) Exocytosis rate

ResultBefore <- ResultAfter <- matrix(,ncol=4, nrow=nrow(Area)/2)
colnames(ResultBefore) <- colnames(ResultAfter) <- c("Area", "Time (s)", "Number of events", "Exocytosis rate")

#Separate condition 1 and 2 (should be adjusted according to the organization of Area/Timer)
ResultBefore[,1] <- Area[seq(2, length(Area),2)]
ResultAfter[,1] <- Area[seq(1, length(Area),2)]
ResultBefore[,2] <- Timer[seq(2, length(Timer),2)]
ResultAfter[,2] <- Timer[seq(1, length(Timer),2)]

for(i in 1:nrow(ResultAfter)){
	ResultBefore[i,3] <- nrow(ExoBefore[[i]])
	ifelse(is.null(nrow(ExoAfter[[i]])), A <- 0, A <- nrow(ExoAfter[[i]]))
	ResultAfter[i,3] <- A
}

ResultBefore[,4] <- ResultBefore[,3]/(ResultBefore[,2]*ResultBefore[,1])
ResultAfter[,4] <- ResultAfter[,3]/(ResultAfter[,2]*ResultAfter[,1])

Result <- cbind(ResultBefore[,4], ResultAfter[,4])
colnames(Result) <- c("Condition 1","Condition 2")
pvalue <- wilcox.test(Result[,1], Result[,2], paired=TRUE)$p.value
Result <- as.data.frame(Result)

p <- ggpaired(Result, cond1 = "Condition 1", cond2 = "Condition 2", fill = "condition", palette = "jco",  point.size = 3, line.size = 1.5)
p <- p + ggtitle("")+theme_bw()+ylab("Exocytosis/µm²/s") +  theme(panel.border = element_blank()) 
p <- p + theme(plot.title = element_text(size=25))+theme(legend.title = element_text(size=18))+theme(legend.text = element_text(size=18))
p <- p + theme(axis.text.y = element_text(size=18))+theme(axis.title.y = element_text(size=18))+theme(axis.text.x = element_text(size=14, angle=90))+theme(axis.title.x=element_blank())
p <- p +  theme(axis.line = element_line(size=1, colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_blank(), panel.background = element_blank())
p <- p + annotate("text", x = 1.5, y = max(Result), label = paste0("p-value : ", round(pvalue, 3)), size=8)
p

##########################################################################
##3) Single cell spatial analysis


#Condition 1 (Before)
ShiftX <- rep(0, length(MaskBefore))
ShiftY <- rep(0, length(MaskBefore))
WindowListBefore <- PatternListBefore <- KfunctionListBefore <-  list()
for(i in 1:length(PathBefore)){

	Result <- ExoBefore[[i]]
	Result[,1] <- round(Result[,1])
	Result[,2] <- round(Result[,2])
	
	#Bin
	Bin <- MaskBefore[[i]]
	Bin2 <- FloodFill(Bin, Row=floor(nrow(Bin)/2+ShiftY[i]), Col=floor(ncol(Bin)/2)+ShiftX[i])
	Bin2 <- RemoveHoles(I=Bin2)
	
	#2D KDE
	w <- owin(xrange=c(1,ncol(Bin2)*PixelSize), yrange=c(0,nrow(Bin2)*PixelSize), mask=(Bin2==1))
	WindowListBefore[[i]] <- w
	PatternListBefore[[i]] <- cbind(Result[,1]*PixelSize, Result[,2]*PixelSize)
	p <- ppp(x=Result[,1]*PixelSize, y=Result[,2]*PixelSize, window=w)
	ds <- density(p)
	plot(ds, main= "",ylab="Y coordinate", xlab="X coordinate", cex.lab=1.5, cex.main=2, cex.axis = 1.2 )
	points(p, pch=20)
	
	Kfunction <- envelope(p, Kest, nsim=100, nrank=1, transform = expression(. -pi*r^2), correction="best", savefuns=TRUE)
	plot(Kfunction, xlab="Distance (µm)", ylab=expression(paste("K(r)-", pi, "r²")), main=paste0("Ripley's K function (Cell", i,")"), lwd=4, cex.lab=1.5, cex.main=2, cex.axis = 1.2, legend=FALSE )
	KfunctionListBefore[[i]] <- cbind(Kfunction$r, Kfunction$obs, Kfunction$lo, Kfunction$hi)
}


#Condition 2 (After)
ShiftX <- rep(0, length(MaskAfter))
ShiftY <- rep(0, length(MaskAfter))
WindowListAfter <- PatternListAfter <- KfunctionListAfter <-  list()
for(i in 1:length(PathAfter)){

	Result <- ExoAfter[[i]]
	Result[,1] <- round(Result[,1])
	Result[,2] <- round(Result[,2])
	
	#Bin
	Bin <- MaskAfter[[i]]
	Bin2 <- FloodFill(Bin, Row=floor(nrow(Bin)/2+ShiftY[i]), Col=floor(ncol(Bin)/2)+ShiftX[i])
	Bin2 <- RemoveHoles(I=Bin2)
	
	#2D KDE
	w <- owin(xrange=c(1,ncol(Bin2)*PixelSize), yrange=c(0,nrow(Bin2)*PixelSize), mask=(Bin2==1))
	WindowListAfter[[i]] <- w
	PatternListAfter[[i]] <- cbind(Result[,1]*PixelSize, Result[,2]*PixelSize)
	p <- ppp(x=Result[,1]*PixelSize, y=Result[,2]*PixelSize, window=w)
	ds <- density(p)
	plot(ds, main= "",ylab="Y coordinate", xlab="X coordinate", cex.lab=1.5, cex.main=2, cex.axis = 1.2 )
	points(p, pch=20)
	
	Kfunction <- envelope(p, Kest, nsim=100, nrank=1, transform = expression(. -pi*r^2), correction="best", savefuns=TRUE)
	plot(Kfunction, xlab="Distance (µm)", ylab=expression(paste("K(r)-", pi, "r²")), main=paste0("Ripley's K function (Cell", i,")"), lwd=4, cex.lab=1.5, cex.main=2, cex.axis = 1.2, legend=FALSE )
	KfunctionListAfter[[i]] <- cbind(Kfunction$r, Kfunction$obs, Kfunction$lo, Kfunction$hi)
}


##########################################################################
##3) Compair Ripley's K functions


#Before
A1 <- AverageRipley(DataList=PatternListBefore, wmList=WindowListBefore, Inhom=FALSE, sigma=0)
RBefore <- A1[[1]]
AverageRipleyBefore <- A1[[2]]
SEMRipleyBefore <- A1[[3]]

#After
A1 <- AverageRipley(DataList=PatternListAfter, wmList=WindowListAfter, Inhom=FALSE, sigma=0)
RAfter <- A1[[1]]
AverageRipleyAfter <- A1[[2]]
SEMRipleyAfter <- A1[[3]]


#Plot Ripley
Ymax <- max(c(abs(AverageRipleyBefore - pi*RBefore^2), abs(AverageRipleyAfter - pi*RAfter^2)))
plot(x=RBefore, y=AverageRipleyBefore - pi*RBefore^2, type="l", main= "", ylab=expression(paste("K(r)-", pi, "r²")), xlab="Distance (µm)", col="red", lwd=3, cex.lab=1.5, cex.main=2, cex.axis = 1.2, ylim=c(-Ymax,Ymax))
lines(x=RAfter, y=AverageRipleyAfter - pi*RAfter^2, col="goldenrod4", lwd=3)
lines(x=RBefore, y=rep(0,length(RBefore)), col="blue", lty=3, lwd=3)
xshade1 <- c(RBefore,rev(RBefore))
yshade1 <- c((AverageRipleyBefore- pi*RBefore^2)-SEMRipleyBefore,rev((AverageRipleyBefore- pi*RBefore^2)+SEMRipleyBefore))
xshade2 <- c(RAfter,rev(RAfter))
yshade2 <- c((AverageRipleyAfter- pi*RAfter^2)-SEMRipleyAfter,rev((AverageRipleyAfter- pi*RAfter^2)+SEMRipleyAfter))
polygon(xshade1,yshade1, border = NA, col=adjustcolor("red", 0.4))
polygon(xshade2,yshade2, border = NA, col=adjustcolor("goldenrod4", 0.4))
legend(0, 200, legend=c("Condition 1", "Condition 2"), col=c("red", "goldenrod4"), lty=1, cex=1, text.font=4, box.lwd=2.5, lwd=4,  bty = "n")

#Test statistique
DataBefore <- list()
for(i in 1:length(PatternListBefore)){
	w <- WindowListBefore[[i]]
	A1 <- PatternListBefore[[i]]
	DataBefore[[i]] <- ppp(x=A1[,1], y=A1[,2], window=w)
}

DataAfter <- list()
for(i in 1:length(PatternListAfter)){
	w <- WindowListAfter[[i]]
	A1 <- PatternListAfter[[i]]
	DataAfter[[i]] <- ppp(x=A1[,1], y=A1[,2], window=w)
}

X <- list()
X[[1]] <- DataBefore
X[[2]] <- DataAfter
p <- studpermu.test(X, nperm=999)
plot(p)
p