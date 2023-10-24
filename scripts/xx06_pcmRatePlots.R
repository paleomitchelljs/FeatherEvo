setwd("~/Dropbox/Research/FeatherEvolution/scripts")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\scripts")
source("xx00_readFormatDat.R")

# this calculates the weight of the ONE RATE model
calcWt <- function(x) 	{
	if (x > 0)	{ # then 2-rate is preferred model
		Ex <- exp(-0.5*x)
		Wt <- Ex / (Ex + 1)
	}
	if (x < 0)	{
		Ex <- exp(0.5*x)
		Wt <- Ex / (Ex + 1)
		Wt <- 1 - Wt
	}
	return(Wt)
}

setwd("~/Dropbox/Research/FeatherEvolution/output")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\output")
flightless <- read.csv("lessMat.csv", row.names=1)
flying <- read.csv("volMat.csv", row.names=1)
oneRate <- read.csv("oneMat.csv", row.names=1)
AICcs <- read.csv("aiccMat.csv", row.names=1)
oAICcs <- read.csv("oaiccMat.csv", row.names = 1)

dAICcs <- oAICcs - AICcs
wtMat <- apply(dAICcs, 2, function(x) sapply(x, calcWt))

#colnames(flightless) <- colnames(ModelDat)[7:ncol(ModelDat)]
#colnames(flying) <- colnames(ModelDat)[7:ncol(ModelDat)]
#colnames(wtMat) <- colnames(ModelDat)[7:ncol(ModelDat)]

colnames(flightless) <- c("mass", "wingR", "tailR", "tarsusR", "rachisWL", "barbDR", "barbLR", "barbAng", "barbL_lead")
colnames(flying) <- c("mass", "wingR", "tailR", "tarsusR", "rachisWL", "barbDR", "barbLR", "barbAng", "barbL_lead")
colnames(wtMat) <- c("mass", "wingR", "tailR", "tarsusR", "rachisWL", "barbDR", "barbLR", "barbAng", "barbL_lead")


makePlot <- function(mat1, mat2, colu, colOpt=1, Alpha=0.5)	{
	if (colOpt == 1)	{
		col1 <- rgb(179/255,88/255,6/255, Alpha)
		b1 <- rgb(179/255,88/255,6/255, 1)

		col2 <- rgb(84/255,39/255,136/255, Alpha)
		b2 <- rgb(84/255,39/255,136/255, 1)

	}
	else {
		col1 <- rgb(84/255,39/255,136/255, Alpha)
		b1 <- rgb(84/255,39/255,136/255, 1)

		col2 <- rgb(179/255,88/255,6/255, Alpha)
		b2 <- rgb(179/255,88/255,6/255, 1)

	}
	xRange <- range(c(mat1[,colu], mat2[,colu]))
	den1 <- hist(mat1[,colu], breaks=seq(from=0, to=max(xRange), by=max(xRange)/8), ylim=c(0,100), col=col1, border=b1, xlab="", ylab="", main=colnames(mat1)[colu])
	den2 <- hist(mat2[,colu], breaks=seq(from=0, to=max(xRange), by=max(xRange)/8), ylim=c(0,100), col=col2, border=b2, add=TRUE, xlab="", ylab="", main="", xaxt="n", yaxt="n")
}


Diff <- flying / flightless

setwd("~/Dropbox/Research/FeatherEvolution/figures")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\figures")
library(vioplot)
pdf("rates_comparison.pdf", height=4.5, width=9)
par(mfrow=c(1,2), mar=c(4,5,1,1), cex.axis=0.5, las=1, mgp=c(2, 0.5, 0), tck=-0.01)
barplot(1-apply(wtMat, 2, median), ylim=c(0,1), ylab="avg Akaike weight for 2-rate BM model", col='black')
vioplot(log(Diff), col='white', colMed2='black', ylab="log flightless rate / volant rate") 
dev.off()

setwd("~/Dropbox/Research/FeatherEvolution/figures")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\figures")
pdf("rates_macro.pdf")
par(mfrow=c(2,2), las=1, bty="l", mar=c(3,3,1,1))
# macro
for (i in 1:4)	{
	makePlot(flightless, mat2=flying, colu = i, colOpt = 1 )
}
dev.off()


pdf("rates_feathers.pdf", height=6, width=12)
par(mfrow=c(3,6), las=1, bty="l", mar=c(3,3,1,1))
# macro
for (i in 5:21)	{
	makePlot(flightless, mat2=flying, colu = i, colOpt = 1 )
}
dev.off()

pdf("rates_featherPCs.pdf")
par(mfrow=c(2,2), las=1, bty="l", mar=c(3,3,1,1))
# macro
for (i in 22:25)	{
	makePlot(flightless, mat2=flying, colu = i, colOpt = 1 )
}
dev.off()

pdf("akaikewt_onerate.pdf")
# macro
barplot(apply(wtMat, 2, mean)[1:4], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
# angles
barplot(apply(wtMat, 2, mean)[5:6], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
# rachis
barplot(apply(wtMat, 2, mean)[7:10], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
# barb
barplot(apply(wtMat, 2, mean)[11:13], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
# barbule len
barplot(apply(wtMat, 2, mean)[14:15], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
# barbule density ti
barplot(apply(wtMat, 2, mean)[16:17], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
# barbule density middle
barplot(apply(wtMat, 2, mean)[18:19], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
# density int
barplot(apply(wtMat, 2, mean)[20:21], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
# PCs
barplot(apply(wtMat, 2, mean)[22:26], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
dev.off()

#########################################
##########       RAW       ##############
#########################################
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\output")
flightless <- read.csv("raw_lessMat.csv", row.names=1)
flying <- read.csv("raw_volMat.csv", row.names=1)
oneRate <- read.csv("raw_oneMat.csv", row.names=1)
AICcs <- read.csv("raw_aiccMat.csv", row.names=1)
oAICcs <- read.csv("raw_oaiccMat.csv", row.names = 1)
dAICcs <- oAICcs - AICcs
wtMat <- apply(dAICcs, 2, function(x) sapply(x, calcWt))

colnames(flightless) <- colnames(RawDat)[6:ncol(RawDat)]
colnames(wtMat) <- colnames(RawDat)[6:ncol(RawDat)]


makePlot <- function(mat1, mat2, colu, colOpt=1, Alpha=0.5)	{
	if (colOpt == 1)	{
		col1 <- rgb(179/255,88/255,6/255, Alpha)
		b1 <- rgb(179/255,88/255,6/255, 1)

		col2 <- rgb(84/255,39/255,136/255, Alpha)
		b2 <- rgb(84/255,39/255,136/255, 1)

	}
	else {
		col1 <- rgb(84/255,39/255,136/255, Alpha)
		b1 <- rgb(84/255,39/255,136/255, 1)

		col2 <- rgb(179/255,88/255,6/255, Alpha)
		b2 <- rgb(179/255,88/255,6/255, 1)

	}
	xRange <- range(c(mat1[,colu], mat2[,colu]))
	den1 <- hist(mat1[,colu], breaks=seq(from=0, to=max(xRange), by=max(xRange)/8), ylim=c(0,100), col=col1, border=b1, xlab="", ylab="", main=colnames(mat1)[colu])
	den2 <- hist(mat2[,colu], breaks=seq(from=0, to=max(xRange), by=max(xRange)/8), ylim=c(0,100), col=col2, border=b2, add=TRUE, xlab="", ylab="", main="", xaxt="n", yaxt="n")
}

setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\figures")
pdf("rates_raw_macro.pdf")
par(mfrow=c(2,2), las=1, bty="l", mar=c(3,3,1,1))
# macro
for (i in 1:4)	{
	makePlot(flightless, mat2=flying, colu = i, colOpt = 1 )
}
dev.off()


pdf("rates_raw_feathers.pdf", height=6, width=12)
par(mfrow=c(3,6), las=1, bty="l", mar=c(3,3,1,1))
# macro
for (i in 5:21)	{
	makePlot(flightless, mat2=flying, colu = i, colOpt = 1 )
}
dev.off()

pdf("akaikewt_raw_onerate.pdf")
# macro
barplot(apply(wtMat, 2, mean)[1:4], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
# angles
barplot(apply(wtMat, 2, mean)[5:6], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
# rachis
barplot(apply(wtMat, 2, mean)[7:10], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
# barb
barplot(apply(wtMat, 2, mean)[11:13], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
# barbule len
barplot(apply(wtMat, 2, mean)[14:15], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
# barbule density ti
barplot(apply(wtMat, 2, mean)[16:17], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
# barbule density middle
barplot(apply(wtMat, 2, mean)[18:19], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
# density int
barplot(apply(wtMat, 2, mean)[20:21], ylab="Akaike weight of 1-rate model", ylim=c(0,1), col='black')
dev.off()
