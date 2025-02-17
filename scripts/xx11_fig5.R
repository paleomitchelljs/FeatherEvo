library(phytools)
setwd("~/Dropbox/Research/FeatherEvolution/scripts")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\scripts")
source("xx00_readFormatDat.R")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\scripts")
source("xx01_sisterPairs.R")
# b0 = intercept
# b1 = volant slope
# b2 = flightless intercept
# b3 = flightless slope
# s1 = volant rate
# s2 = flightless rate
addLine <- function(x1, xint, xb1, xb2, xb3, color=c("red","red"), Lty=1, Lwd=2)	{
	Start <- min(x1)
	Stop <- max(x1)
	yStart1 <- (Start*xb1) + xint
	yStart2 <- (Start*xb1) + (Start*xb3) + xb2 + xint

	yStop1 <- (Stop*xb1) + xint
	yStop2 <- (Stop*xb1) + xb2 + (xb3*Stop) + xint

	segments(Start,yStart1,Stop,yStop1, col=color[1], lty=Lty, lwd=Lwd)
	segments(Start,yStart2,Stop,yStop2, col=color[2], lty=Lty, lwd=Lwd)
}

addSis_old <- function(trait, Xlim=c(-3,4), Shift=2)	{
	Vloc <- -1
	Lloc <- 1
	Fsize <- 0.8
	Ylim <- range(c(sisMat_vol[,trait],sisMat_less[,trait]), na.rm=T, finite=T)

	#par(las=1, mgp=c(1.5, 0.5, 0), tck=-0.01)
	plot(1, 1, type="n", xaxt="n", ylim=Ylim, xlim=Xlim, xlab="", ylab="")
	axis(1, at=c(Vloc-Shift*Vloc, Lloc+Shift*Lloc), labels=c("Vol.","Nonvol."))
	Npairs <- max(sisMat_vol[,1])
	Slope <- c()
	for (j in 1:Npairs)	{
		Slope[j] <- sisMat_vol[j,trait] - sisMat_less[j,trait]
		Col <- ifelse(Slope[j] > 0, "black", "red")
		segments(Vloc, sisMat_vol[j,trait], Lloc, sisMat_less[j,trait], col=Col, lwd=1.25)
	}

	# place volant species tips
	vY <- spreadlabels(sisMat_vol[,trait], fsize=Fsize, rng = range(c(sisMat_less[,trait], sisMat_vol[,trait])))
	silent <- sapply(1:length(vY), function(x) segments(Vloc, sisMat_vol[x,trait], Vloc-max(Fsize*strheight(rownames(sisMat_vol))), vY[x], col='gray70', lty=1, lwd=0.5))
	text(rep(Vloc, nrow(sisMat_vol))-(0.5*Fsize*strheight("a")), vY, abbreviate(rownames(sisMat_vol)), cex=Fsize, pos=2, font=3)

	# place flightless species tips
	lY <- spreadlabels(sisMat_less[,trait], fsize=Fsize, rng = range(c(sisMat_less[,trait], sisMat_vol[,trait])))
	silent <- sapply(1:length(lY), function(x) segments(Lloc, sisMat_less[x,trait], Lloc+max(Fsize*strheight(rownames(sisMat_less))), lY[x], col='gray70', lty=1, lwd=0.5))
	text(rep(Lloc, nrow(sisMat_less))+(0.5*Fsize*strheight("a")), lY, abbreviate(rownames(sisMat_less)), cex=Fsize, pos=4, font=3)

}


addSis <- function(trait, Xlim=c(-1.25,2), Shift=1, Xlab1=-1, Xlab2=0.5, addLabs = TRUE, xshift = 0.1)	{
	Vloc <- -1
	Lloc <- 0.5
	Fsize <- 0.75
	Ylim <- range(c(sisMat_vol[,trait],sisMat_less[,trait]), na.rm=T, finite=T)

	#par(las=1, mgp=c(1.5, 0.5, 0), tck=-0.01)
	plot(1, 1, type="n", xaxt="n", ylim=Ylim, xlim=Xlim, xlab="", ylab="", bty="l")
	axis(1, at=c(Xlab1, Xlab2), labels=c("Vol.","Nonvol."))
	Slope <- c()
	Col <- c()
	for (j in 1:nrow(sisMat_vol))	{
		Slope[j] <- sisMat_vol[j,trait] - sisMat_less[j,trait]
		Col[j] <- ifelse(Slope[j] > 0, "black", "red")
		segments(Vloc, sisMat_vol[j,trait], Lloc, sisMat_less[j,trait], col=Col[j], lwd=1.25)
	}

	if (addLabs)	{
	# place flightless species labels
		lY <- spreadlabels(sisMat_less[,trait], fsize=Fsize, rng = range(c(sisMat_less[,trait], sisMat_vol[,trait])))
		silent <- sapply(1:nrow(sisMat_less), function(x) segments(Lloc, sisMat_less[x,trait], Lloc+max(Fsize*strheight(rownames(sisMat_less))), lY[x], col='gray70', lty=1, lwd=0.5))
		text(rep(Lloc, nrow(sisMat_less))+(0.2*Fsize*strheight("a")), lY, abbreviate(rownames(sisMat_less)), cex=Fsize, pos=4, font=3, xpd=NA, col=Col)
#		text(rep(Lloc, nrow(sisMat_less))+(0.5*Fsize*strheight("a")), lY, sapply(rownames(sisMat_less), function(x) strsplit(x, " ")[[1]][1]), cex=Fsize, pos=4, font=3, xpd=NA, col = Col)
	}
}

makeReg <- function(x, y, flight, regression, Ndraw = 1e2, Xlab="", Ylab="", Xline = 1.5, Yline = 2, labcex = 1, legLoc = "topleft")	{
	Opac <- 10 / Ndraw
	Cols <- c(rgb(241/255,163/255,64/255,Opac),rgb(153/255,142/255,195/255,Opac))
	plot(x, y, pch=c(22,21)[flight+1], bg=Cols[flight+1], xlim=range(pretty(x)), ylim=range(pretty(y)), cex=1.5, xlab="", ylab="", type="n", bty='l')
	mtext(Xlab, side = 1, line = Xline, cex = labcex)
	mtext(Ylab, side = 2, line = Yline, las = 0, cex = labcex)
	for (i in 1:Ndraw)	{
		Draw <- sample(1:length(regression$b0), 1)
		addLine(x, 
			regression$b0[Draw], 
			regression$b1[Draw], 
			regression$b2[Draw], 
			regression$b3[Draw], 
			color=Cols
		)
	}

	addLine(x, 
		mean(regression$b0), 
		mean(regression$b1), 
		mean(regression$b2), 
		mean(regression$b3), 
		color=c(rgb(241/255,163/255,64/255,1),rgb(153/255,142/255,195/255,1))
		)

	points(x, y, pch=c(22,21)[flight+1], bg=c(rgb(241/255,163/255,64/255,1),rgb(153/255,142/255,195/255,1))[flight+1])
	legend("topleft", cex=0.75, bty="n", legend=c("Vol.", "Nonvol."), pch=c(15,16), col=c(rgb(241/255,163/255,64/255,1),rgb(153/255,142/255,195/255,1)), text.col=c(rgb(241/255,163/255,64/255,1),rgb(153/255,142/255,195/255,1)), xpd=NA)

#	Tf <- diff(pretty(y))[1]
#	legLoc <- c(min(pretty(x)), max(pretty(y))+0.5*Tf)
#	legend(legLoc, bty="n", legend=c("Vol.", "Nonvol."), pch=c(15,16), col=c(rgb(241/255,163/255,64/255,1),rgb(153/255,142/255,195/255,1)), xpd=NA)
}

makeHist <- function(Trait, Xlab="", Ylab="", Xline = 1.5, Yline = 2, labcex = 1)	{
#	Trait$s1 <- Trait$s1^2
#	Trait$s2 <- Trait$s2^2
	Xlim <- range(c(log(Trait$s1), log(Trait$s2)))
	Bs <- pretty(Xlim)
	Ylim <- c(0, 1)
	Alpha <- 0.25
	H <- hist(log(Trait$s1), plot=F, breaks=Bs)
	H$counts <- H$counts / sum(H$counts)
	plot(H, border=rgb(0,0,0,0), col=rgb(179/255,88/255, 6/255, Alpha), xlim=range(Bs), ylim=Ylim, ylab="", yaxt="n", xlab="", main="")
	mtext(Xlab, side = 1, line = Xline, cex = labcex)
	mtext(Ylab, side = 2, line = Yline, las = 0, cex = labcex)

	lines(H$breaks, c(H$counts, 0), lwd=2, col=rgb(179/255,88/255, 6/255, 1), type="s")
	axis(2)

	Alpha <- 0.25
	H <- hist(log(Trait$s2), plot=F, breaks=Bs)
	H$counts <- H$counts / sum(H$counts)
	plot(H, border=rgb(0,0,0,0), col=rgb(84/255,39/255,136/255, Alpha), xlim=range(Bs), add=T, main="", xlab="", ylab="")
	lines(H$breaks, c(H$counts, 0), lwd=2, col=rgb(84/255,39/255,136/255, 1), type="s")
}

rownames(sisMat_vol)[grep("occipitalis occipitalis", rownames(sisMat_vol))] <- "Podiceps occipitalis"

RawDat$barbuleDR <- RawDat$barbule_density_proximal / RawDat$barbule_density_distal

setwd("~/Dropbox/Research/FeatherEvolution/output")
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\output")
barb_ang <- read.csv("barb_ang_MM.csv")
barb_ang <- barb_ang[floor(0.5*nrow(barb_ang)):nrow(barb_ang),]
barb_dens <- read.csv("barb_dens_MM.csv")
barb_dens <- barb_dens[floor(0.5*nrow(barb_dens)):nrow(barb_dens),]
barb_len <- read.csv("barb_len_MM.csv")
barb_len <- barb_len[floor(0.5*nrow(barb_len)):nrow(barb_len),]
barbule_dens <- read.csv("barbule_dens_MM.csv")
barbule_dens <- barbule_dens[floor(0.5*nrow(barbule_dens)):nrow(barbule_dens),]
wing <- read.csv("wing_MM.csv")
wing <- wing[floor(0.5*nrow(wing)):nrow(wing),]


WR <- log(wing$s1) - log(wing$s2)
length(which(WR > 0)) / length(WR)
mean(WR)

BAR <- log(barb_ang$s1) - log(barb_ang$s2)

BDR <- log(barb_dens$s1) - log(barb_dens$s2)
length(which(BDR > 0)) / length(BDR)
mean(BDR)
mean(barb_dens$lam)

## Make the plot
setwd("~/Dropbox/Research/FeatherEvolution/figures")
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\figures")

pdf("Fig5_sislabs_onlyPos.pdf", height=6.65, width=6.65)
par(mfrow=c(4,4), mar=c(3,4,0.5,0.5), mgp=c(2,0.35,0), tck=-0.01, las=1, oma=c(0, 0, 4, 0))
Ndraw <- 1e2
Opac <- 10 / Ndraw
LWD <- 2
Cols <- c(rgb(241/255,163/255,64/255,1),rgb(153/255,142/255,195/255,1))
Labcex <- 0.75
Labs <- c("Vol.", "Nonvol.")
OMACex <- 1
OMALine <- 1
boxplot(ModelDat$barbDR~ModelDat$flightless, boxwex=0.25, border=Cols, col="white", lty=1, ylab="", xlab="", xaxt="n", yaxt="n", lwd=LWD, bty='l')
mtext("Direct Comparisons", side = 3, line = OMALine, font=2, outer = F, cex = OMACex)
mtext("Barb density ratio (t/l)", side = 2, line=par("mgp")[1], cex=Labcex, las=0)
axis(1, at=c(1, 2), label=Labs)
axis(2)

addSis("barbDR")
mtext("Barb density ratio (t/l)", side=2, line=par("mgp")[1], cex=Labcex, las=0)
mtext("Sister-Taxa\nComparisons", side = 3, line = OMALine, font=2, outer = F, cex = OMACex)

makeReg(RawDat$barb_density_t, RawDat$barb_density_l, flight = RawDat$flightless, regression = barb_dens, Xlab="Trail barb density (#/mm)", Ylab = "Lead barb density (#/mm)", labcex = Labcex, legLoc=c(2, 4))
mtext("Phylogenetic\nRegressions", side = 3, line = OMALine, font=2, outer = F, cex = OMACex)

makeHist(barb_dens, Xlab = expression(paste("log est. rate (", gamma, ")", sep="")), Ylab = "Proportion", labcex = Labcex)
mtext("Rate\nComparisons", side = 3, line = OMALine, font=2, outer = F, cex = OMACex)


### Barb LR
boxplot(ModelDat$barbLR~ModelDat$flightless, boxwex=0.25, border=Cols, col="white", lty=1, ylab="", xlab="", xaxt="n", yaxt="n", lwd=LWD, bty='l')
mtext("Barb length ratio (t/l)", side = 2, line=par("mgp")[1], cex=Labcex, las=0)
axis(1, at=c(1, 2), label=Labs)
axis(2)

addSis("barbLR")
mtext("Barb length ratio (t/l)", side = 2, line=par("mgp")[1], cex=Labcex, las=0)

makeReg(RawDat$barb_length_t, RawDat$barb_length_l, flight = RawDat$flightless, regression = barb_len, Xlab="Trail barb length (mm)", Ylab = "Lead barb length (mm)", labcex = Labcex, legLoc = c(10, 100))
makeHist(barb_len, Xlab = expression(paste("log est. rate (", gamma, ")", sep="")), Ylab = "Proportion", labcex = Labcex)

### Barb Angle
boxplot(ModelDat$barbAng~ModelDat$flightless, boxwex=0.25, border=Cols, col="white", lty=1, ylab="", xlab="", xaxt="n", yaxt="n", lwd=LWD, bty='l')
mtext("Barb angle ratio (t/l)", side = 2, line=par("mgp")[1], cex=Labcex, las=0)
axis(1, at=c(1, 2), label=Labs)
axis(2)

addSis("barbAng")
mtext("Barb angle ratio (t/l)", side = 2, line=par("mgp")[1], cex=Labcex, las=0)

makeReg(RawDat$barb_angle_t, RawDat$barb_angle_l, flight = RawDat$flightless, regression = barb_ang, Xlab="Trail barb angle", Ylab = "Lead barb angle", labcex = Labcex, legLoc=c(10, 60))
makeHist(barb_ang, Xlab = expression(paste("log est. rate (", gamma, ")", sep="")), Ylab = "Proportion", labcex = Labcex)


### Wing
boxplot(RawDat$wingR~RawDat$flightless, boxwex=0.25, border=Cols, col="white", lty=1, ylab="", xlab="", xaxt="n", yaxt="n", lwd=LWD, bty='l')
mtext("Rel. wing size (length/mass)", side = 2, line=par("mgp")[1], cex=Labcex, las=0)
axis(1, at=c(1, 2), label=Labs)
axis(2)

# placeholder for sister-pair comparison plot
addSis("wingR")
mtext("Wing:Mass", side = 2, line=par("mgp")[1], cex=Labcex, las=0)

makeReg(RawDat$mass, RawDat$wing, flight = RawDat$flightless, regression = wing, Xlab="Mass (log g)", Ylab = "Wing length (log mm)", labcex = Labcex, legLoc = c(3, 7.5))
makeHist(wing, Xlab = expression(paste("log est. rate (", gamma, ")", sep="")), Ylab = "Proportion", labcex = Labcex)

dev.off()
