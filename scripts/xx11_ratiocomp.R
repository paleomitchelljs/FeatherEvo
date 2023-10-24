library(phytools)
setwd("~/Dropbox/Research/FeatherEvolution/scripts")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\scripts")
source("xx00_readFormatDat.R")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\scripts")
source("xx01_sisterPairs.R")
# s1 = flightless
# s2 = flighted
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

addSis <- function(trait, Xlim=c(-0.5,1.5), Shift=0.2)	{
	Vloc <- 0.3
	Lloc <- 0.7
	Fsize <- 0.5
	Ylim <- range(c(sisMat_vol[,trait],sisMat_less[,trait]), na.rm=T, finite=T)

	par(las=1, mgp=c(1.5, 0.5, 0), tck=-0.01)
	plot(1, 1, type="n", xaxt="n", ylim=Ylim, xlim=Xlim, xlab="", ylab="")
	axis(1, at=c(Vloc-Shift*Vloc, Lloc+Shift*Lloc), labels=c("volant","flightless"))
	for (j in 1:Npairs)	{
		Slope <- sisMat_vol[j,trait] - sisMat_less[j,trait]
		Col <- ifelse(Slope > 0, "#fc8d59", "#91bfdb")
		segments(Vloc, sisMat_vol[j,trait], Lloc, sisMat_less[j,trait], col=Col, lwd=2)
	}

	# place volant species tips
	vY <- spreadlabels(sisMat_vol[,trait], fsize=Fsize)
	silent <- sapply(1:length(vY), function(x) segments(Vloc, sisMat_vol[x,trait], Vloc-max(Fsize*strheight(rownames(sisMat_vol))), vY[x], col='gray50', lty=3))
	text(rep(Vloc, nrow(sisMat_vol))-(0.5*Fsize*strheight("a")), vY, rownames(sisMat_vol), cex=Fsize, pos=2, font=3)

	# place flightless species tips
	lY <- spreadlabels(sisMat_less[,trait], fsize=Fsize)
	silent <- sapply(1:length(lY), function(x) segments(Lloc, sisMat_less[x,trait], Lloc+max(Fsize*strheight(rownames(sisMat_less))), lY[x], col='gray50', lty=3))
	text(rep(Lloc, nrow(sisMat_less))+(0.5*Fsize*strheight("a")), lY, rownames(sisMat_less), cex=Fsize, pos=4, font=3)

}

Alpha <- 0.75
Cols <- c(rgb(179/255,88/255,6/255, Alpha),rgb(84/255,39/255,136/255, Alpha))

RawDat$barbuleDR <- RawDat$barbule_density_proximal / RawDat$barbule_density_distal

setwd("~/Dropbox/Research/FeatherEvolution/output")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\output")
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


hist(log(wing[,"s1"]/wing[,"s2"]))
hist(log(barb_ang[,"s1"]/barb_ang[,"s2"]))
hist(log(barb_len[,"s1"]/barb_len[,"s2"]))
hist(log(barb_dens[,"s1"]/barb_dens[,"s2"]))
hist(log(barbule_dens[,"s1"]/barbule_dens[,"s2"]))


Ndraw <- 1e2
Opac <- 10 / Ndraw
## Tri-plot
setwd("~/Dropbox/Research/FeatherEvolution/figures")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\figures")

pdf("triplots.pdf", height=12, width=9)
par(mfrow=c(4,3), mar=c(4,5,1,1), mgp=c(2,0.5,0), tck=-0.01, las=1)

#pdf("barbule_asymmetry_triplot.pdf", height=3, width=9)
#par(mfrow=c(1,3), mar=c(4,5,1,1), mgp=c(2,0.5,0), tck=-0.01, las=1)
boxplot(RawDat$barbDR~RawDat$flightless, boxwex=0.25, border=Cols, col="white", lty=1, ylab="barb density ratio (trail / lead)", xlab="", xaxt="n", yaxt="n")
axis(1, at=c(1, 2), label=c("volant", "flightless"))
axis(2)

# placeholder for sister-pair comparison plot
addSis("barbDR", Xlim=c(-1,2), Shift=0.5)
plot(RawDat$barb_density_t, RawDat$barb_density_l, pch=c(22,21)[RawDat$flightless+1], bg=Cols[RawDat$flightless+1], xlim=c(0,8), ylim=c(0,8), cex=1.5, xlab="trail barb density", ylab="lead barb density", type="n")
#abline(0,1,lty=3,lwd=2)
for (i in 1:Ndraw)	{
	Draw <- sample(1:length(barb_dens$b0), 1)
	addLine(RawDat$barb_density_t, 
		barb_dens$b0[Draw], 
		barb_dens$b1[Draw], 
		barb_dens$b2[Draw], 
		barb_dens$b3[Draw], 
		color=c(rgb(241/255,163/255,64/255,Opac),rgb(153/255,142/255,195/255,Opac))
	)
}

addLine(RawDat$barb_density_t, 
	mean(barb_dens$b0), 
	mean(barb_dens$b1), 
	mean(barb_dens$b2), 
	mean(barb_dens$b3), 
	color=c(rgb(241/255,163/255,64/255,1),rgb(153/255,142/255,195/255,1))
	)

points(RawDat$barb_density_t, RawDat$barb_density_l, pch=c(22,21)[RawDat$flightless+1], bg=Cols[RawDat$flightless+1])

#dev.off()

#par(mfrow=c(1,3)
#pdf("barblen_asymmetry_triplot.pdf", height=3, width=9)
#par(mfrow=c(1,3), mar=c(4,5,1,1), mgp=c(2.25,0.5,0), tck=-0.01, las=1)
boxplot(RawDat$barbLR~RawDat$flightless, boxwex=0.25, border=Cols, col="white", lty=1, ylab="barb length ratio (trailing / leading)", xlab="", xaxt="n", yaxt="n")
axis(1, at=c(1, 2), label=c("volant", "flightless"))
axis(2)

# placeholder for sister-pair comparison plot
addSis("barbLR", Xlim=c(-1,2), Shift=0.5)
plot(RawDat$barb_length_t, RawDat$barb_length_l, pch=c(22,21)[RawDat$flightless+1], bg=Cols[RawDat$flightless+1], xlim=c(0,80), ylim=c(0,80), cex=1.5, xlab="barb length (trailing)", ylab="barb length (leading)", type="n")
#abline(0,1,lty=3,lwd=2)

for (i in 1:Ndraw)	{
	Draw <- sample(1:length(barbule_dens$b0), 1)
	addLine(RawDat$barb_length_t, 
		barb_len$b0[Draw], 
		barb_len$b1[Draw], 
		barb_len$b2[Draw], 
		barb_len$b3[Draw], 
		color=c(rgb(241/255,163/255,64/255,Opac),rgb(153/255,142/255,195/255,Opac))
	)
}

addLine(RawDat$barb_length_t, 
	mean(barb_len$b0), 
	mean(barb_len$b1), 
	mean(barb_len$b2), 
	mean(barb_len$b3), 
	color=c(rgb(241/255,163/255,64/255,1),rgb(153/255,142/255,195/255,1))
	)

points(RawDat$barb_length_t, RawDat$barb_length_l, pch=c(22,21)[RawDat$flightless+1], bg=Cols[RawDat$flightless+1])

#dev.off()

#par(mfrow=c(1,3)
#pdf("barblen_asymmetry_triplot.pdf", height=3, width=9)
#par(mfrow=c(1,3), mar=c(4,5,1,1), mgp=c(2.25,0.5,0), tck=-0.01, las=1)
boxplot(RawDat$barbAng~RawDat$flightless, boxwex=0.25, border=Cols, col="white", lty=1, ylab="barb angle ratio (trailing / leading)", xlab="", xaxt="n", yaxt="n")
axis(1, at=c(1, 2), label=c("volant", "flightless"))
axis(2)

# placeholder for sister-pair comparison plot
addSis("barbAng", Xlim=c(-1,2), Shift=0.5)
plot(RawDat$barb_angle_t, RawDat$barb_angle_l, pch=c(22,21)[RawDat$flightless+1], bg=Cols[RawDat$flightless+1], xlim=c(0,80), ylim=c(0,80), cex=1.5, xlab="barb angle (trailing)", ylab="barb angle (leading)", type="n")
#abline(0,1,lty=3,lwd=2)

for (i in 1:Ndraw)	{
	Draw <- sample(1:length(barb_ang$b0), 1)
	addLine(RawDat$barb_angle_t, 
		barb_ang$b0[Draw], 
		barb_ang$b1[Draw], 
		barb_ang$b2[Draw], 
		barb_ang$b3[Draw], 
		color=c(rgb(241/255,163/255,64/255,Opac),rgb(153/255,142/255,195/255,Opac))
	)
}

addLine(RawDat$barb_angle_t, 
	mean(barb_ang$b0), 
	mean(barb_ang$b1), 
	mean(barb_ang$b2), 
	mean(barb_ang$b3), 
	color=c(rgb(241/255,163/255,64/255,1),rgb(153/255,142/255,195/255,1))
	)

points(RawDat$barb_angle_t, RawDat$barb_angle_l, pch=c(22,21)[RawDat$flightless+1], bg=Cols[RawDat$flightless+1])

#dev.off()

#pdf("wing-massRatio_triplot.pdf", height=3, width=9)
#par(mfrow=c(1,3), mar=c(4,5,1,1), mgp=c(2.25,0.5,0), tck=-0.01, las=1)
boxplot(RawDat$wingR~RawDat$flightless, boxwex=0.25, border=Cols, col="white", lty=1, ylab="wing-mass ratio", xlab="", xaxt="n", yaxt="n")
axis(1, at=c(1, 2), label=c("volant", "flightless"))
axis(2)

# placeholder for sister-pair comparison plot
addSis(3, Xlim=c(-1,2), Shift=0.5)
plot(RawDat$mass, RawDat$wing, pch=c(22,21)[RawDat$flightless+1], bg=Cols[RawDat$flightless+1], xlim=c(2,12), ylim=c(3.5,8), cex=1.5, xlab="mass", ylab="wing", type="n")
#abline(0,1,lty=3,lwd=2)
for (i in 1:Ndraw)	{
	Draw <- sample(1:length(wing$b0), 1)
	addLine(RawDat$mass, 
		wing$b0[Draw], 
		wing$b1[Draw], 
		wing$b2[Draw], 
		wing$b3[Draw], 
		color=c(rgb(241/255,163/255,64/255,Opac),rgb(153/255,142/255,195/255,Opac))
	)
}

addLine(RawDat$mass, 
	mean(wing$b0), 
	mean(wing$b1), 
	mean(wing$b2), 
	mean(wing$b3), 
	color=c(rgb(241/255,163/255,64/255,1),rgb(153/255,142/255,195/255,1))
	)

points(RawDat$mass, RawDat$wing, pch=c(22,21)[RawDat$flightless+1], bg=Cols[RawDat$flightless+1])

dev.off()


#par(mfrow=c(1,3)
#pdf("barblen_asymmetry_triplot.pdf", height=3, width=9)
#par(mfrow=c(1,3), mar=c(4,5,1,1), mgp=c(2.25,0.5,0), tck=-0.01, las=1)
boxplot(RawDat$barbAng~RawDat$flightless, boxwex=0.25, border=Cols, col="white", lty=1, ylab="barb length ratio (leading / trailing)", xlab="", xaxt="n", yaxt="n")
axis(1, at=c(1, 2), label=c("volant", "flightless"))
axis(2)

# placeholder for sister-pair comparison plot
addSis(12, Xlim=c(-1,2), Shift=0.5)
plot(RawDat$barb_angle_t, RawDat$barb_angle_l, pch=c(22,21)[RawDat$flightless+1], bg=Cols[RawDat$flightless+1], xlim=c(0,80), ylim=c(0,80), cex=1.5, xlab="barb length (trailing)", ylab="barb length (leading)", type="n")
#abline(0,1,lty=3,lwd=2)

for (i in 1:Ndraw)	{
	Draw <- sample(1:length(barb_ang$b0), 1)
	addLine(RawDat$barb_angle_t, 
		barb_ang$b0[Draw], 
		barb_ang$b1[Draw], 
		barb_ang$b2[Draw], 
		barb_ang$b3[Draw], 
		color=c(rgb(241/255,163/255,64/255,Opac),rgb(153/255,142/255,195/255,Opac))
	)
}

addLine(RawDat$barb_angle_t, 
	mean(barb_ang$b0), 
	mean(barb_ang$b1), 
	mean(barb_ang$b2), 
	mean(barb_ang$b3), 
	color=c(rgb(241/255,163/255,64/255,1),rgb(153/255,142/255,195/255,1))
	)

points(RawDat$barb_angle_t, RawDat$barb_angle_l, pch=c(22,21)[RawDat$flightless+1], bg=Cols[RawDat$flightless+1])

dev.off()

library(vioplot)
pdf("triplot_rates.pdf", height = 4, width = 8)
par(mar=c(4,5,1,1), mgp=c(2,0.5,0), tck=-0.01, las=1)

vioplot(log(wing$s1) - log(wing$s2), log(barb_len$s1) - log(barb_len$s2), log(barb_dens$s1) - log(barb_dens$s2), log(barb_ang$s1) - log(barb_ang$s2), names=c("wing", "barb len", "barb dens", "barb ang"), col='white', ylab="log rate ratio (flightless / volant)")
dev.off()

par(mfrow=c(2,2), boxwex=0.25)
hist(wing$s1, xlim=c(0,0.005), main="", xlab="", ylab="", col=Cols[1])
hist(barbule_dens$s1, xlim=c(0,10), main="", xlab="", ylab="", col=Cols[1])
hist(wing$s2, xlim=c(0,0.005), main="", xlab="", ylab="", col=Cols[2])
hist(barbule_dens$s2, xlim=c(0,10), main="", xlab="", ylab="", col=Cols[2])

boxplot(wing$s1, wing$s2, barbule_dens$s1, barbule_dens$s2)

