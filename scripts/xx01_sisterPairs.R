library(phytools)setwd("~/Dropbox/Research/FeatherEvolution/scripts")
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\scripts")
source("spreadlabels.R")
source("xx00_readFormatDat.R")

setwd("~/Dropbox/Research/FeatherEvolution/output")
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\output")
# Set Mesitornis as flightless since it may as well be
# specifically check barbule length + rachis width for semiaquatic, barb length in terrestrial365

#BaseCors <- cor(Flight, Feathers, use="complete")

# Make sure to specifically test rachis width, barbule length, and 
#### Sister Lineages
Npairs <- max(Pairs$pair, na.rm=T)

# Do all traits
sCol <- which(colnames(Dat) == "mass")
eCol <- ncol(Dat)
nCols <- length(sCol:eCol)


Dat <- ModelDat
# Choose traits
Dat$barbuleDR <- log(Dat$barbule_density_proximal) - log(Dat$barbule_density_distal)
Alt <- c("mass", "wingR", "tarsusR", "barbDR", "barbLR", "barbL_lead", "barbAng", "rachisWL", "tailR", "barbL_trail", "barbuleL_distrail", "barbAng", "rachis_width", "barb_density_l", "barb_density_t", "barb_length_l", "barb_length_t", "barb_angle_l", "barb_angle_t", "barbuleDR")
nCols <- length(Alt)

sisMat_vol <- matrix(0, nrow=Npairs, ncol=nCols+1)
sisMat_less <- matrix(0, nrow=Npairs, ncol=nCols+1)

vNames <- c()
lNames <- c()
for (i in 1:Npairs)	{
	Choose <- which(Pairs$pair == i)
	NamesFlight <- Pairs[Choose,c(1,2)]
#	sisMat_vol[i,] <- c(i,unlist(RawDat[which(RawDat$species == NamesFlight[which(NamesFlight[,2]==1),1]),sCol:eCol]))
	sisMat_vol[i,] <- c(i,unlist(Dat[which(Dat$species == NamesFlight[which(NamesFlight[,2]==1),1]),Alt]))

	vNames[i] <- NamesFlight[which(NamesFlight[,2]==1),1]

#	sisMat_less[i,] <- c(i,unlist(RawDat[which(RawDat$species == NamesFlight[which(NamesFlight[,2]==0),1]),sCol:eCol]))
	sisMat_less[i,] <- c(i,unlist(Dat[which(Dat$species == NamesFlight[which(NamesFlight[,2]==0),1]),Alt]))

	lNames[i] <- NamesFlight[which(NamesFlight[,2]==0),1]

}

rownames(sisMat_vol) <- vNames
rownames(sisMat_less) <- lNames

#colnames(sisMat_vol) <- c("pair", colnames(RawDat)[sCol:eCol])
#colnames(sisMat_less) <- c("pair", colnames(RawDat)[sCol:eCol])

colnames(sisMat_vol) <- c("pair", Alt)
colnames(sisMat_less) <- c("pair", Alt)

write.csv(sisMat_vol, "sisMat_vol_alt.csv")
write.csv(sisMat_less, "sisMat_less_alt.csv")


notsis <- setdiff(Tree$tip.label, Pairs[!is.na(Pairs$pair),1])
sisterTree <- drop.tip(Tree, notsis)

setwd("~/Dropbox/Research/FeatherEvolution/figures")
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\figures")
pdf("sister_comparisons_alt.pdf", height=4, width=8)
Vloc <- 0.35
Lloc <- 0.65
Fsize <- 0.3

Names <- c("", "mass", "wing:mass", "tarsus:mass", "barb density t:l", "barb length t:l", "barb len:rachis width (l)", "barb angle t:l")
par(mfrow=c(2,4), las=1, mgp=c(1.5, 0.5, 0), tck=-0.01, mar=c(4,3,2,1))
plot(sisterTree, cex=0.45, no.margin=F)
for (trait in 2:8)	{
	Ylim <- range(c(sisMat_vol[,trait],sisMat_less[,trait]), na.rm=T, finite=T)

#	par(las=1, mgp=c(1.5, 0.5, 0), tck=-0.01)
	plot(1, 1, type="n", xlim=c(0, 1), xaxt="n", ylim=Ylim, xlab="", ylab="", main=Names[trait])
	axis(1, at=c(Vloc, Lloc), labels=c("vol.","not vol."))
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
dev.off()

################# 
#################
#################
#################
#pdf("sis_div_time_cors.pdf")
#for (i in 2:ncol(siscompMat))	{
#	par(las=0)
#	plot(siscompMat[,1]/1e6, siscompMat[,i], ylab=colnames(siscompMat)[i], xlab="age (Ma)")
#}
#dev.off()
