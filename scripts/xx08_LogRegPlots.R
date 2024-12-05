library(phytools)
library(brms)
library(bayesplot)

setwd("~/Dropbox/Research/FeatherEvolution/scripts")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\scripts")
source("spindlePlot.R")
library(vioplot)

source("xx00_readFormatDat.R")


getPval <- function(x)	{
	mu <- mean(x)
	if (mu > 0)	{
		out <- 1 - length(which(x > 0))/length(x)
	}
	if (mu < 0)	{
		out <- 1 - length(which(x < 0))/length(x)		
	}
	return(out)
}

setwd("~/Dropbox/Research/FeatherEvolution/output")
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\output")

model_logreg_phy <- readRDS("model_logreg_phy.rds")
model_logreg_phy_FULL <- readRDS("model_logreg_phy_FULL.rds")
model_logreg_phy_macro <- readRDS("model_logreg_phy_macro.rds")
model_logreg_phy_len <- readRDS("model_logreg_phy_len.rds")
model_logreg_phy_len_red <- readRDS("model_logreg_phy_len_red.rds")
model_logreg_phy_len_pcm_traits <- readRDS("model_logreg_phy_len_pcm_traits.rds")


resTable <- summary(model_logreg_phy_len_pcm_traits)
resTable <- resTable$fixed
resTable <- apply(resTable, 2, round, digits=2)
write.csv(resTable, "resTable.csv")

setwd("~/Dropbox/Research/FeatherEvolution/figures")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\figures")
pdf("postpredcheck_pcm.pdf")
#pp_check(model_logreg_phy, ndraws=100)
#pp_check(model_logreg_phy_FULL, ndraws=100, type="scatter_avg")
#pp_check(model_logreg_phy_macro, ndraws=100)
pp_check(model_logreg_phy_len_red, ndraws=100)
pp_check(model_logreg_phy_len_red, ndraws=100, type="scatter_avg")
dev.off()

loo_FULL <- loo(model_logreg_phy_FULL)
loo_pc <- loo(model_logreg_phy)
loo_mac <- loo(model_logreg_phy_macro)
loo_len <- loo(model_logreg_phy_len)
loo_len_red <- loo(model_logreg_phy_len_red)

phy_post <- posterior_samples(model_logreg_phy)
phy_FULL_post <- posterior_samples(model_logreg_phy_FULL)
phy_macro_post <- posterior_samples(model_logreg_phy_macro)
phy_len_post <- posterior_samples(model_logreg_phy_len)
phy_red_post <- posterior_samples(model_logreg_phy_len_red)
phy_pcm_post <- posterior_samples(model_logreg_phy_len_pcm_traits)

par(mfrow=c(1,2), mar=c(4,4,1,1), mgp=c(1.5, 0.5, 0), tck=-0.01)
# PC logreg
plot(1, 1, xlim=c(0,5), ylim=c(-25,25), type="n", xaxt="n", xlab="", ylab="")
vioplot(phy_post$b_mass, at=1, add=TRUE, col='#fb8072')
vioplot(phy_post$b_wing, at=2, add=TRUE, col='#fb8072')
vioplot(phy_post$b_tail, at=3, add=TRUE, col='#fb8072')
vioplot(phy_post$b_tarsus, at=4, add=TRUE, col='#fb8072')
abline(h=0, col='gray50')
axis(1, at=c(1:4), labels=c("mass", "wing", "tail", "tarsus"), tck=-0.03)

plot(1, 1, xlim=c(0,5), ylim=c(-2,2), type="n", xaxt="n", xlab="", ylab="")
vioplot(phy_post$b_PC1, at=1, add=TRUE, col="#8dd3c7")
vioplot(phy_post$b_PC2, at=2, add=TRUE, col="#8dd3c7")
vioplot(phy_post$b_PC3, at=3, add=TRUE, col="#8dd3c7")
vioplot(phy_post$b_PC4, at=4, add=TRUE, col="#8dd3c7")
abline(h=0, col='gray50')
axis(1, at=c(1:4), labels=c("PC1", "PC2", "PC3", "PC4"), tck=-0.03)

####################
### no PC logreg ###
####################
###############################################################
Mat <- phy_FULL_post

pdf("model_logreg_phy_FULL_regCoefficients.pdf", height=4, width=12)
Pvals <- apply(Mat, 2, getPval)
Thresh <- c(1, 0.1, 0.05)

Cols <- c("#e0ecf4", "#9ebcda", "#8856a7")
Buffer <- -1

Chosen <- c("b_mass", "b_wing", "b_tail", "b_tarsus")
par(mfrow=c(1,2), mar=c(5,3,0.5,0.5), mgp=c(1.5, 0.5, 0), tck=-0.01, las=1)

plot(1, 1, xlim=c(0,length(Chosen)+1), ylim=c(-40,40), type="n", xaxt="n", xlab="", ylab="")
abline(h=0, col='gray50')
for (i in 1:length(Chosen))	{
	vioplot(Mat[,Chosen[i]], at=i, add=TRUE, col=Cols[sum(Thresh>Pvals[Chosen[i]])])
}
axis(1, at=c(1:length(Chosen)), labels=gsub("b\\_", "", Chosen), tck=-0.02)
legend("bottomleft", bty="n", pch=21, pt.bg=Cols, legend=c("p > 0.1", "p < 0.1", "p < 0.05"), pt.cex=2)

Chosen <- c("b_rachis_width", "b_rachis_length", "b_barb_length_t", "b_barb_length_l", "b_barb_angle_t", "b_barb_angle_l")
plot(1, 1, xlim=c(0,length(Chosen)+1), ylim=c(-10,10), type="n", xaxt="n", xlab="", ylab="")
abline(h=0, col='gray50')
for (i in 1:length(Chosen))	{
	vioplot(Mat[,Chosen[i]], at=i, add=TRUE, col=Cols[sum(Thresh>Pvals[Chosen[i]])])
}
axis(1, at=c(1:length(Chosen)), labels=c("","","","","",""), tck=-0.01)
ToMark <- 1:length(Chosen)
Names <- gsub("b\\_", "", Chosen)
silent <- sapply(1:length(ToMark), function(x) text(ToMark[x], par("usr")[3]+Buffer, Names[x], adj = 1, srt = 35, xpd=NA, cex=1))

dev.off()


###############################################################
setwd("~/Dropbox/Research/FeatherEvolution/figures")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\figures")
Mat <- phy_len_post

pdf("model_logreg_phy_len_regCoefficients.pdf", height=4, width=12)

Pvals <- apply(Mat, 2, getPval)
Thresh <- c(1, 0.1, 0.05)

Cols <- c("#e0ecf4", "#9ebcda", "#8856a7")
Buffer <- -1

Chosen <- c("b_mass", "b_wing", "b_tail", "b_tarsus")
par(mfrow=c(1,2), mar=c(5,3,0.5,0.5), mgp=c(1.5, 0.5, 0), tck=-0.01, las=1)

plot(1, 1, xlim=c(0,length(Chosen)+1), ylim=c(-40,40), type="n", xaxt="n", xlab="", ylab="")
abline(h=0, col='gray50')
for (i in 1:length(Chosen))	{
	vioplot(Mat[,Chosen[i]], at=i, add=TRUE, col=Cols[sum(Thresh>Pvals[Chosen[i]])])
}
axis(1, at=c(1:length(Chosen)), labels=gsub("b\\_", "", Chosen), tck=-0.02)
legend("bottomleft", bty="n", pch=21, pt.bg=Cols, legend=c("p > 0.1", "p < 0.1", "p < 0.05"), pt.cex=2)

Chosen <- c("b_rachis_width", "b_rachis_length", "b_barb_length_t", "b_barb_length_l", "b_barb_angle_t", "b_barb_angle_l")
plot(1, 1, xlim=c(0,length(Chosen)+1), ylim=c(-10,10), type="n", xaxt="n", xlab="", ylab="")
abline(h=0, col='gray50')
for (i in 1:length(Chosen))	{
	vioplot(Mat[,Chosen[i]], at=i, add=TRUE, col=Cols[sum(Thresh>Pvals[Chosen[i]])])
}
axis(1, at=c(1:length(Chosen)), labels=c("","","","","",""), tck=-0.01)
ToMark <- 1:length(Chosen)
Names <- gsub("b\\_", "", Chosen)
silent <- sapply(1:length(ToMark), function(x) text(ToMark[x], par("usr")[3]+Buffer, Names[x], adj = 1, srt = 35, xpd=NA, cex=1))

dev.off()

##############################################################
setwd("~/Dropbox/Research/FeatherEvolution/figures")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\figures")

Mat <- phy_red_post
pdf("model_logreg_phy_red_regCoefficients.pdf", height=4, width=12)

Pvals <- apply(Mat, 2, getPval)
Thresh <- c(1, 0.1, 0.05)

Cols <- c("#e0ecf4", "#9ebcda", "#8856a7")
Buffer <- -1

Chosen <- c("b_mass", "b_wing", "b_tarsus", "b_tail")
par(mfrow=c(1,2), mar=c(5,3,0.5,0.5), mgp=c(1.5, 0.5, 0), tck=-0.01, las=1)

plot(1, 1, xlim=c(0,length(Chosen)+1), ylim=c(-20,20), type="n", xaxt="n", xlab="", ylab="")
abline(h=0, col='gray50')
for (i in 1:length(Chosen))	{
	vioplot(Mat[,Chosen[i]], at=i, add=TRUE, col=Cols[sum(Thresh>Pvals[Chosen[i]])])
}
axis(1, at=c(1:length(Chosen)), labels=gsub("b\\_", "", Chosen), tck=-0.02)
legend("bottomleft", bty="n", pch=21, pt.bg=Cols, legend=c("p > 0.1", "p < 0.1", "p < 0.05"), pt.cex=2)

Chosen <- c("b_barb_length_t", "b_barb_length_l", "b_barb_density_t", "b_barb_density_l", "b_barb_angle_t", "b_barb_angle_l")
plot(1, 1, xlim=c(0,length(Chosen)+1), ylim=c(-15,15), type="n", xaxt="n", xlab="", ylab="")
abline(h=0, col='gray50')
for (i in 1:length(Chosen))	{
	vioplot(Mat[,Chosen[i]], at=i, add=TRUE, col=Cols[sum(Thresh>Pvals[Chosen[i]])])
}
axis(1, at=c(1:length(Chosen)), labels=c("","","","","",""), tck=-0.01)
ToMark <- 1:length(Chosen)
Names <- gsub("b\\_", "", Chosen)
Names <- c("barb length (t)", "barb length (l)", "barb density (t)", "barb density (l)", "barb angle (t)", "barb angle (l)")
silent <- sapply(1:length(ToMark), function(x) text(ToMark[x], par("usr")[3]+Buffer, Names[x], adj = 1, srt = 35, xpd=NA, cex=1))

dev.off()

##############################################################
setwd("~/Dropbox/Research/FeatherEvolution/figures")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\figures")

Mat <- phy_pcm_post
pdf("model_logreg_phy_pcm_regCoefficients.pdf", height=4, width=8)

Pvals <- apply(Mat, 2, getPval)
Thresh <- c(1, 0.1, 0.05)

Cols <- c("#e0ecf4", "#9ebcda", "#8856a7")
Buffer <- -1

Chosen <- c("b_mass", "b_wingR", "b_tarsusR", "b_tailR")
par(mfrow=c(1,2), mar=c(5,3,0.5,0.5), mgp=c(1.5, 0.5, 0), tck=-0.01, las=1)

plot(1, 1, xlim=c(0,length(Chosen)+1), ylim=c(-15,15), type="n", xaxt="n", xlab="", ylab="")
abline(h=0, col='gray50')
for (i in 1:length(Chosen))	{
	vioplot(Mat[,Chosen[i]], at=i, add=TRUE, col=Cols[sum(Thresh>Pvals[Chosen[i]])])
}
axis(1, at=c(1:length(Chosen)), labels=c("", "", "", ""), tck=-0.02)
ToMark <- 1:length(Chosen)
Names <- c("mass", "wing:mass", "tarsus:mass", "tail:mass")
silent <- sapply(1:length(ToMark), function(x) text(ToMark[x], par("usr")[3]+Buffer, Names[x], adj = 1, srt = 35, xpd=NA, cex=1))
legend("topleft", bty="n", pch=21, pt.bg=Cols, legend=c("p > 0.1", "p < 0.1", "p < 0.05"), pt.cex=2)

Chosen <- c("b_rachisWL", "b_barbDR", "b_barbLR", "b_barbAng", "b_barbL_lead")
plot(1, 1, xlim=c(0,length(Chosen)+1), ylim=c(-15,15), type="n", xaxt="n", xlab="", ylab="")
abline(h=0, col='gray50')
for (i in 1:length(Chosen))	{
	vioplot(Mat[,Chosen[i]], at=i, add=TRUE, col=Cols[sum(Thresh>Pvals[Chosen[i]])])
}
axis(1, at=c(1:length(Chosen)), labels=c("","","","",""), tck=-0.01)
ToMark <- 1:length(Chosen)
Names <- gsub("b\\_", "", Chosen)
Names <- c("rachis width:length", "barb density t:l", "barb length t:l", "barb angle t:l", "barb length:rachis (l)")
silent <- sapply(1:length(ToMark), function(x) text(ToMark[x], par("usr")[3]+Buffer, Names[x], adj = 1, srt = 35, xpd=NA, cex=1))

dev.off()
