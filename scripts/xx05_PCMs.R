library(phytools)
library(OUwie)

setwd("~/Dropbox/Research/FeatherEvolution/scripts")
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\scripts")
source("xx00_readFormatDat.R")

#### Load simmap tree, 
setwd("~/Dropbox/Research/FeatherEvolution/output")
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\output")
#strees <- read.simmap("simtrees.tre", format="phylip")
strees <- readRDS("strees.rds")

dwtrees <- sapply(strees, drop.tip.simmap, tip="Casuarius unappendiculatus", simplify=F)
class(dwtrees) <- "multiPhylo"

Rows <- which(ModelDat$species != "Casuarius unappendiculatus")
Traits <- ModelDat[Rows,c("mass", "wingR", "tarsusR", "rachisWL", "barbDR", "barbLR", "barbAng", "barbL_lead")]
rownames(Traits) <- ModelDat$species

lessMat <- matrix(0, nrow=100, ncol=ncol(Traits))
volMat <- matrix(0, nrow=100, ncol=ncol(Traits))
oneMat  <- matrix(0, nrow=100, ncol=ncol(Traits))
aiccMat  <- matrix(0, nrow=100, ncol=ncol(Traits))
oaiccMat  <- matrix(0, nrow=100, ncol=ncol(Traits))

#use OUwie to fit BMS and BM
for (i in 1:length(strees))	{
	for (j in 1:ncol(Traits))	{
		useTree <- dwtrees[[i]]
		useDat <- data.frame(rownames(Traits), RawDat$flightless, Traits[,j], rep(MErr[colnames(Traits)[j]], nrow(Traits)))

		BMs <- OUwie(useTree, useDat, model="BMS", simmap.tree=TRUE, root.station = FALSE, mserr="known")
		lessMat[i,j] <- BMs$solution[2,"1"]
		volMat[i,j] <- BMs$solution[2,"0"]
		aiccMat[i,j] <- BMs$AICc

		# don't need to do the one-rate model for every simmap tree, just for every trait.
		if (i == 1)	{
			BM <- OUwie(useTree, useDat, model="BM1", simmap.tree=TRUE, root.station = FALSE, mserr="known")
			oneMat[i,j] <- BM$solution[2,1]
			oaiccMat[i,j] <- BM$AICc
		}
		else if (i != 1)	{
			oneMat[i,j] <- oneMat[1, j]
			oaiccMat[i, j] <- oaiccMat[1, j]
		}
	}
	cat("\n\n Done with tree ", i, "\n\n")
}

setwd("~/Dropbox/Research/FeatherEvolution/output")
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\output")
write.csv(lessMat, "lessMat.csv")
write.csv(volMat, "volMat.csv")
write.csv(oneMat, "oneMat.csv") 
write.csv(aiccMat, "aiccMat.csv") 
write.csv(oaiccMat, "oaiccMat.csv") 


##############################################################################
##############################################################################
##############################################################################
##############################################################################
setwd("~/Dropbox/Research/FeatherEvolution/output")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\output")

raw_lessMat <- matrix(0, nrow=100, ncol=ncol(rawTraits))
raw_volMat <- matrix(0, nrow=100, ncol=ncol(rawTraits))
raw_oneMat <- matrix(0, nrow=100, ncol=ncol(rawTraits))
raw_aiccMat <- matrix(0, nrow=100, ncol=ncol(rawTraits))
raw_oaiccMat <- matrix(0, nrow=100, ncol=ncol(rawTraits))
for (i in 1:length(strees))	{
	for (j in 1:ncol(rawTraits))	{
		useTree <- strees[[i]]
		useDat <- RawDat[Rows,c(1,3,6+j)]
		if (colnames(Traits)[j] == "dW")	{
			useTree <- dwtrees[[i]]
			useDat <- RawDat[ModelDat[Rows,1]%in%useTree$tip.label,c(1,2,6+j)]
		}

		rBMs <- OUwie(useTree, useDat, model="BMS", simmap.tree=TRUE, root.station = FALSE, ub=Inf)
		raw_lessMat[i,j] <- rBMs$solution[2,1]
		raw_volMat[i,j] <- rBMs$solution[2,2]
		raw_aiccMat[i,j] <- rBMs$AICc

		rBM <- OUwie(useTree, useDat, model="BM1", simmap.tree=TRUE, root.station = FALSE, ub=Inf)
		raw_oneMat[i,j] <- rBM$solution[2,1]
		raw_oaiccMat[i,j] <- rBM$AICc

	}
	cat("\n\n Done with tree ", i, "\n\n")
}

write.csv(raw_lessMat, "raw_lessMat.csv")
write.csv(raw_volMat, "raw_volMat.csv")
write.csv(raw_oneMat, "raw_oneMat.csv")
write.csv(raw_aiccMat, "raw_aiccMat.csv")
write.csv(raw_oaiccMat, "raw_oaiccMat.csv")
