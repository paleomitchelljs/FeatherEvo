library(phytools)
Store <- getwd()
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\phy")
setwd("~/Dropbox/Research/FeatherEvolution/phy")

Tree <- read.tree("fulltree.tre")
Tree$tip.label <- gsub("\\_", " ", Tree$tip.label)
Tree <- drop.tip(Tree, grep("Casuar", Tree$tip.label))
A <- ape::vcv.phylo(Tree)

degreeToRad <- pi / 180
setwd("~/Dropbox/Research/FeatherEvolution/data")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\data")

Pairs <- read.csv("sisterpairs.csv")

RawDat <- read.csv("primary.csv")
# Drop Casuarius
RawDat <- RawDat[-29,]
# remove the _m tags as all used data are from middle
colnames(RawDat) <- gsub("\\_m", "", colnames(RawDat))

RawDat$wingR <- log(RawDat$wing) - log(RawDat$mass)
RawDat$tarsusR <- log(RawDat$tarsus) - log(RawDat$mass)
RawDat$tailR <- log(RawDat$tail) - log(RawDat$mass)

RawDat$rachisWL <- log(RawDat$rachis_width) - log(RawDat$rachis_length)
RawDat$barbDR <- log(RawDat$barb_density_t) - log(RawDat$barb_density_l)
RawDat$barbLR <- log(RawDat$barb_length_t) - log(RawDat$barb_length_l)

RawDat$barbAng <- RawDat$barb_angle_t - RawDat$barb_angle_l # Feo et al 2015

RawDat$barbL_lead <- log(RawDat$barb_length_l) - log(RawDat$rachis_length) #
RawDat$barbL_trail <- log(RawDat$barb_length_t) - log(RawDat$rachis_length) #
RawDat$barbuleL_distrail <- RawDat$barbule_length_t / RawDat$rachis_length #

RawDat$mass <- log(RawDat$mass)
RawDat$wing <- log(RawDat$wing)
RawDat$tarsus <- log(RawDat$tarsus)
RawDat$tail <- log(RawDat$tail)

# flip flightless/volant for ease of interpreting graphs
RawDat$flightless <- rep(0, nrow(RawDat))
RawDat$flightless[which(RawDat$flight == 0)] <- 1

#colnames(RawDat)[1] <- "species"

# Set Mesitornis as flightless since it may as well be
Fly <- RawDat$flight
Feathers <- RawDat[,c("barb_angle_l", "barb_angle_t", "rachis_width", "rachis_length", "barb_length_l", "barb_length_t", "barb_density_l", "barb_density_t", "barbule_length_t", "barbule_density_distal", "barbule_density_proximal")]
rownames(Feathers) <- RawDat$species
BaseCors <- cor(Fly, RawDat[,11:32], use="complete")
colnames(BaseCors)[order(abs(BaseCors))]

Macro <- RawDat[,c("mass", "wingR", "tarsusR", "tailR")]

# Phylo PCA? Set barb angle to 0 instead of NA for Casuaris
Scores <- phyl.pca( Tree, Feathers, method = "lambda" , mode="cov")

# Adjust A matrix -- UNNECESSARY PER https://discourse.mc-stan.org/t/phylogenetic-signal-in-brms/16457/3
#DiagA <- diag(A)
#A <- A*Scores$lambda
#diag(A) <- DiagA

centerscale <- function(x, Scale=F)	{
	out <- (x - mean(x, na.rm=T, finite=T))
	if (Scale)	{
		out <- out / sd(x, na.rm=T)
	}
	return(out)
}
########## 
ModelDat <- RawDat
ModelDat[,1] <- RawDat[,2]
ModelDat <- ModelDat[,-ncol(ModelDat)]
ModelDat[,2] <- RawDat$flightless
colnames(ModelDat)[1] <- "species"
colnames(ModelDat)[2] <- "flightless"
ModelDat[,7:ncol(ModelDat)] <- apply(ModelDat[,7:ncol(ModelDat)], 2, centerscale)

ModelDat <- cbind(ModelDat, Scores$S)

setwd(Store)
######## Create a derived-trait ModelDat variant that looks at log(wing/body) etc ratios and includes asymmetry. Additionally, check pairwise correlations with mass to look for feather traits to adjust
######## Install JAGS
######## Run the Fuentes model (https://onlinelibrary.wiley.com/doi/am-pdf/10.1111/evo.13899)
######## Need to modify Fuentes to allow separate intercepts, too

