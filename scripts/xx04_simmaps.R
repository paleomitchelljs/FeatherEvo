library(phytools)
library(OUwie)

#setwd("~/Dropbox/Research/FeatherEvolution/scripts")
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\scripts")
source("xx00_readFormatDat.R")

#### Use simmap to map regimes onto tree, save simmap tree
Fly <- setNames(ModelDat$flightless, ModelDat$species)

Model <- matrix(c(0,0,1,0),2)
rownames(Model) <- colnames(Model) <- c("0", "1")

#### Load simmap tree, use OUwie to fit BMS and various OU models
strees <- make.simmap(Tree, Fly, model = Model, pi = setNames(c(1, 0), c("0","1")), nsim=1e2, Q="mcmc")

#setwd("~/Dropbox/Research/FeatherEvolution/phy")
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\phy")
write.simmap(strees, file="simtrees.tre", version = 1.5, format = "nexus")
saveRDS(strees, "strees.rds") # kludge since write.simmap was doing something annoying


#obj <- densityMap(strees, fsize=c(0.5,1), leg.text="P(Volant)", plot=FALSE, outline=TRUE)
#n<-length(obj$cols)
#obj$cols[1:n]<-grey(0:(n-1)/(n-1))
#pdf("densitymap.pdf")
#plot(obj, outline=TRUE, fsize=c(0.5, 1))
#dev.off()



library(phytools)
library(OUwie)

setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution")
source("xx00_readFormatDat.R")

#### Use simmap to map regimes onto tree, save simmap tree
Fly <- setNames(ModelDat$flightless, ModelDat$species)

#### Load simmap tree, use OUwie to fit BMS and various OU models
strees <- make.simmap(Tree, Fly, model = "ARD", pi = setNames(c(1, 0), c("0","1")), nsim=1e2, Q="mcmc")

setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\phy")
pdf("reversible_flightless.pdf")
plot(strees[[1]])
dev.off()
