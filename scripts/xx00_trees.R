library(phytools)
library(brms)
library(bayesplot)
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution")

Dat <- read.csv("primary.csv")
Pairs <- read.csv("sisterpairs.csv")
#Time <- Dat[,"Max..time.flightless...yr."]

Gruiformes <- Dat[which(Dat[,"Phylogenetic..group"] == "Gruiformes"),]

# Read in tree
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\phy")
cranes <- read.tree("garcia_gruiformes.newick")

# Reformat names
Rename <- sapply(cranes$tip.label, function(x) paste(unlist(strsplit(x, "_"))[c(1,2)], collapse=" "))
cranes$tip.label <- Rename

# Force the species names in the data file to match tree
ToFix <- setdiff(Gruiformes[,1], Rename)
ToFix

#Gruiformes[which(Gruiformes[,1] == ToFix[1]),1] <- Rename[grep("cuvieri", Rename)]
#Gruiformes[which(Gruiformes[,1] == ToFix[2]),1] <- Rename[grep("mortierii", Rename)]
#Gruiformes[which(Gruiformes[,1] == ToFix[3]),1] <- "Dryolimnas aldabranus"
# Dryolimnas aldabranus needs to be sister to Dryolimnas cuvieri with divergence time of 0.13 Ma
#Gruiformes[which(Gruiformes[,1] == ToFix[4]),1] <- Rename[grep("philippensis", Rename)]
#Gruiformes[which(Gruiformes[,1] == ToFix[5]),1] <- Rename[grep("sylvestris", Rename)]
#Gruiformes[which(Gruiformes[,1] == ToFix[6]),1] <- Rename[grep("pusilla", Rename)]
#Gruiformes[which(Gruiformes[,1] == ToFix[7]),1] <- Rename[grep("tabuensis", Rename)]
#Gruiformes[which(Gruiformes[,1] == ToFix[8]),1] <- Rename[grep("chloropus", Rename)]
#Gruiformes[which(Gruiformes[,1] == ToFix[9]),1] <- Rename[grep("ventralis", Rename)]

# Prune tree
gruTree <- drop.tip(cranes, setdiff(cranes$tip.label, Gruiformes[,1]))
write.tree(gruTree, "crane_tree_renamed.tre")
gruTree2 <- read.tree("crane_tree_renamed_alda.tre")


