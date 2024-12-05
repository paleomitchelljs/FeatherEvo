###
library(phytools)
library(brms)
library(bayesplot)setwd("~/Dropbox/Research/FeatherEvolution")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution")

Dat <- read.csv("primary.csv")
Dat_sp <- Dat[,1]
Dat_sp <- gsub(" ", "\\_", Dat_sp)
Dat_gen <- sapply(Dat[,1], function(x) unlist(strsplit(x, " "))[1])

setwd("~/Dropbox/Research/FeatherEvolution/phy")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\phy")
tree <- read.nexus("prum_etal_2015_vegavis.tre")

StripTips_sp <- sapply(tree$tip.label, function(x) paste(unlist(strsplit(x, "\\_"))[4:5], collapse=" "))
StripTips_gen <- sapply(tree$tip.label, function(x) paste(unlist(strsplit(x, "\\_"))[4], collapse=" "))

Bind <- sapply(StripTips_gen, function(x) ifelse(x%in%Dat_gen, Dat_sp[which(Dat_gen==x)], 0))

setdiff(Dat_sp, Bind)

#write.csv(cbind(StripTips_sp,Bind), "PrumSpecies.csv")
Guide <- read.csv('PrumSpecies_hand.csv')

tree_re <- tree

ToDrop <- c()
for (i in 1:Ntip(tree_re))	{
	Row <- which(Guide[,1] == tree_re$tip.label[i])
	if (Guide[Row,"Bind"] == 0)	{
		ToDrop <- c(ToDrop, i)
	}
	else	{
		tree_re$tip.label[i] <- Guide[Row,"Bind"]
	}
}

tree_re <- drop.tip(tree_re, ToDrop)

# read the gruiformes tree from xx00_trees.R file
gruTree2 <- read.tree("crane_tree_renamed_alda.tre")

# get length of gruiforme tree
gruH <- max(nodeHeights(gruTree2))   

# get the Rallus tip & it's length
bindTip <- grep("Rallus", tree_re$tip.label)
bindLength <- tree_re$edge.length[which(tree_re$edge[,2]==bindTip)]

# Bind two trees together
btree <- bind.tree(tree_re, gruTree2, where=bindTip, position=gruH)
btree <- drop.tip(btree, grep("(=Gallirallus)", btree$tip.label))

write.tree(btree, "wholetree.tre")

# Still to fix:
setdiff(Dat_sp, btree$tip.label)

# Have to hand-add taxa past here
# Prum backbone: https://doi.org/10.1038/nature15697
# Rails: DOI: 10.1016/j.ympev.2021.107106
# Penguins: https://doi.org/10.1073/pnas.2006659117 (Aptenodytes:22,(Spheniscus:15,Eudyptes:15):7);
# Strix nebulosa sister to Strix varia at 15.8Ma from Jetz et al 2012
# Xenicus lyalli 32.9 Ma sister based on spreadsheet
# Anas chlorotis + aucklandii 0.73 Ma from Mitchell (2013)
# Pinguinus based on Smith & Clarke 2015
# Mergus-Anas split 19.7Ma based on Fulton et al 2012 Proc B
# Anas - Tachyeres split 11.9 based on Fulton 2012
# Nestor-Strigops from Mitchell 2016 Ancient mito
# Podiceps - Rollandia 7.79Ma split from Ogawa et al 2015
# (Podiceps,Rollandia) - Podilymbus 28.5Ma split from Jetz et al 2012
# Podilymbus_gigas/podiceps split set to 0.959999Ma due to ambiguity, older than Tac split (discussion between Saitta and M van Tuinen)
# Livezey 1989 suggested M australis diverged right after Lophodytes, so 5.7Ma was used for Mergus split
#btree <- read.tree("fulltree.tre")
#plot(btree, cex=0.5)
#setdiff(Dat_sp, btree$tip.label)