library(phytools)
library(brms)
library(bayesplot)
options(mc.cores = parallel::detectCores())

setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution")
Dat <- read.csv("primary.csv")
Pairs <- read.csv("sisterpairs.csv")

Time <- Dat[,"Max..time.flightless...yr."]
Flight <- Dat[,2]
Feathers <- Dat[,9:21]

# correlated traits from naiive analyses
Traits <- c(1,2,3,4,7,8,10,11)
Dat2 <- Feathers[,Traits]
Scores <- princomp(Dat2, cor=FALSE, scores = TRUE )

########## PCA 
ModelDat <- data.frame(
	species = Dat[,1], 
	clade = Dat[,4],

	age = Dat[,22],
	swim = factor(Dat[,5], ordered=TRUE), 
	mass = log(Dat[,6]), 
	flight = factor(Flight, ordered=FALSE), 
	pc1 = Scores$scores[,1]/100,
	pc2 = Scores$scores[,2]/30,
	pc3 = Scores$scores[,3]/10,
	pc4 = Scores$scores[,4]/10
	)

Form <- bf(
		pc1 ~ mo(swim) + mass + flight + (1 | clade)#,
#		sigma ~ flight + (1 | clade)
	) #+ set_rescor(TRUE)

#### Check prior behavior
Priors <- get_prior(Form, data=ModelDat)

#### Once satisfied with prior, run model
model_pcBirds <- brm(
					Form,
					data = ModelDat,
					family = gaussian(), #lognormal(link = "identity", link_sigma = "log"), # barb density can't be negative, so can't use gaussian response. gamma() would be for heteroskedatic
					prior = c(
						prior(normal(0, 0.5), "b"),
    					prior(normal(0, 0.5), "Intercept"),
						prior(student_t(3, 0, 1), "sd"),
						prior(student_t(3, 0, 1), "sigma")
  						),
					cores = 4,
					iter = 5e3,
					control = list(adapt_delta = 0.9, max_treedepth = 12)
					)

pdf("allbirds_clade.pdf")
plot(conditional_effects(model_pcBirds), points = TRUE) 
dev.off()

pp_check(model_pcBirds)

pdf("allbirds_pc1_mcmc.pdf")
mcmc_plot(model_pcBirds)
dev.off()

me_post_flight <- conditional_effects(model_pcBirds, effects="flight")

conditions <- data.frame(clade = unique(ModelDat$clade))
rownames(conditions) <- unique(ModelDat$clade)

me_post_clade <- marginal_effects(
  model_pcBirds, conditions = conditions, 
  re_formula = NULL, method = "predict"
)
p1 <- plot(me_post_clade, ncol = 3, points = TRUE, plot = FALSE)


add_ic(model_pcBirds) <- "loo"
summary(model_pcBirds)

######### Gruiformes
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\phy")
gruTree2 <- read.tree("crane_tree_renamed_alda.tre")
Gruiformes <- ModelDat[which(ModelDat[,"clade"] == "Gruiformes"),]
Rename <- sapply(gruTree2$tip.label, function(x) paste(unlist(strsplit(x, "_"))[c(1,2)], collapse=" "))
gruTree2$tip.label <- Rename

A <- ape::vcv.phylo(gruTree2)

Form <- bf(
		pc1 ~ mass + flight + (1 | gr(species, cov = A))
	) #+ set_rescor(TRUE)

bprior <- prior(normal(0,0.5), class = b) + 
          prior(normal(0,0.01), class = sd) + 
          prior(normal(0, 0.5), class = Intercept) +
          prior(student_t(3,0,2), class = sigma)

#### Once satisfied with prior, run model
model_gru <- brm(
					Form,
					data = Gruiformes,
					family = gaussian(), #lognormal(link = "identity", link_sigma = "log"), # barb density can't be negative, so can't use gaussian response. gamma() would be for heteroskedatic
					prior = bprior,
					data2 = list(A = A),
					cores = 4,
					iter = 1e4,
					control = list(adapt_delta = 0.999, max_treedepth = 12)
					)

pdf("gruiformes_cond_nosig_pc1.pdf")
plot(conditional_effects(model_gru), points = TRUE) 
dev.off()

pp_check(model_gru)

pdf("gruiformes_mcmc_nosig_pc1.pdf")
mcmc_plot(model_gru)
dev.off()


Form <- bf(
		pc2 ~ mass + flight + (1 | gr(species, cov = A)),
		sigma ~ flight
	) #+ set_rescor(TRUE)

bpriorsig <- prior(normal(0,0.5), class = b) + 
          prior(normal(0,0.01), class = sd) + 
          prior(normal(0, 0.5), class = Intercept)

#### Once satisfied with prior, run model
model_gru_sig <- brm(
					Form,
					data = Gruiformes,
					family = gaussian(), #lognormal(link = "identity", link_sigma = "log"), # barb density can't be negative, so can't use gaussian response. gamma() would be for heteroskedatic
					prior = bpriorsig,
					data2 = list(A = A),
					cores = 4,
					iter = 1e4,
					control = list(adapt_delta = 0.999, max_treedepth = 12)
					)

mcmc_plot(model_gru_sig)
