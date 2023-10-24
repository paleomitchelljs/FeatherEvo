library(phytools)
library(brms)
library(bayesplot)
options(mc.cores = parallel::detectCores())

setwd("~/Dropbox/Research/FeatherEvolution/scripts")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\scripts")
source("xx00_readFormatDat.R")

######### All birds
## Logistic regression

################
################

Form <- bf( flightless ~ mass + wing + tarsus + rachis_width + rachis_length + barbDR + barbLR + barbAng + barbuleL_distrail + barbule_density_distal, family = bernoulli(link = "logit") )

Priors <- get_prior(Form, data=RawDat)
prior <- Priors

prior$prior[which(prior$class == "Intercept")] <- "student_t(3, 0, 2)"
prior$prior[which(prior$class == "b")] <- "student_t(3, 0, 1)"

model_logreg <- brm(
					Form,
					data = RawDat,
					family = bernoulli(link = "logit"), #lognormal(link = "identity", link_sigma = "log"), # barb density can't be negative, so can't use gaussian response. gamma() would be for heteroskedatic
					prior = prior,
					cores = 4,
					warmup = 5e4,
					iter = 1e5,
					control = list(adapt_delta = 0.9, max_treedepth = 12)
					)

setwd("~/Dropbox/Research/FeatherEvolution/output")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\output")
saveRDS(model_logreg, "model_logreg.rds")
pp_check(model_logreg, ndraws=100)
pp_check(model_logreg, ndraws=100, type = "scatter_avg")

bayes_R2(model_logreg)
summary(model_logreg)

par(mfrow=c(1,3))
mcmc_plot(model_logreg, variable = c("b_mass", "b_wing", "b_tarsus"))
mcmc_plot(model_logreg, variable = c("b_rachis_width", "b_rachis_length"))
mcmc_plot(model_logreg, variable = c("b_barbDR", "b_barbLR", "b_barbAng"))
