library(phytools)
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\scripts")
source("xx00_readFormatDat.R")

library(brms)
library(bayesplot)
options(mc.cores = parallel::detectCores())


#### TEST DENSITY ASYMMETRY
BF <- bf( barbuleDR ~ flightless + (1 | gr(species, cov = A)), sigma ~ flightless + (1 | gr(species, cov = A)) )

prior <- get_prior(BF, data=RawDat, data2 = list(A = A))
prior$prior[which(prior$class == "Intercept")] <- "student_t(3, 0, 5)"
prior$prior[which(prior$class == "b")] <- "student_t(3, 0, 1)"
prior$prior[which(prior$group == "species")] <- "student_t(3, 0, 0.001)"

ITER <- 1e5
brblFit <- brm(BF, 
			family =  gaussian(), 
			data=RawDat, 
			data2 = list(A = A), 
			prior=prior, 
			cores = 4,
			warmup = ITER / 2,
			save_pars = save_pars(all = TRUE),
			iter = ITER,
			control = list(adapt_delta = 0.999, max_treedepth = 12)
			)

hyp <- c("exp(sigma_Intercept) = 0",
         "exp(sigma_Intercept + sigma_flightless) = 0")
hypothesis(brblFit, hyp)

#### TEST LENGTH ASYMMETRY
BF <- bf( barbLR ~ flightless + (1 | gr(species, cov = A)), sigma ~ flightless)

prior <- get_prior(BF, data=RawDat, data2 = list(A = A))
prior$prior[which(prior$class == "Intercept")] <- "student_t(3, 0, 5)"
prior$prior[which(prior$class == "b")] <- "student_t(3, 0, 1)"
prior$prior[which(prior$group == "species")] <- "student_t(3, 0, 0.001)"

ITER <- 1e5
brbFit <- brm(BF, 
			family =  gaussian(), 
			data=RawDat, 
			data2 = list(A = A), 
			prior=prior, 
			cores = 4,
			warmup = ITER / 2,
			save_pars = save_pars(all = TRUE),
			iter = ITER,
			control = list(adapt_delta = 0.999, max_treedepth = 12)
			)

hyp <- c("exp(sigma_Intercept) = 0",
         "exp(sigma_Intercept + sigma_flightless) = 0")
hypothesis(brbFit, hyp)

