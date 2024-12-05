library(phytools)
library(brms)
library(bayesplot)
options(mc.cores = parallel::detectCores())
setwd("~/Dropbox/Research/FeatherEvolution/scripts")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\scripts")
source("xx00_readFormatDat.R")

setwd("~/Dropbox/Research/FeatherEvolution/output")
#setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution\\output")

###################################
#### Logistic regression Phylo ####
###################################

##################################################
## RAN 12 NOV FOR TABLE 2
Form <- bf( flightless ~ mass + wing + tarsus + tail + barb_density_l + barb_density_t + rachis_width + rachis_length + barb_length_l + barb_length_t + barbule_density_distal + barb_angle_t + barb_angle_l + (1 | gr(species, cov = A)), family = bernoulli() )

Priors <- get_prior(Form, data=ModelDat, data2 = list(A = A))
prior <- Priors

prior$prior[which(prior$class == "Intercept")] <- "student_t(3, 0, 3)"
prior$prior[which(prior$class == "b")] <- "student_t(3, 0, 2)"
prior$prior[which(prior$group == "species")] <- "student_t(3, 0, 1e-1)"

ITER <- 1e4
model_logreg_phy_len <- brm(
					Form,
					data = ModelDat,
					data2 = list(A = A),
					family = bernoulli(link = "logit"), #lognormal(link = "identity", link_sigma = "log"), # barb density can't be negative, so can't use gaussian response. gamma() would be for heteroskedatic
					prior = prior,
					cores = 4,
					warmup = ITER / 2,
					iter = ITER,
					save_pars = save_pars(all = TRUE),
					control = list(adapt_delta = 0.99, max_treedepth = 12)
					)
saveRDS(model_logreg_phy_len, "model_logreg_phy_len.rds")

write.csv(summary(model_logreg_phy_len)$fixed, "table2.csv")
#####################################################

Form <- brms:::bf( flightless ~ mass * (wing + tail + tarsus) + rachis_width * rachis_length + barb_density_l * barb_density_t + barb_length_l * barb_length_t + barb_angle_t * barb_angle_l + barbule_density_distal + (1 | gr(species, cov = A)), family = bernoulli() )

Priors <- get_prior(Form, data=ModelDat, data2 = list(A = A))
prior <- Priors

prior$prior[which(prior$class == "Intercept")] <- "student_t(3, 0, 2)"
prior$prior[which(prior$class == "b")] <- "student_t(3, 0, 1)"
prior$prior[which(prior$class == "sd")] <- "student_t(3, 0, 1.5)"
prior$prior[which(prior$group == "species")] <- "student_t(2, 0, 0.5)"
prior$prior[which(prior$coef == "mass")] <- "student_t(3, 2, 1)"
prior$prior[which(prior$coef == "wing")] <- "student_t(3, -2, 1)"
prior <- prior[-1,]

ITER <- 1e5
model_logreg_phy_len_red <- brm(
					Form,
					data = ModelDat,
					data2 = list(A = A),
					family = bernoulli(link = "logit"), #lognormal(link = "identity", link_sigma = "log"), # barb density can't be negative, so can't use gaussian response. gamma() would be for heteroskedatic
					prior = prior,
					cores = 4,
					warmup = ITER / 2,
					iter = ITER,
					save_pars = save_pars(all = TRUE),
					control = list(adapt_delta = 0.99, max_treedepth = 12)
					)
saveRDS(model_logreg_phy_len_red, "model_logreg_phy_len_red.rds")

#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################


Form <- bf( flightless ~ mass + wing + tarsus + tail + PC1 + PC2 + PC3 + PC4 + (1 | gr(species, cov = A)), family = bernoulli(link = "logit") )

Priors <- get_prior(Form, data=ModelDat, data2 = list(A = A))
prior <- Priors

prior$prior[which(prior$class == "Intercept")] <- "student_t(3, 0, 2)"
prior$prior[which(prior$class == "b")] <- "student_t(3, 0, 1)"
prior$prior[which(prior$group == "species")] <- "student_t(3, 0, 0.01)"


ITER <- 1e5
model_logreg_phy <- brm(
					Form,
					data = ModelDat,
					data2 = list(A = A),
					family = bernoulli(link = "logit"), #lognormal(link = "identity", link_sigma = "log"), # barb density can't be negative, so can't use gaussian response. gamma() would be for heteroskedatic
					prior = prior,
					cores = 4,
					warmup = ITER / 2,
					save_pars = save_pars(all = TRUE),
					iter = ITER,
					control = list(adapt_delta = 0.99, max_treedepth = 12)
					)
saveRDS(model_logreg_phy, "model_logreg_phy.rds")

#summary(model_logreg_phy)
#pdf("model_PC.pdf", height = 4, width = 8)
#par(mfrow=c(1,2))
#mcmc_plot(model_logreg_phy, variable = c("b_mass", "b_wing", "b_tarsus"))
#mcmc_plot(model_logreg_phy, variable = c("b_pc1", "b_pc2", "b_pc3"))
#dev.off()

#mcmc_plot(ranef(model_logreg_phy)$clade)

##################################################
##################################################
### FULL TESTING Added 5 June 2023 to address all of Evan's hypotheses at once
### Modified againt 10 Oct
##################################################
##################################################
# 5 June
# Form <- bf( flightless ~ mass * (wing + tarsus + tail) + rachis_middle_width * (barb_length_lm + barbule_length_distal_tm) + barbDR + barbLR + WL + mo(swim) + (1 | gr(species, cov = A)), family = bernoulli(link = "logit") )

Form <- bf( flightless ~ mass * (wing + tarsus) + rachis_width * rachis_length + barbDR + barbLR + barbAng + barbuleL_distrail + barbule_density_distal + (1 | gr(species, cov = A)), family = bernoulli(link = "logit") )

Priors <- get_prior(Form, data=ModelDat, data2 = list(A = A))
prior <- Priors

prior$prior[which(prior$class == "Intercept")] <- "student_t(3, 0, 2)"
prior$prior[which(prior$class == "b")] <- "student_t(3, 0, 1)"
prior$prior[which(prior$class == "sd")] <- "student_t(3, 0, 1.5)"
prior$prior[which(prior$group == "species")] <- "student_t(3, 0, 1)"

prior <- prior[-1,]

ITER <- 1e5
model_logreg_phy_nopc <- brm(
					Form,
					data = ModelDat,
					data2 = list(A = A),
					family = bernoulli(link = "logit"), #lognormal(link = "identity", link_sigma = "log"), # barb density can't be negative, so can't use gaussian response. gamma() would be for heteroskedatic
					prior = prior,
					cores = 4,
					warmup = ITER / 2,
					iter = ITER,
					save_pars = save_pars(all = TRUE),
					control = list(adapt_delta = 0.99, max_treedepth = 12)
					)
saveRDS(model_logreg_phy_nopc, "model_logreg_phy_FULL.rds")

#####################################################
#### Logistic regression Phylo Traits MACRO ONLY ####
#####################################################
Form <- bf( flightless ~ mass + wing + tarsus + tail + (1 | gr(species, cov = A)), family = bernoulli(link = "logit") )

Priors <- get_prior(Form, data=ModelDat, data2 = list(A = A))
prior <- Priors

prior$prior[which(prior$class == "Intercept")] <- "student_t(3, 0, 5)"
prior$prior[which(prior$class == "b")] <- "student_t(3, 0, 1)"
prior$prior[which(prior$group == "species")] <- "student_t(3, 0, 0.01)"

ITER <- 1e5
model_logreg_phy_macro <- brm(
					Form,
					data = ModelDat,
					data2 = list(A = A),
					family = bernoulli(link = "logit"), #lognormal(link = "identity", link_sigma = "log"), # barb density can't be negative, so can't use gaussian response. gamma() would be for heteroskedatic
					prior = prior,
					cores = 4,
					warmup = ITER / 2,
					iter = ITER,
					save_pars = save_pars(all = TRUE),					
					control = list(adapt_delta = 0.99, max_treedepth = 12)
					)
saveRDS(model_logreg_phy_macro, "model_logreg_phy_macro.rds")

#summary(model_logreg_phy_nopc)
#par(mfrow=c(1,2))
#mcmc_plot(model_logreg_phy_nopc, variable = c("b_mass", "b_wing", "b_tarsus"))
#dev.off()


Form <- brms:::bf( flightless ~ mass + wingR + tailR + tarsusR + rachisWL + barbDR + barbLR + barbAng + barbL_lead + (1 | gr(species, cov = A)), family = bernoulli() )

Priors <- get_prior(Form, data=ModelDat, data2 = list(A = A))
prior <- Priors

prior$prior[which(prior$class == "Intercept")] <- "student_t(3, 0, 2)"
prior$prior[which(prior$class == "b")] <- "student_t(3, 0, 1)"
prior$prior[which(prior$class == "sd")] <- "student_t(3, 0, 1.5)"
prior$prior[which(prior$group == "species")] <- "student_t(2, 0, 0.5)"
prior$prior[which(prior$coef == "mass")] <- "student_t(3, 2, 1)"
prior$prior[which(prior$coef == "wingR")] <- "student_t(3, -2, 1)"
prior <- prior[-1,]

ITER <- 1e5
model_logreg_phy_len_pcm_traits <- brm(
					Form,
					data = ModelDat,
					data2 = list(A = A),
					family = bernoulli(link = "logit"), #lognormal(link = "identity", link_sigma = "log"), # barb density can't be negative, so can't use gaussian response. gamma() would be for heteroskedatic
					prior = prior,
					cores = 4,
					warmup = ITER / 2,
					iter = ITER,
					save_pars = save_pars(all = TRUE),
					control = list(adapt_delta = 0.99, max_treedepth = 12)
					)
saveRDS(model_logreg_phy_len_pcm_traits, "model_logreg_phy_len_pcm_traits.rds")
