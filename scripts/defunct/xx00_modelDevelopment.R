library(phytools)
library(brms)
library(bayesplot)
options(mc.cores = parallel::detectCores())
setwd("C:\\Users\\jonsm\\Dropbox\\Research\\FeatherEvolution")

Dat <- read.csv("primary.csv")
Pairs <- read.csv("sisterpairs.csv")

Time <- Dat[,"Age"]
Flight <- Dat[,2]
Feathers <- Dat[,9:21]
Macro <- Dat[,c("Mass","WingL","TarsusL")]
BaseCors <- cor(Flight,Feathers, use="complete")

##### Sanity check -- sister species pair correlations
sisCors <- c()
for (j in 1:ncol(Feathers))	{

	Trait <- j
	Pts <- matrix(0, ncol=2, nrow=max(Pairs$pair, na.rm=T))

	#plot(1, 1, type="n", xlab="", ylab="", ylim=c(0, max(Feathers[,Trait], na.rm=T)), xlim=c(0,1), main=colnames(Feathers)[j])
	plot(1, 1, type="n", xlab="age of split", ylab="trait difference", xlim=c(1e4,1e7), ylim=c(-1*mean(Feathers[,Trait]), mean(Feathers[,Trait])), log="x")
	Npairs <- max(Pairs$pair, na.rm=T)
	for (i in 1:Npairs)	{
		Choose <- which(Pairs$pair == i)
		Measures <- Feathers[Choose,Trait]
		DivT <- max(Time[Choose], na.rm=T) 
		if (length(DivT) == 1)	{ 
			points(DivT, diff(Measures))
			Pts[i,] <- c(DivT, diff(Measures))
		}
		else	{
			Pts[i,] <- c(NA, diff(Measures))
		}
	#	segments(Flight[Choose[1]], Measures[1], Flight[Choose[2]], Measures[2], lwd=2)
	}
	sisCors[j] <- cor(Pts[-7,1], Pts[-7,2], use="complete")
}
#################

plot(abs(sisCors), abs(BaseCors[1,]), xlim=c(0, 1), ylim=c(0, 1))
text(abs(sisCors), abs(BaseCors[1,]), 1:14)
abline(v=0.25, lty=2)
abline(h=0.25, lty=2)

# correlated traits from naiive analyses
Traits <- c(1,2,3,7,8,10,11)
Dat2 <- Feathers[,Traits]
Scores <- princomp(Dat2, cor=FALSE, scores = TRUE )

fScale <- function(x, fudge = 1e-5)	{
	x <- x + fudge
	x <- x / mean(x)
	return(log(x))
}

ModelDat <- data.frame(
	species = Dat[,1], 
	clade = Dat[,4],
	age = Dat[,22],
	swim = factor(Dat[,5], ordered=TRUE), 
	mass = log(Dat[,6]), 
	flight = factor(Flight, ordered=FALSE), 
	Iwidth = Dat2[,1] - median(Dat2[,1]), 
	Ilength = Dat2[,2] - median(Dat2[,2]), 
	Mwidth = log(Dat2[,3]), 
	Mlength = log(Dat2[,4]),
	bWidth = Dat2[,5] - median(Dat2[,5]), 
	DbLength = Dat2[,6] - median(Dat2[,6]), 
	DbDense = Dat2[,7] - median(Dat2[,7]), 
	PbDense = Dat2[,8] - median(Dat2[,8])
	)

Flighted <- which(as.character(ModelDat$flight) == "1")
Flightless <- which(as.character(ModelDat$flight) == "0")

FlightedByClade <- tapply(ModelDat[Flighted,"DbDense"], ModelDat[Flighted, "clade"], mean, na.omit=T)
FlightlessByClade <- tapply(ModelDat[Flightless,"DbDense"], ModelDat[Flightless, "clade"], mean, na.omit=T)

Range <- range(pretty(range(c(FlightedByClade,FlightlessByClade))))
par(mar=c(4,4,2,1), las=1, mgp=c(2, 0.5, 0), tck=-0.01)
plot(FlightedByClade[names(FlightlessByClade)], FlightlessByClade, xlim=Range, ylim=Range, pch=16, xlab="distal barb density (flightless)", ylab="distal barb density (volant)", main="mean by clade")
abline(0, 1, lwd=2, lty=2)

# Test model construction on Dbdense (distal barb density) as it has highest absolute 'naive' correlation. Work out exploratory analyses there, then see if they hold throughout. This prevents "fishing" in all of the variables, while still allowing exploration and proper model tuning.
# https://rpubs.com/jwesner/gamma_glm
# https://magesblog.com/post/2018-08-02-use-domain-knowledge-to-review-prior-predictive-distributions/

Form <- bf(
		DbDense ~ mo(swim) + mass + Mwidth + flight + (1 | clade)#,
#		sigma ~ flight + (1 | clade)
	)


#### Check prior behavior
Priors <- get_prior(Form, data=ModelDat)

#mcmc_areas(as.array(prior_model_simple)) 
#conditions <- data.frame(clade = unique(ModelDat$clade))
#rownames(conditions) <- unique(ModelDat$clade)

#me_prior <- marginal_effects(
#  prior_model_simple, conditions = conditions, 
#  re_formula = NULL, method = "predict"
#)
#p1 <- plot(me_prior, ncol = 3, points = TRUE, plot = FALSE)


#### Once satisfied with prior, run model
model_simple <- brm(
					Form,
					data = ModelDat,
					family = gaussian(), #lognormal(link = "identity", link_sigma = "log"), # barb density can't be negative, so can't use gaussian response. gamma() would be for heteroskedatic
					prior = c(
						prior(normal(0, 3), "b"),
    					prior(normal(0, 3), "Intercept"),
						prior(student_t(3, 0, 0.5), "sd"),
						prior(student_t(3, 0, 5), "sigma")
  						),
					cores = 4,
					iter = 5e3,
					control = list(adapt_delta = 0.9, max_treedepth = 12)
					)


pp_check(model_simple)

mcmc_plot(model_simple)

me_post_flight <- conditional_effects(model_simple, effects="flight")

conditions <- data.frame(clade = unique(ModelDat$clade))
rownames(conditions) <- unique(ModelDat$clade)

me_post_clade <- marginal_effects(
  model_simple, conditions = conditions, 
  re_formula = NULL, method = "predict"
)
p1 <- plot(me_post_clade, ncol = 3, points = TRUE, plot = FALSE)


add_ic(model_simple) <- "loo"
summary(model_simple)


########## PCA 
ModelDat <- data.frame(
	species = Dat[,1], 
	clade = Dat[,4],
	age = Dat[,22],
	swim = factor(Dat[,5], ordered=TRUE), 
	mass = log(Dat[,6]), 
	flight = factor(Flight, ordered=FALSE), 
	pc1 = Scores$scores[,1],
	pc2 = Scores$scores[,2],
	pc3 = Scores$scores[,3],
	pc4 = Scores$scores[,4]
	)

Form <- bf(
		mvbind(pc1,pc2) ~ mo(swim) + mass + flight + (1 | clade),
		sigma ~ flight + (1 | clade)
	) + set_rescor(TRUE)

#### Check prior behavior
Priors <- get_prior(Form, data=ModelDat)

#### Once satisfied with prior, run model
model_simple <- brm(
					Form,
					data = ModelDat,
					family = gaussian(), #lognormal(link = "identity", link_sigma = "log"), # barb density can't be negative, so can't use gaussian response. gamma() would be for heteroskedatic
					cores = 4,
					iter = 5e3,
					control = list(adapt_delta = 0.9, max_treedepth = 12)
					)


pp_check(model_simple, resp="pc2")

mcmc_plot(model_simple)

me_post_flight <- conditional_effects(model_simple, effects="flight")

conditions <- data.frame(clade = unique(ModelDat$clade))
rownames(conditions) <- unique(ModelDat$clade)

me_post_clade <- marginal_effects(
  model_simple, conditions = conditions, 
  re_formula = NULL, method = "predict"
)
p1 <- plot(me_post_clade, ncol = 3, points = TRUE, plot = FALSE)
