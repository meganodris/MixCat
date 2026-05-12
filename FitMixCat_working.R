rm(list=ls())
library(rstan)
library(bayesplot)
library(ggplot2)


source("C:/Users/megan/Documents/GitHub/MixCat/utils.R")
setwd("C:/Users/megan/Documents/GitHub/MixCat/data")
simdata <- readRDS("SimulatedData.RDS")



df <- simdata$simdf
df$group <- as.numeric(as.factor(df$age_group))
ageL <- c(0,10,20,30,40,50,60,70,80)
ageU <- c(9,19,29,39,49,59,69,79,90)

# data inputs
moddata <- list(N = nrow(df),
                y = df$titer,
                ageG = df$group,
                NAgeG = length(ageL),
                ageL = ageL,
                ageU = ageU,
                n_age = as.integer(tabulate(df$group, nbins = length(ageL))),
                y_fit = seq(min(df$titer)-1, max(df$titer)+1, 0.2),
                predL = length(seq(min(df$titer)-1, max(df$titer)+1, 0.2)),
                prior_means=c(2,3,1,1),
                prior_sds=c(1,1,1,1))


# fit
setwd("C:/Users/megan/Documents/GitHub/MixCat/StanModels")
fit <- stan(file='MixtureCatalytic.stan', data=moddata, chains=3, cores=1,
            iter=1000, warmup=500, refresh=10)
fitmix <- stan(file='Mixture.stan', data=moddata, chains=3, cores=1,
            iter=1000, warmup=500, refresh=10)

# check convergence
draws <- rstan::extract(fit)
draws <- rstan::extract(fitmix)

color_scheme_set("mix-blue-red")
mcmc_trace(fit, regex_pars=c("lam","mu","sd"))
mcmc_trace(fitmix, regex_pars=c("lam","mu","sd"))


# outputs
plot_seroprev(draws, moddata)
plot_dist_fit(draws, moddata)
plot_mean_titer(draws, moddata)
plot_titer_pred(draws, moddata)
plot_prob_seropos(draws, moddata)

extract_seroprev(draws, moddata)
sero_ests <- extract_seroprev_mix(draws, moddata)
sero_ests
ggplot(sero_ests, aes(age_label, seroprev))+ geom_point()+
  geom_linerange(aes(ymin=criL, ymax=criU))+ theme_minimal()+ ylim(0,NA)

plot_seroprev_mix(draws, moddata)
