library(rstan)


#----- Function to fit mixture catalytic model -----#
## intakes log transformed data
## should it specify priors?

fit_mixcat <- function(data){


  # age labels & groups
  age_groups <- paste(data$age_min, data$age_max, sep="-")
  NAgeG <- length(unique(data$age_min))
  ageG <- as.numeric(as.factor(data$age_min))
  ageL <- sort(unique(data$age_min))
  ageU <- sort(unique(data$age_max))

  # titer values for prediction
  y_fit <- seq(min(log_titer), max(log_titer), 0.5)
  predL <- length(y_fit)


  # data inputs
  data <- list(y=log_titer,
               N=length(log_titer),
               ageG=ageG,
               NAgeG=NAgeG,
               ageL=ageL,
               ageU=ageU,
               y_fit=y_fit,
               predL=predL)


  # setup stan model
  mod <- stan_model(paste0(file_path, "cat_mix.stan"))


  # fit model
  fit <- sampling(mod, data=data, chains=3, iter=5000, warmup=3000, refresh=10)


  # check convergence
  color_scheme_set("mix-blue-red")
  trace_plot <- mcmc_trace(fit, regex_pars=c("lambda","mu0","mu1","sd0","sd1"))

  # return results
  results <- list(mod_fit=fit,
                  trace_plot=trace_plot)

}


