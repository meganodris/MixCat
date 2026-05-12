# MixCat

MixCat includes code for fitting **mixture** and **mixture-catalytic** models to antibody titer data using Stan. MixCat estimates age-specific seroprevalence using two-component Gaussian mixture models, assuming the population antibody titer distribution is composed of seropositive and seronegative individuals.

## Models

| Stan model | Overview |
|:------------|:----------------------------------------------------------|
| `Mixture` | Age-specific seroprevalence is estimated directly from a mixture model (no endemicity assumption) |
| `MixtureCatalytic` | Age-specific seroprevalence and force of infection (FOI) are jointly inferred from a mixture-catalytic model (assumes endemic transmission) |

## When to use each model

In general, these models should only be used when the population antibody titer distribution has a bimodal pattern (i.e. some separation of antibody responses between positive and negative individuals) OR if the approximate range of antibody responses for positive or negative individuals is known from control samples, which can be used to inform model priors.

**`Mixture`** makes no assumption about pathogen transmission dynamics and can therefore be used in all scenarios. Age-specific seroprevalence estimates from this model can give an agnostic indication of past transmission dynamics (e.g. endemic or epidemic) for pathogens that induce life-long serological responses.

**`MixtureCatalytic`** can be used when endemic transmission is a reasonable assumption and/or when a trend of increasing seroprevalence by age was observed from the **`mixture`** model.

Both models support the estimation of age-specific seroprevalence / FOI by sub-populations (e.g. urban/rural locations, sex, socioeconomic status).

#### Assumptions

-   Both models assume the population to be composed of seronegative and seropositive individuals. Model estimates may be unreliable if true seroprevalence is close to 0% or 100%.

-   Both models assume no antibody waning or seroreversion. Seroprevalence will be underestimated for pathogens/antigens where antibodies wane over short timescales.

-   **`MixtureCatalytic`** assumes a constant endemic force of infection (FOI), i.e. no differences in the risk of infection by age or over time. For pathogens with seasonal or annual fluctuations in transmission intensity, the FOI parameter will represent a long-term average annual estimate of this transmission intensity.

## Setup

Clone the repository and source the helper functions:

``` r
# install.packages(c("rstan", "bayesplot", "tidyverse"))

source("utils.R")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
```

## Quick start

``` r
# Prepare Stan input
model_data <- make_model_data(
  titer     = df$titer,
  age_group = df$age_group_int,
  ageL      = c(0, 10, 20, 30, 40, 50, 60, 70, 80),
  ageU      = c(9, 19, 29, 39, 49, 59, 69, 79, 90)
)

# Fit the catalytic model
fit <- stan(
  file   = "StanModels/MixtureCatalytic.stan",
  data   = model_data,
  chains = 3, iter = 2000, warmup = 2000
)

# Plot seroprevalence curve with 95% CrI
draws <- rstan::extract(fit)
plot_seroprev(draws, model_data)
```

See `vignettes/introduction.Rmd` for a full worked example including convergence checks, model comparison, and all available plot and extract functions.
