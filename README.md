
# MixCat

MixCat fits **mixture** and **mixture-catalytic** models to antibody
titer data using Stan. Both models represent the population titer
distribution as a two-component Gaussian mixture, where the two
components correspond to seronegative (μ₀, σ₀) and seropositive (μ₀ +
μ₁, σ₁) individuals. Age-specific seroprevalence is estimated either
directly from the mixture components (`Mixture`) or jointly with a force
of infection (FOI) parameter via a catalytic model (`MixtureCatalytic`).
See the introductory
[vignette](https://raw.githack.com/meganodris/MixCat/main/vignettes/Vignette.html)
for a full worked example including the **`Mixture`** and
**`Mixture Catalytic`** models, convergence diagnostics, and all
available extract and plot functions.

## Models

| Stan model | Overview |
|:---|:---|
| `Mixture` | Age-specific seroprevalence estimated directly from the mixture (no endemicity assumption) |
| `MixtureCatalytic` | Age-specific seroprevalence and FOI jointly inferred from a mixture-catalytic model (assumes endemic transmission) |

Both models support stratified estimates by sub-population
(e.g. urban/rural, sex) via an optional `group` argument.

## When to use each model

These models should only be used when the population antibody titer
distribution shows a bimodal pattern (i.e. some separation between
seropositive and seronegative individuals), or when the approximate
titer range for each group is known from control samples and can be used
to inform priors.

**`Mixture`** makes no assumption about transmission dynamics and can be
used in all scenarios. Age-specific seroprevalence estimates from this
model can give an agnostic indication of past transmission dynamics
(e.g. endemic or epidemic) for pathogens that induce life-long
serological responses.

**`MixtureCatalytic`** can be used when endemic transmission is a
reasonable assumption, for example when a trend of increasing
seroprevalence by age was observed from the **`Mixture`** model.

#### Assumptions

- Both models assume the population is composed of seronegative and
  seropositive individuals. Estimates may be unreliable when true
  seroprevalence is close to 0% or 100%.
- Both models assume no antibody waning or seroreversion. Seroprevalence
  will be underestimated for pathogens where antibodies wane over short
  timescales.
- **`MixtureCatalytic`** assumes a constant endemic FOI with no
  differences in infection risk by age or over time. For pathogens with
  seasonal or annual fluctuations, the FOI represents a long-term
  average.

## Setup

Clone the repository and source the helper functions:

``` r
# install.packages(c("rstan", "bayesplot", "tidyverse"))
library(rstan)
library(bayesplot)
library(tidyverse)

source("utils.R")
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
```

## Data preparation

MixCat requires three inputs:

- **`titer`**: numeric vector of antibody titer values (log-transform
  first for measurements on a large scale)
- **`age_group`**: integer vector assigning each individual to an age
  group (youngest group = 1)
- **`group`** *(optional)*: integer vector identifying sub-population
  groups for stratified estimates

`make_model_data()` assembles the Stan input list. The `prior_means` and
`prior_sds` arguments specify Normal priors for μ₀, μ₁, σ₀, σ₁ — choose
values consistent with the scale of your titer data. Note that μ₁ is the
*increment* from seronegative to seropositive, not the absolute
seropositive mean.

``` r
df <- load_example_data()
head(df)
```

    ##   age_group location     titer age_group_int location_int
    ## 1     80-90    Rural 1.6964506             9            1
    ## 2     80-90    Rural 3.0761887             9            1
    ## 3     20-29    Rural 0.4808274             3            1
    ## 4     70-79    Urban 4.1101071             8            2
    ## 5     50-59    Rural 6.5925526             6            1
    ## 6     40-49    Rural 4.9717100             5            1

``` r
ageL <- c(0, 10, 20, 30, 40, 50, 60, 70, 80)
ageU <- c(9, 19, 29, 39, 49, 59, 69, 79, 90)

model_data <- make_model_data(
  titer       = df$titer,
  age_group   = df$age_group_int,
  ageL        = ageL,
  ageU        = ageU,
  group       = df$location_int,
  prior_means = c(2, 3, 1, 1),
  prior_sds   = c(1, 1, 1, 1)
)
```

## Model fitting

``` r
fit <- stan(
  file    = "StanModels/MixtureCatalytic.stan",
  data    = model_data,
  chains  = 3, cores  = 1,
  iter    = 1500, warmup = 500
)

draws <- rstan::extract(fit)
```

## Results

FOI and seroprevalence estimates can be extracted with the `extract_foi`
and `extract_sero` functions:

``` r
extract_foi(draws, model_data, group_labels = c("Rural", "Urban"))
extract_sero(draws, model_data, group_labels = c("Rural", "Urban"))
```

The functions `extract_sero_age_mix` and `extract_sero_age_cat` return
age-specific seroprevalence estimates for the **`Mixture`** and
**`MixtureCatalytic`** models respectively. These estimates can also be
visualised using the corresponding plotting functions,
e.g. `plot_sero_age_cat`.

``` r
extract_sero_age_cat(draws, model_data, group_labels = c("Rural", "Urban"))
plot_sero_age_cat(draws, model_data, group_labels = c("Rural", "Urban"))
```

![](vignettes/Vignette_files/figure-html/seroprev-cat-1.png)<!-- -->

Model fit to the observed titer distributions can be visualised using
the `plot_dist_fit` and `plot_mean_titer` functions.

``` r
plot_dist_fit(draws, model_data)
```

![](vignettes/Vignette_files/figure-html/dist-fit-cat-1.png)<!-- -->

``` r
plot_mean_titer(draws, model_data, group_labels = c("Rural", "Urban"))
```

![](vignettes/Vignette_files/figure-html/dist-fit-cat-2.png)<!-- -->

See the introductory
[vignette](https://raw.githack.com/meganodris/MixCat/main/vignettes/Vignette.html)
for a full worked example including the **`Mixture`** and
**`Mixture Catalytic`** models, convergence diagnostics, and all
available extract and plot functions.
