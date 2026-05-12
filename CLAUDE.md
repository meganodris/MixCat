# MixCat

R package for fitting **mixture catalytic models** to antibody titer data using Stan (via rstan).

## What it does

Jointly estimates:
- **Force of infection (FOI)** via a catalytic model → smooth age-specific seroprevalence
- **Titer distributions** for seronegative and seropositive individuals via a two-component Gaussian mixture

## Stan models (`StanModels/`)

| File | Purpose |
|---|---|
| `MixtureCatalytic.stan` | Single population, endemic transmission, single FOI |
| `Mixture.stan` | Age-specific seroprevalence estimated directly (no catalytic assumption) |
| `MixtureCatalyticGrouped.stan` | Multiple populations with separate FOI per group (e.g. urban/rural) |

Pre-compiled `.rds` files exist for `MixtureCatalytic` and `Mixture`.

## Key files

- `utils.R` — all user-facing functions: `make_model_data()`, `make_model_data_grouped()`, and paired `extract_*` / `plot_*` functions for seroprevalence, distribution fit, mean titer, individual predictions, and P(seropositive)
- `catmix_functions.R` — older wrapper `fit_mixcat()`; likely to be removed or replaced before publication
- `simulate.R` — script that generated `data/SimulatedData.RDS` (N=500, λ=0.07, μ₀=2, μ₁=2.5, σ₀=0.7, σ₁=0.9)
- `FitMixCat_working.R` — working/scratch fitting script
- `vignettes/introduction.Rmd` — primary user-facing vignette demonstrating the full workflow

## Model parameterisation notes

- Seropositive mean is parameterised as `mu0 + mu1` (mu1 is the *increment*, not the absolute mean)
- `pC` is a 2×N array of log-component probabilities; `log_lik` uses `log_sum_exp` for numerical stability
- FOI prior: `lambda ~ exponential(4)` (mean ≈ 0.25)
- Mixture priors passed in as `prior_means` / `prior_sds` vectors (length 4: μ₀, μ₁, σ₀, σ₁)
- `sero[a]` in the Stan model is the average seroprevalence *within* an age group (interval-average of the catalytic curve), not evaluated at a single age

## Publication readiness checklist

- [ ] Decide fate of `catmix_functions.R` (retire or integrate)
- [ ] R package scaffolding (DESCRIPTION, NAMESPACE, roxygen docs)
- [ ] Unit/integration tests
- [ ] `MixtureCatalyticGrouped` vignette / worked example
- [ ] README
- [ ] Harden `simulate.R` (remove hardcoded `setwd`)
