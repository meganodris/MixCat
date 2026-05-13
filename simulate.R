# Simulate antibody titer data from a mixture-catalytic model with two
# location strata (Urban / Rural) having different forces of infection.
#
# The population antibody titer distribution is a two-component Gaussian
# mixture: seronegative individuals ~ N(mu0, sd0) and seropositive
# individuals ~ N(mu0 + mu1, sd1), where mu1 is the increment above
# the seronegative mean.
#
# Serostatus is assigned via a catalytic model: an individual of age a
# in stratum g is seropositive with probability 1 - exp(-lambda_g * a).


# --- Parameters --------------------------------------------------------------

N <- 400  # total sample size (split roughly equally between strata)

# Titer distribution parameters (shared across strata)
mu0 <- 2    # seronegative mean
mu1 <- 2.5  # seropositive increment (seropositive mean = mu0 + mu1)
sd0 <- 0.7  # seronegative SD
sd1 <- 0.9  # seropositive SD

# Force of infection per stratum (annual probability of infection per susceptible population)
lam_urban <- 0.06
lam_rural <- 0.02


# --- Sample population -------------------------------------------------------

set.seed(42)

age      <- floor(runif(N, 0, 90))
age_group <- cut(age, breaks = c(-Inf, 9, 19, 29, 39, 49, 59, 69, 79, Inf),
                 labels = paste(c(0, 10, 20, 30, 40, 50, 60, 70, 80),
                                c(9, 19, 29, 39, 49, 59, 69, 79, 90), sep = "-"))
location <- sample(c("Urban", "Rural"), N, replace = TRUE)
lam      <- ifelse(location == "Urban", lam_urban, lam_rural)


# --- Simulate serostatus and antibody titers ------------------------------------------

# Seropositive probability follows catalytic model: P(pos | age, stratum) = 1 - exp(-lambda * age)
seropos <- rbinom(N, 1, prob = 1 - exp(-lam * age))

titer <- ifelse(seropos == 0,
                rnorm(N, mu0, sd0),
                rnorm(N, mu0 + mu1, sd1))

simdf <- data.frame(age       = age,
                    age_group = age_group,
                    location  = location,
                    status    = seropos,
                    titer     = titer)


# --- Save --------------------------------------------------------------------

output <- list(
  simdf    = simdf[, c("age_group", "location", "titer")],
  truepars = data.frame(
    lam_urban = lam_urban,
    lam_rural = lam_rural,
    mu0       = mu0,
    mu1       = mu1,
    sd0       = sd0,
    sd1       = sd1
  )
)

saveRDS(output, "data/SimulatedData.RDS")
