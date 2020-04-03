# Ideally would like need the number of cases and deaths stratified by whether
# when detected they were:
#       - in unsupervised quarantine
#       - in supervised quarantine
#       - not in quarantine
# But more plausibly wee could estimate this by whether they were imports, known
# contacts, or unknown transmission)

# Estimate variation in the reporting rate between states, using a hierarchical
# model to share information. Might expect that NSW and Vic are getting overwhelmed.
# Use number of tests as a covariate?

# get cumulative underestimation parameter for each state from these case data
# get total number of imports/local/community cases per state from current national sitrep.
# get split of the 24 deaths between imports/local/community

# probability that a case will have either died or recovered by a given day
# post-hospotalisation
probability <- function(day, mean_hdt = 13, median_hdt = 9.1) {
  
  # parameters of lognormal delay distribution
  mu_hdt <- log(median_hdt)
  sigma_hdt <- sqrt(2*(log(mean_hdt) - mu_hdt))
  
  # probability that the delay between hospitalisation and death is 'day' days
  plnorm(day + 1, mu_hdt, sigma_hdt) - plnorm(day, mu_hdt, sigma_hdt)
  
}

# compute the underestimation correction for CFR
underestimation <- function(daily_cases, total_deaths){
  
  total_cases <- sum(daily_cases)
  
  # get cases we woud have known about
  cumulative_known <- 0
  
  # Sum over cases up to time tt
  for(day in seq_along(daily_cases)){
    
    known_day <- 0 # number of cases with known outcome at time ii
    
    # distribute cases across dayas by which we would know if they died or survived
    for(delay in 0:(day - 1)){
      known_delay <- (daily_cases[day - delay]*probability(delay))
      known_day <- known_day + known_delay
    }
    
    cumulative_known <- cumulative_known + known_day # Tally cumulative known
  }
  
  # underestimation of CFR due to unknown outcomes
  cumulative_known / total_cases
  
}

# remotes::install_github("RamiKrispin/coronavirus")

# subset the data to Aus states
library(coronavirus)
library(dplyr)
aus <- coronavirus %>%
  filter(Country.Region == "Australia") %>%
  transmute(state = Province.State, date, count = cases, type) %>%
  group_by(state)

# for each state get the total number of cases and deaths
totals <- aus %>%
  summarise(total_deaths = sum(count[type == "death"]),
            total_cases = sum(count[type == "confirmed"]))

# for each state, get the underestimation parameter and append to totals, for use in modelling
data <- aus %>%
  left_join(totals) %>%
  filter(type == "confirmed") %>%
  summarise(underestimation = underestimation(count, total_deaths)) %>%
  left_join(totals)

library(greta)

n <- nrow(data)

# Distribution over plausible baseline CFR values from China study. The 95% CIs
# are symmetric around the estimate, so we assume it's a an approximately
# Gaussian distribution, truncated to allowable values. 
baseline_cfr <- normal(1.38, 0.077, dim = n, truncation = c(0, 1))

# A separate reporting rate for each country, with all reporting rates a priori
# equally as likely.
sigma_psi <- normal(0, 1, truncation = c(0, Inf))
mu_psi <- normal(0, 1)
z_raw <- normal(0, 1, dim = n)
z <- mu_psi + z_raw * sigma_psi
psi <- iprobit(z)
sim <- calculate(psi[1], nsim = 10000)

# visualise prior:
# hist(sim$`psi[1]`)

# Observation model. Equivalent to, but more numerically stable than:
#   nCFR <- (baseline_cfr / 100) * underestimation / psi
log_nCFR <- log(baseline_cfr) - log(100) + log(data$underestimation) - log(psi)
nCFR <- exp(log_nCFR)
distribution(data$total_deaths) <- binomial(data$total_cases, nCFR)

set.seed(2020-04-02)
m <- model(psi)
draws <- mcmc(m, chains = 10, n_samples = 3000)

# check convergence before continuing
coda::gelman.diag(draws)
sry <- summary(draws)

str(draws)
draws_mat <- as.matrix(draws)
colnames(draws_mat) <- data$state
library(bayesplot)
bayesplot::color_scheme_set("blue")
bayesplot::mcmc_intervals(draws_mat) + ggplot2::theme_minimal()

reporting_rate <- data.frame(estimate = colMeans(draws_mat),
           lower = apply(draws_mat, 2, quantile, 0.025),
           upper = apply(draws_mat, 2, quantile, 0.975))


knitr::kable(round(reporting_rate, 3))
  
  data %>%
  mutate(cCFR = (total_deaths / total_cases) / underestimation)

# summarise estimates
allTogetherClean2$underreporting_estimate <- sry$statistics[, "Mean"]
allTogetherClean2$lower <- sry$quantiles[, "2.5%"]
allTogetherClean2$upper <- sry$quantiles[, "97.5%"]



