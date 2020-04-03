# time-series estimates of psi


# probability that a case will have either died or recovered by a given day
# post-hospitalisation
probability <- function(day, mean_hdt = 13, median_hdt = 9.1) {
  
  # parameters of lognormal delay distribution
  mu_hdt <- log(median_hdt)
  sigma_hdt <- sqrt(2*(log(mean_hdt) - mu_hdt))
  
  # probability that the delay between hospitalisation and death is 'day' days
  plnorm(day + 1, mu_hdt, sigma_hdt) - plnorm(day, mu_hdt, sigma_hdt)
  
}

# compute the (non-cumulative) number of cases for which we would know the outcome
cases_known_outcome <- function(daily_cases){
  
  n_days <- length(daily_cases)
  days <- seq_len(n_days)
  
  # get a probability of delaying each of these number of days (starting from 0)
  delay_probs <- probability(days - 1)
  
  # effective number of cases we would have known about
  cases_known <- rep(0, n_days)
  
  # disaggregate these cases across subsequent days
  for(day in days){
    
    days_ahead <- seq_len(n_days - day + 1) - 1
    day_assign <- day + days_ahead
    
    # get the number of cases on each subsequent day (probabilities indexed from 1)
    new_partial_cases <- daily_cases[day] * delay_probs[days_ahead + 1]
    
    # assign them
    cases_known[day_assign] <- cases_known[day_assign] + new_partial_cases
    
  }
  
  cases_known
  
}

# download the latest Johns Hopkins data (it's in the package directly, so need to reinstall regularly)
remotes::install_github("RamiKrispin/coronavirus")
library(coronavirus)

# subset the data to Aus states
library(dplyr)
aus <- coronavirus %>%
  filter(Country.Region == "Australia") %>%
  transmute(state = Province.State, date, count = cases, type) %>%
  group_by(state)

# For each of the states, get the time series of cases, deaths, and expected
# number of cases with known outcomes (remove any negative cases or deaths).
aus_timeseries <- aus %>%
  tidyr::pivot_wider(names_from = type, values_from = count) %>%
  select(-recovered, cases = confirmed, deaths = death) %>%
  mutate(cases = pmax(0, cases),
         deaths = pmax(0, deaths)) %>%
  mutate(cases_known_outcome = cases_known_outcome(cases))

# get wide form versions of the deaths, and cases with known outcomes
death_table <- aus_timeseries %>%
  select(-cases, -cases_known_outcome) %>%
  tidyr::pivot_wider(names_from = state, values_from = deaths)

death_matrix <- death_table %>%
  select(-date) %>%
  as.matrix

cases_known_table <- aus_timeseries %>%
  select(-cases, -deaths) %>%
  tidyr::pivot_wider(names_from = state, values_from = cases_known_outcome)

cases_known_matrix <- cases_known_table %>%
  select(-date) %>%
  as.matrix

# check the dates match
stopifnot(identical(death_table$date, cases_known_table$date))

# check there are no deaths on days without cases that have known outcomes
stopif <- function(expr) {stopifnot(!expr)}
stopif(any(death_matrix > 0 & cases_known_matrix == 0))

# build model for contemporary observed number of deaths and expected number of
# deaths:
#   deaths_t ~ Poisson(expected_deaths_t)
#   expected_deaths_t = cases_known_outcomes_t * CFR_t
#   CFR_t = baseline_CFR / psi_t

library(greta)

n_states <- ncol(death_matrix)
n_times <- nrow(death_matrix)

# Distribution over plausible baseline CFR values from China study. The 95% CIs
# are symmetric around the estimate, so we assume it's a an approximately
# Gaussian distribution, truncated to allowable values. 
baseline_cfr_perc <- normal(1.38, 0.077, dim = n_states, truncation = c(0, 1))

# A separate reporting rate for each country, with all reporting rates a priori
# equally as likely.
sigma_psi <- normal(0, 1, truncation = c(0, Inf))
mu_psi <- normal(0, 1)
z_raw <- normal(0, 1, dim = n_states)
z <- mu_psi + z_raw * sigma_psi
psi <- iprobit(z)

# visualise prior:
# sim <- calculate(psi[1], nsim = 10000)
# hist(sim$`psi[1]`)

# expand out psi to matrix, to mock-up temporally-varying model
log_psi <- log(psi)
log_psi_t <- sweep(zeros(n_times, n_states), 2, log_psi, "+")

# compute CFR for each state and time
log_baseline_cfr <- log(baseline_cfr_perc) - log(100)
log_cfr_t <- sweep(-log_psi_t, 2, log_baseline_cfr, "+")

# define sampling distribution, subsetting to where cases_known_matrix > 0
# do the exponentiation down here, so greta can use log version in poisson density
some_cases_known <- which(cases_known_matrix > 0)

log_expected_deaths <- log_cfr_t[some_cases_known] +
  log(cases_known_matrix)[some_cases_known]
expected_deaths <- exp(log_expected_deaths)

observed_deaths <- death_matrix[some_cases_known]
distribution(observed_deaths) <- poisson(expected_deaths)

set.seed(2020-04-02)
m <- model(psi)
draws <- mcmc(m, chains = 10, n_samples = 3000)

# check convergence before continuing
coda::gelman.diag(draws)
sry <- summary(draws)

draws_mat <- as.matrix(draws)
colnames(draws_mat) <- colnames(death_matrix)
library(bayesplot)
bayesplot::color_scheme_set("blue")
bayesplot::mcmc_intervals(draws_mat) + ggplot2::theme_minimal()

reporting_rate <- data.frame(estimate = colMeans(draws_mat),
                             lower = apply(draws_mat, 2, quantile, 0.025),
                             upper = apply(draws_mat, 2, quantile, 0.975))


knitr::kable(round(reporting_rate, 3))
