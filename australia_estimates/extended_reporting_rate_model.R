# Ideally would like the number of cases and deaths stratified by whether
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

# compute the cumulative underestimation correction for CFR
cumulative_underestimation <- function(daily_cases){
  
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


# # check this against the cumulative version
# daily_cases <- rpois(20, exp(1:20/10))
# for(i in seq_along(daily_cases)) {
#   daily_cases_sub <- daily_cases[1:i]
#   new <- sum(cases_known_outcome(daily_cases_sub)) / sum(daily_cases_sub)
#   old <- cumulative_underestimation(daily_cases_sub)
#   print(max(abs(new - old)) < 1e-9)
# }

# download the latest Johns Hopkins data (it's in the package directly, so need to reinstall regularly)
remotes::install_github("RamiKrispin/coronavirus")
library(coronavirus)

# subset the data to Aus states
library(dplyr)
aus <- coronavirus %>%
  filter(Country.Region == "Australia") %>%
  transmute(state = Province.State, date, count = cases, type) %>%
  group_by(state)



# how to handle timesteps with zero cases?

# This binomial sampling is pretty weird. Instead, we'd want to estimate for a
# pool of cases whether or not they subsequently died. But we can't do that with
# integers because we don't have the outcomees for individuals.

# Instead, could do the Poisson approximation to the binomial, so that the
# number of cases with known outcome is the rate.

#   deaths_t ~ Poisson(cases_known_outcomes * CFR)
#   cCFR_t = baseline_CFR / psi_t



# get (non-cumulative) known cases at each time t, divide these by the
# (non-cumulative) cases counts at that time to underestimation at each time t.
# Them include with contemporary case counts in the following model:

#   deaths_t ~ Binomial(cases_t, nCFR_t)
#   nCFR_t = cCFR_t * underestimation_t
#   cCFR_t = baseline_CFR / psi_t
#   underestimation_t = cases_known_t / cases_t

# rewrite the cases_known estimation method to be non-cumulative, then compute
# the cumulative sum, and test against the existing implementation.

# write out the greta timeseries model using a sort of hierarchical GP for
# probit(psi) using greta.gp? same kernel (so only one inversion per iteration),
# multiple draws (one for each state, plus a national average), then add the
# national one to each of the others.

# or just have a national-level timeseries?
# can we include a regression component?

# can be back-calculate a national estimate?
# weighted sum of psis, using known cases as weights?


# given the 


# for each state get the total number of cases and deaths
totals <- aus %>%
  summarise(total_deaths = sum(count[type == "death"]),
            total_cases = sum(count[type == "confirmed"]))

# for each state, get the underestimation parameter and append to totals, for use in modelling
data <- aus %>%
  left_join(totals) %>%
  filter(type == "confirmed") %>%
  summarise(underestimation = cumulative_underestimation(count)) %>%
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

# visualise prior:
# sim <- calculate(psi[1], nsim = 10000)
# hist(sim$`psi[1]`)

# Observation model:
#   total_deaths ~ Poisson(expected_deaths)
#   expected_deaths = cases_known_outcomes * CFR
#   CFR = baseline_CFR / psi
data$cases_known_outcomes <- data$underestimation * data$total_cases
log_expected_deaths <- log(data$cases_known_outcomes) + log(baseline_cfr) - log(100) - log(psi)
expected_deaths <- exp(log_expected_deaths)
distribution(data$total_deaths) <- poisson(expected_deaths)



# # Observation model. Equivalent to, but more numerically stable than:
# #   nCFR <- (baseline_cfr / 100) * underestimation / psi
# log_nCFR <- log(baseline_cfr) - log(100) + log(data$underestimation) - log(psi)
# nCFR <- exp(log_nCFR)
# distribution(data$total_deaths) <- binomial(data$total_cases, nCFR)

set.seed(2020-04-02)
m <- model(psi)
draws <- mcmc(m, chains = 10, n_samples = 3000)

# check convergence before continuing
coda::gelman.diag(draws)
sry <- summary(draws)

draws_mat <- as.matrix(draws)
colnames(draws_mat) <- data$state
library(bayesplot)
bayesplot::color_scheme_set("blue")
bayesplot::mcmc_intervals(draws_mat) + ggplot2::theme_minimal()

reporting_rate <- data.frame(estimate = colMeans(draws_mat),
                             lower = apply(draws_mat, 2, quantile, 0.025),
                             upper = apply(draws_mat, 2, quantile, 0.975))


knitr::kable(round(reporting_rate, 3))
