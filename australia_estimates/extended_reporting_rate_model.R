# Time-series estimates of reporting rates for Australian states

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

# a timeseries of reporting rates for each state, modelled as hierarchical Gaussian processes
library(greta.gp)

# squared-exponential GP kernels with unknown parameters for the national and
# state-level processes
national_lengthscale_raw <- normal(0, 1, truncation = c(0, Inf))
national_lengthscale <- national_lengthscale_raw * 100
national_sigma <- normal(0, 1, truncation = c(0, Inf))
national_kernel <- rbf(lengthscales = national_lengthscale,
                       variance = national_sigma ^ 2)

state_lengthscale_raw <- normal(0, 1, truncation = c(0, Inf))
state_lengthscale <- state_lengthscale_raw * 100
state_sigma <- normal(0, 1, truncation = c(0, Inf))
state_kernel <- rbf(lengthscales = state_lengthscale,
                       variance = state_sigma ^ 2)

# a set of inducing points at which to estimate the GPs (subset of regressors
# approximation)
# put an inducing point on the last time point (most recent date), but otherwise
# space them out
times <- seq_len(n_times)
n_inducing <- 10
inducing_points <- seq(min(times), max(times), length.out = n_inducing + 1)[-1]

tol <- 1e-6
# GP for the national mean effect
mu <- greta.gp::gp(times, national_kernel, tol = tol)

# GPs for the state deviations (manually defined as multiple GPs at once isn't
# yet possible in greta.gp)
v <- normal(0, 1, dim = c(n_inducing, n_states))
Kmm <- state_kernel(inducing_points)
Kmm <- Kmm + diag(n_inducing) * tol
Lm <- t(chol(Kmm))
Kmn <- state_kernel(inducing_points, times)
A <- forwardsolve(Lm, Kmn)
z_state <- t(A) %*% v

# add the mean effect on
z <- sweep(z_state, 1, mu, "+")

# convert to probabilities
psi <- iprobit(z)

# # visualise prior:
# nsim <- 100
# sims <- calculate(psi, nsim = nsim)[[1]]
# plot(sims[1, , 1] ~  times, type = "n", ylim = c(0, 1))
# for(i in seq_len(nsim)) lines(sims[i, , 1] ~ times, lwd = 0.2)

# Distribution over plausible baseline CFR values from China study. The 95% CIs
# are symmetric around the estimate, so we assume it's a an approximately
# Gaussian distribution, truncated to allowable values. 
baseline_cfr_perc <- normal(1.38, 0.077, dim = n_states, truncation = c(0, 1))

# compute CFR for each state and time
log_baseline_cfr <- log(baseline_cfr_perc) - log(100)
log_psi <- log(psi)
log_cfr <- sweep(-log_psi, 2, log_baseline_cfr, "+")

# define sampling distribution, subsetting to where cases_known_matrix > 0
# do the exponentiation down here, so greta can use log version in poisson density
some_cases_known <- which(cases_known_matrix > 0)

log_expected_deaths <- log_cfr[some_cases_known] +
  log(cases_known_matrix)[some_cases_known]
expected_deaths <- exp(log_expected_deaths)

observed_deaths <- death_matrix[some_cases_known]
distribution(observed_deaths) <- poisson(expected_deaths)

set.seed(2020-04-02)
m <- model(national_lengthscale,
           national_sigma,
           state_lengthscale,
           state_sigma)
draws <- mcmc(m, chains = 10, n_samples = 3000, one_by_one = TRUE)

# check convergence before continuing
coda::gelman.diag(draws)
coda::effectiveSize(draws)

png("reporting_rate_timeseries_by_state.png",
    width = 1200,
    height = 1800,
    pointsize = 30)
par(mfrow = c(4, 2))

state_names <- colnames(death_matrix)
for (i in seq_len(n_states)) {
  
  # subset to the time after the state's first case
  start <- which(cases_known_matrix[, i] > 0)[1]
  index <- start:n_times
  
  # predict each state's timeseries
  draws <- calculate(psi[index, i], values = draws)
  draws_mat <- as.matrix(draws)
  mean <- colMeans(draws_mat)
  ci <- apply(draws_mat, 2, quantile, c(0.025, 0.975))
  iqr <- apply(draws_mat, 2, quantile, c(0.25, 0.75))
  
  times_plot <- times[index]
  # subset times to when this state first saw a case
  
  plot(mean ~ times_plot, type = "n",
       ylim = c(0, 1),
       xlim = range(times),
       axes = FALSE, xlab = "date", ylab = "reporting rate")
  polygon(x = c(times_plot, rev(times_plot)),
          y = c(ci[1, ], rev(ci[2, ])),
          col = blues9[2], lty = 0)
  polygon(x = c(times_plot, rev(times_plot)),
          y = c(iqr[1, ], rev(iqr[2, ])),
          col = blues9[3], lty = 0)
  lines(mean ~ times_plot, lwd = 4, col = blues9[6])
  axis(2, las = 2)
  axis(1, at = inducing_points, labels = min(aus_timeseries$date) + inducing_points - 1)
  
  title(main = state_names[i])
  abline(v = times_plot[1], lwd = 1.5, col = "red")
  text(x = times_plot[1], y = 1.08,
       labels = "first case reported",
       xpd = NA, col = "red", cex = 0.8)
  
}
dev.off()

# calculate latest estimates
draws <- calculate(psi[72, ], values = draws)
draws_mat <- as.matrix(draws)
colnames(draws_mat) <- state_names
library(bayesplot)
library(ggplot2)
bayesplot::color_scheme_set("blue")
bayesplot::mcmc_intervals(draws_mat, point_est = "mean", prob = 0.5, prob_outer = 0.95) +
  ggplot2::xlim(0, 1) +
  ggplot2::theme_minimal() +
  ggplot2::ggtitle(paste("estimated reporting rates on",
                         max(aus_timeseries$date)))
ggplot2::ggsave("latest_reporting_rates.png", scale = 0.8)

reporting_rate <- data.frame(mean = colMeans(draws_mat),
                             lower_50 = apply(draws_mat, 2, quantile, 0.25),
                             upper_50 = apply(draws_mat, 2, quantile, 0.75),
                             lower_95 = apply(draws_mat, 2, quantile, 0.025),
                             upper_95 = apply(draws_mat, 2, quantile, 0.975))

knitr::kable(round(reporting_rate, 3))
