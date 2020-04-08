# Time-series estimates of reporting rates for Australian states

# given a 'recipient' matrix with some missing values, fill in the missing
# values with the corresponding elements in the 'donor' matrix of the same
# dimensions, multiplied by the scalar 'correction' factor (e.g. a probabilisitc
# reassignment)
impute_values <- function(recipient, donor, correction) {
  missing <- which(is.na(recipient), arr.ind = TRUE)
  if(inherits(correction, "greta_array")) {
    recipient[missing] <- 0  
    recipient <- as_data(recipient)
  }
  recipient[missing] <- donor[missing] * correction
  recipient
}

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
  
  # temporarily omit NAs when filling these in
  missing <- is.na(daily_cases)
  daily_cases[missing] <- 0 
  
  # disaggregate these cases across subsequent days
  for(day in days){
    
    days_ahead <- seq_len(n_days - day + 1) - 1
    day_assign <- day + days_ahead
    
    # get the number of cases on each subsequent day (probabilities indexed from 1)
    new_partial_cases <- daily_cases[day] * delay_probs[days_ahead + 1]
    
    # assign them
    cases_known[day_assign] <- cases_known[day_assign] + new_partial_cases
    
  }
  
  cases_known[missing] <- NA
  cases_known
  
}

# build a circulant matrix of prob, with masked lower values, then
# matrix-multiply a date-by-state matrix of daily cases to get the cases with
# known cases on that date
cases_known_outcome_matrix <- function (daily_cases) {
  
  n_days <- NROW(daily_cases)
  days <- seq_len(n_days)
  
  # get a probability of delaying each of these number of days (starting from 0)
  delay_probs <- probability(days - 1)
  
  # build an upper-triangular circulant matrix (contribution of each day's cases
  # to each other's)
  mat <- matrix(0, n_days, n_days)
  indices <- col(mat) - row(mat) + 1
  mask <- indices > 0
  mat[mask] <- delay_probs[indices[mask]]
  
  # matrix-multiply to disaggregate them
  t(mat) %*% daily_cases
  
}

# fit a Subset-of-Regressors-approximated hierarchical GP with rbf + bias
# kernels for both the mean and states
hierarchical_gp <- function (times,
                             n_states,
                             lengthscale_mu,
                             lengthscale_z,
                             sigma_mu,
                             sigma_z,
                             sigma_mu_intercept = 0.7,
                             sigma_z_intercept = 0.7,
                             n_inducing = 5,
                             tol = 1e-4) {
  
  # GP kernels
  kernel_mu <-
    rbf(lengthscales = lengthscale_mu, variance = sigma_mu ^ 2) +
    bias(sigma_mu_intercept ^ 2)
  
  kernel_z <-
    rbf(lengthscales = lengthscale_z, variance = sigma_z ^ 2) + 
    bias(sigma_z_intercept ^ 2)
  
  # inducing points  
  inducing_points <- seq(min(times), max(times), length.out = n_inducing + 1)[-1]
  
  # mean GP
  mu <- greta.gp::gp(times, inducing = inducing_points, kernel_mu, tol = tol)
  
  # GPs for the states deviations (manually defined as multiple GPs at once isn't
  # yet possible in greta.gp)
  v <- normal(0, 1, dim = c(n_inducing, n_states))
  Kmm <- kernel_z(inducing_points)
  Kmm <- Kmm + diag(n_inducing) * tol
  Lm <- t(chol(Kmm))
  Kmn <- kernel_z(inducing_points, times)
  A <- forwardsolve(Lm, Kmn)
  z_deviations <- t(A) %*% v
  
  # add the mean effect onto the deviations and return
  sweep(z_deviations, 1, mu, "+")
  
}

# fit a Subset-of-Regressors-approximated temporal GP with rbf + bias kernel
temporal_gp <- function (times,
                         lengthscale,
                         sigma,
                         sigma_intercept = 0.7,
                         n_inducing = 5,
                         tol = 1e-4) {
  
  kernel <-
    rbf(lengthscales = lengthscale, variance = sigma ^ 2) +
    bias(sigma_intercept ^ 2)
  
  inducing_points <- seq(min(times), max(times), length.out = n_inducing + 1)[-1]
  
  greta.gp::gp(times, inducing = inducing_points, kernel, tol = tol)
  
}

# default hyperpriors for hierarchical GPs
default_hypers <- function() {
  list(lengthscale_mu = lognormal(4, 0.5),
       lengthscale_z = lognormal(4, 0.5),
       sigma_mu = lognormal(-2, 0.5),
       sigma_z = lognormal(-2, 0.5)
  )
}

# get date-by-state matrix of 'var' from data
date_state_matrix <- function(var, data) {
  matrix <- data %>%
    select(state, date, !!var) %>%
    tidyr::pivot_wider(names_from = state, values_from = !!var) %>%
    select(-date) %>%
    as.matrix
  rownames(matrix) <- as.character(unique(data$date))
  matrix
}

# download the latest Johns Hopkins data (it's in the package directly, so need to reinstall regularly)
remotes::install_github("RamiKrispin/coronavirus", upgrade = "never")
library(coronavirus)

# load data from covid19data.org.au (scraped by Chris Baker) with cases by
# source for some states.
source("load_data.R")
source_data <- load_data()

# subset the case data to Aus states
library(dplyr)
aus <- coronavirus %>%
  filter(Country.Region == "Australia") %>%
  transmute(state = Province.State, date, count = cases, type) %>%
  group_by(state)

# expand the type of count for each state, attach the sources where known, and
# compute cases with known outcomes
aus_timeseries <- aus %>%
  tidyr::pivot_wider(names_from = type, values_from = count) %>%
  select(-recovered, cases = confirmed, deaths = death) %>%
  left_join(source_data) %>%
  mutate_at(c("cases", "deaths", "overseas", "known_local", "unknown_local", "other"),
            ~ pmax(., 0))
  

# get date-by-state matrices for daily case counts by state for each source

# deaths and all cases (as in Johns Hopkins data)
deaths <- date_state_matrix("deaths", aus_timeseries)
cases <- date_state_matrix("cases", aus_timeseries)

# cases from each source, with some NAs
unknown_local <- date_state_matrix("unknown_local", aus_timeseries)
known_local <- date_state_matrix("known_local", aus_timeseries)
overseas <- date_state_matrix("overseas", aus_timeseries)
other <- date_state_matrix("other", aus_timeseries)

# impute missing values

# where no information is available (Queensland, or outside dates when sources
# were reported), disaggregate into the fraction that are local with unknown
# source and those that are 'other'; then disaggregate the 'other' into those
# that are overseas-acquired, vs. local with known source. For WA, 'other' is
# known, so just the second step.

library(greta.gp)

# estimating those parameters directly as with this model, the numbers are so
# large there's essentially no uncertainty on them. Could relax the model a bit
# later (e.g. hierarchical model by state on the proportions), and infer them

# p_other_is_overseas <- uniform(0, 1)
# p_all_is_unknown_local <- uniform(0, 1)

overseas_vs_other <- aus_timeseries %>%
  ungroup %>%
  select(overseas, known_local) %>%
  na.omit() %>%
  mutate(other = overseas + known_local) %>%
  summarise(successes = sum(overseas),
            trials = sum(other))

# distribution(overseas_vs_other$successes) <- binomial(overseas_vs_other$trials, p_other_is_overseas)

p_other_is_overseas <- overseas_vs_other$successes / overseas_vs_other$trials

unknown_local_vs_all <- aus_timeseries %>%
  ungroup %>%
  select(unknown_local, cases) %>%
  na.omit() %>%
  summarise(successes = sum(unknown_local),
            trials = sum(cases))

# distribution(unknown_local_vs_all$successes) <- binomial(unknown_local_vs_all$trials, p_all_is_unknown_local)
p_all_is_unknown_local <- unknown_local_vs_all$successes / unknown_local_vs_all$trials

# fill in values
unknown_local <- impute_values(unknown_local, cases, p_all_is_unknown_local)
other <- impute_values(other, cases, 1 - p_all_is_unknown_local)
overseas <- impute_values(overseas, other, p_other_is_overseas)
known_local <- impute_values(known_local, other, 1 - p_other_is_overseas)

# now disaggregate these cases using the delay distribution
unknown_local <- cases_known_outcome_matrix(unknown_local)
known_local <- cases_known_outcome_matrix(known_local)
overseas <- cases_known_outcome_matrix(overseas)

n_states <- ncol(deaths)
n_times <- nrow(deaths)
times <- seq_len(n_times)

# define hierarchical probit-GPs for reporting rates of unknown local, and known
# local cases. Assume reporting rate for overseas-acquired cases is perfect.

# unknown_local_hypers <- default_hypers()
# unknown_local_reporting_z <- hierarchical_gp(
#   times = times,
#   n_states = n_states,
#   lengthscale_mu = unknown_local_hypers$lengthscale_mu,
#   lengthscale_z = unknown_local_hypers$lengthscale_z,
#   sigma_mu = unknown_local_hypers$sigma_mu,
#   sigma_z = unknown_local_hypers$sigma_z
# )
unknown_local_gp_lengthscale <- lognormal(4, 0.5)
unknown_local_gp_sigma <- lognormal(-2, 0.5)
unknown_local_sigma <- normal(0, 0.5, truncation = c(0, Inf))
unknown_local_intercepts_raw <- normal(0, 1, dim = c(1, n_states))
unknown_local_intercepts <- unknown_local_intercepts_raw * unknown_local_sigma
unknown_local_gp <- temporal_gp(
  times = times,
  lengthscale = unknown_local_gp_lengthscale,
  sigma = unknown_local_gp_sigma
)
unknown_local_z <- kronecker(
  unknown_local_gp,
  unknown_local_intercepts,
  FUN = "+"
)
unknown_local_reporting <- iprobit(unknown_local_z)

# known_local_hypers <- default_hypers()
# known_local_reporting_z <- hierarchical_gp(
#   times = times,
#   n_states = n_states,
#   lengthscale_mu = known_local_hypers$lengthscale_mu,
#   lengthscale_z = known_local_hypers$lengthscale_z,
#   sigma_mu = known_local_hypers$sigma_mu,
#   sigma_z = known_local_hypers$sigma_z
# )
known_local_gp_lengthscale <- lognormal(4, 0.5)
known_local_gp_sigma <- lognormal(-2, 0.5)
known_local_sigma <- normal(0, 0.5, truncation = c(0, Inf))
known_local_intercepts_raw <- normal(0, 1, dim = n_states)
known_local_intercepts <- known_local_intercepts_raw * known_local_sigma
known_local_gp <- temporal_gp(
  times = times,
  lengthscale = known_local_gp_lengthscale,
  sigma = known_local_gp_sigma
)
known_local_z <- kronecker(
  known_local_gp,
  t(known_local_intercepts),
  FUN = "+"
)
known_local_reporting <- iprobit(known_local_z)

overseas_reporting <- ones(n_times, n_states)

# # visualise priors:
# nsim <- 10000
# sims <- calculate(unknown_local_reporting[n_times, 1],
#                   nsim = nsim)[[1]]
# plot(sims[1, , 1] ~  times, type = "n", ylim = c(0, 1))
# for(i in seq_len(nsim)) lines(sims[i, , 1] ~ times, lwd = 0.2)

# find zeros in the imputed daily case counts, so we can skip them
# this needs to happen here for some reason
known_cases <- (overseas + unknown_local + known_local) > 0
some_cases_known <- which(known_cases > 0)

# divide the reported cases by the reporting rates to get the expected total
# number of of cases in each source/state/time, then sum to get the expected
# total in each state/time
overseas_total_cases <- overseas[some_cases_known] / overseas_reporting[some_cases_known]
known_local_total_cases <- known_local[some_cases_known] / known_local_reporting[some_cases_known]
unknown_local_total_cases <- unknown_local[some_cases_known] / unknown_local_reporting[some_cases_known]
total_cases <- overseas_total_cases + known_local_total_cases + unknown_local_total_cases

# Distribution over plausible baseline CFR values from China study. The 95% CIs
# are symmetric around the estimate, so we assume it's a an approximately
# Gaussian distribution, truncated to allowable values. Define separately for each state
true_cfr_mean <- 1.38
true_cfr_sigma <- 0.077
baseline_cfr_perc <- normal(true_cfr_mean, true_cfr_sigma, truncation = c(0, 100))
log_baseline_cfr <- log(baseline_cfr_perc) - log(100)

# expected deaths each state and time
log_expected_deaths <- log(total_cases) + log_baseline_cfr

# sampling distribution
expected_deaths <- exp(log_expected_deaths)
observed_deaths <- deaths[some_cases_known]
distribution(observed_deaths) <- poisson(expected_deaths)

set.seed(2020-04-02)
m <- model(
  baseline_cfr_perc,
  unknown_local_gp_lengthscale,
  unknown_local_gp_sigma,
  unknown_local_sigma,
  known_local_gp_lengthscale,
  known_local_gp_sigma,
  known_local_sigma
)

n_chains <- 10

draws <- mcmc(
  m,
  chains = n_chains,
  n_samples = 1000,
  # initial_values = inits,
  one_by_one = TRUE
)

# check convergence before continuing
r_hats <- coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
n_eff <- coda::effectiveSize(draws)
max(r_hats)
min(n_eff)

# compute different ascertainment rates
# some of these will have non-finite values because of missing data (no cases with known outcomes)
# find a better way of computing them?
overseas_total_cases_all <- overseas / overseas_reporting
known_local_total_cases_all <- known_local / known_local_reporting
unknown_local_total_cases_all <- unknown_local / unknown_local_reporting
total_cases_all <- overseas_total_cases_all + known_local_total_cases_all + unknown_local_total_cases_all

# compute the ratio of total to observed cases in state and nationally, and trace these
total_cases_detected <- overseas + known_local + unknown_local

# overall reporting rate (all sources), by state and nationwide
combined_reporting <- total_cases_detected / total_cases_all
combined_reporting_national <- rowSums(total_cases_detected) / rowSums(total_cases_all)

# reporting rates for each source nationwide
unknown_local_reporting_national <- rowSums(unknown_local) / rowSums(unknown_local_total_cases_all)
known_local_reporting_national <- rowSums(known_local) / rowSums(known_local_total_cases_all)

unknown_local_reporting_national_draws <- calculate(unknown_local_reporting_national, values = draws)
unknown_local_reporting_national_draws_mat <- as.matrix(unknown_local_reporting_national_draws)
plot(colMeans(unknown_local_reporting_national_draws_mat) ~ times, type = "l", ylim = c(0, 1))
lines(apply(unknown_local_reporting_national_draws_mat, 1, quantile, 0.025), lty = 2)

png("reporting_rate_timeseries_by_state.png",
    width = 1200,
    height = 1800,
    pointsize = 30)
par(mfrow = c(4, 2),
    mar = c(5, 4, 4, 3))

state_names <- colnames(death_matrix)
for (i in seq_len(n_states)) {
  
  # subset to the time after the state's first case
  start <- which(cases_known_matrix[, i] > 0)[1]
  index <- start:n_times
  
  # predict each state's timeseries
  draws <- calculate(reporting_rate_smooth[index, i], values = draws)
  draws_mat <- as.matrix(draws)
  mean <- colMeans(draws_mat)
  ci <- apply(draws_mat, 2, quantile, c(0.025, 0.975))
  iqr <- apply(draws_mat, 2, quantile, c(0.25, 0.75))
  
  times_plot <- times[index]
  # subset times to when this state first saw a case
  
  plot(mean ~ times_plot, type = "n",
       ylim = c(0, 1),
       xlim = range(times),
       axes = FALSE,
       xlab = "date of symptomatic case report",
       ylab = "probability of detection")
  polygon(x = c(times_plot, rev(times_plot)),
          y = c(ci[1, ], rev(ci[2, ])),
          col = blues9[2], lty = 0)
  polygon(x = c(times_plot, rev(times_plot)),
          y = c(iqr[1, ], rev(iqr[2, ])),
          col = blues9[3], lty = 0)
  lines(mean ~ times_plot, lwd = 4, col = blues9[6])
  axis(2, las = 2)
  
  # subtract 13 days from dates to reflect the date at which symptomatic would
  # have been detected
  first_date <- min(aus_timeseries$date - 13)
  
  axis(1, at = inducing_points, labels = first_date + inducing_points - 1)
  
  
  title(main = state_names[i])
  abline(v = times_plot[1], lwd = 1.5, col = "red")
  text(x = times_plot[1], y = 1.08,
       labels = "first case reported",
       xpd = NA, col = "red", cex = 0.8)
  
  # add dates of recorded deaths in upper rug plot 
  death_times <- times_plot[death_matrix[times_plot, i] > 0]
  if (length(death_times) > 0) {
    rug(death_times, side = 1, lwd = 2)
    text(x = max(times_plot), y = -0.02,
         labels = "deaths",pos = 4,
         xpd = NA, cex = 0.8)
  }
  
}
dev.off()

# calculate latest estimates
draws <- calculate(reporting_rate_smooth[n_times, ], values = draws)
draws_mat <- as.matrix(draws)
colnames(draws_mat) <- state_names
library(bayesplot)
library(ggplot2)
bayesplot::color_scheme_set("blue")
bayesplot::mcmc_intervals(draws_mat, point_est = "mean", prob = 0.5, prob_outer = 0.95) +
  ggplot2::xlim(0, 1) +
  ggplot2::theme_minimal() +
  ggplot2::ggtitle(paste("estimated reporting rates for symptomatic cases\non",
                         max(aus_timeseries$date - 13),
                         "(the latest available data)"))
ggplot2::ggsave("latest_reporting_rates.png", scale = 1)

reporting_rate <- data.frame(mean = colMeans(draws_mat),
                             lower_50 = apply(draws_mat, 2, quantile, 0.25),
                             upper_50 = apply(draws_mat, 2, quantile, 0.75),
                             lower_95 = apply(draws_mat, 2, quantile, 0.025),
                             upper_95 = apply(draws_mat, 2, quantile, 0.975))

write.csv(reporting_rate,
          "latest_reporting_rates.csv",
          row.names = FALSE)
knitr::kable(round(reporting_rate, 3))

# compute a national estimate
weights <- colSums(cases_known_matrix)
weights <- weights / sum(weights)
national_reporting_rate_estimate <- reporting_rate_smooth %*% as.matrix(weights)
draws_national <- calculate(national_reporting_rate_estimate[n_times], values = draws)
draws_national_vec <- as.matrix(draws_national)[, 1]
c(mean = mean(draws_national_vec), quantile(draws_national_vec, c(0.25, 0.75, 0.025, 0.975)))

