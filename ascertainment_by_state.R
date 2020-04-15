# Time-series estimates of reporting rates for Australian states

# load misc functions
source("functions/load_data.R")
source("functions/modelling_functions.R")
source("functions/output_functions.R")

# get timeseries of cases and deaths, by state and - where available - source
aus_timeseries <- load_aus_timeseries() %>%
  mutate(cases_known_outcome = cases_known_outcome(cases))

# deaths and all cases (as in Johns Hopkins data)
death_matrix <- date_state_matrix("deaths", aus_timeseries)
cases_known_matrix <- date_state_matrix("cases_known_outcome", aus_timeseries)

# build model for contemporary observed number of deaths and expected number of
# deaths:
#   deaths_t ~ Poisson(expected_deaths_t)
#   expected_deaths_t = cases_known_outcomes_t * CFR_t
#   CFR_t = baseline_CFR / reporting_rate_t

library(greta.gp)

n_states <- ncol(death_matrix)
n_times <- nrow(death_matrix)

# a timeseries of reporting rates for each state, modelled as hierarchical Gaussian processes

# squared-exponential GP kernels (plus intercepts) with unknown parameters for
# the national and state-level processes. lognormal prior for lengthscales to
# reduce prior probability hof high temporal change
national_lengthscale <- lognormal(4, 0.5)
national_sigma <- lognormal(-1, 1)
national_temporal <- rbf(lengthscales = national_lengthscale,
                         variance = national_sigma ^ 2)
national_intercept <- bias(0.5)
national_kernel <- national_intercept + national_temporal

state_lengthscale <- lognormal(4, 0.5)
state_sigma <- lognormal(-1, 1)
state_temporal <- rbf(lengthscales = state_lengthscale,
                       variance = state_sigma ^ 2)
state_intercept <- bias(0.5)
state_kernel <- state_intercept + state_temporal

# IID gaussian kernel to represent observation error (overdispersion)
sigma_obs <- normal(0, 0.5, truncation = c(0, Inf))
observation_kernel <- white(sigma_obs ^ 2)

state_observed_kernel <- state_kernel + observation_kernel

# a set of inducing points at which to estimate the GPs (subset of regressors
# approximation)
# put an inducing point on the last time point (most recent date), but otherwise
# space them out
times <- seq_len(n_times)
n_inducing <- 5
inducing_points <- seq(min(times), max(times), length.out = n_inducing + 1)[-1]

# GP for the national mean effect - add jitter to help with matrix inversion
tol <- 1e-6
mu <- greta.gp::gp(times, inducing = inducing_points, national_kernel, tol = tol)

# GPs for the state deviations (manually defined as multiple GPs at once isn't
# yet possible in greta.gp)
v <- normal(0, 1, dim = c(n_inducing, n_states))
Kmm <- state_observed_kernel(inducing_points)
Lm <- t(chol(Kmm))
Kmn <- state_observed_kernel(inducing_points, times)
A <- forwardsolve(Lm, Kmn)
z_state <- t(A) %*% v

# add the mean effect on
z <- sweep(z_state, 1, mu, "+")

# convert to probabilities
reporting_rate <- iprobit(z)

# # visualise prior:
# nsim <- 300
# sims <- calculate(reporting_rate, nsim = nsim)[[1]]
# plot(sims[1, , 1] ~  times, type = "n", ylim = c(0, 1))
# for(i in seq_len(nsim)) lines(sims[i, , 1] ~ times, lwd = 0.2)


# Distribution over plausible baseline CFR values from China study. The 95% CIs
# are symmetric around the estimate, so we assume it's a an approximately
# Gaussian distribution, truncated to allowable values. 
true_cfr_mean <- 1.38
true_cfr_sigma <- 0.077
baseline_cfr_perc <- normal(true_cfr_mean, true_cfr_sigma, dim = n_states, truncation = c(0, 100))

# compute CFR for each state and time
log_baseline_cfr <- log(baseline_cfr_perc) - log(100)
log_reporting_rate <- log(reporting_rate)
log_cfr <- sweep(-log_reporting_rate, 2, log_baseline_cfr, "+")

# define sampling distribution, subsetting to where cases_known_matrix > 0
# do the exponentiation down here, so greta can use log version in poisson density
some_cases_known <- which(cases_known_matrix > 0)

log_expected_deaths <- log_cfr[some_cases_known] +
  log(cases_known_matrix)[some_cases_known]
expected_deaths <- exp(log_expected_deaths)
observed_deaths <- death_matrix[some_cases_known]
distribution(observed_deaths) <- poisson(expected_deaths)

set.seed(2020-04-02)
m <- model(reporting_rate)

n_chains <- 50

inits <- replicate(
  n_chains,
  initials(
    national_lengthscale = rlnorm(1, 4, 0.5),
    national_sigma = rlnorm(1, -1, 1),
    state_lengthscale = rlnorm(1, 4, 0.5),
    state_sigma = rlnorm(1, -1, 1),
    baseline_cfr_perc = max(0.001, min(99.999,
                                       rnorm(1, true_cfr_mean, true_cfr_sigma)
    ))
  ),
  simplify = FALSE
)

draws <- mcmc(
  m,
  sampler = hmc(Lmin = 15, Lmax = 20),
  chains = n_chains,
  n_samples = 1000,
  one_by_one = TRUE
)

# check convergence before continuing
r_hats <- coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
n_eff <- coda::effectiveSize(draws)
max(r_hats)
min(n_eff)

# now compute reporting rates *without* additional observation error (possible
# overdispersion due to clumped death rates)
Kmm_smooth <- state_kernel(inducing_points)
Lm_smooth <- t(chol(Kmm_smooth))
Kmn_smooth <- state_kernel(inducing_points, times)
A_smooth <- forwardsolve(Lm_smooth, Kmn_smooth)
z_state_smooth <- t(A_smooth) %*% v
z_smooth <- sweep(z_state_smooth, 1, mu, "+")
reporting_rate_smooth <- iprobit(z_smooth)


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

