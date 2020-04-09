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
  
  # create a matrix of the correction factors to subset
  blank <- zeros(nrow(donor), ncol(donor))
  correction_mat <- sweep(blank, 2, correction, FUN = "+")
  
  recipient[missing] <- donor[missing] * correction_mat[missing]
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
  
  # mean GP
  mu <- temporal_gp(
    times = times,
    lengthscale = lengthscale_mu,
    sigma = sigma_mu,
    sigma_intercept = sigma_mu_intercept,
    n_inducing = n_inducing,
    tol = tol
  )
  
  # GPs for the states deviations (manually defined because multiple GPs from
  # the same kernel isn't yet possible in greta.gp)
  kernel_z <-
    rbf(lengthscales = lengthscale_z, variance = sigma_z ^ 2) + 
    bias(sigma_z_intercept ^ 2)
  
  inducing_points <- seq(min(times), max(times), length.out = n_inducing + 1)[-1]
  
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

# given date-by-state matrices of the observed number of cases and the reporting
# rates, compute a national average of the reporting rate
national_average <- function (reporting_rate, observed_cases) {
  
  # compute state weights from timeseries totals on the number of observed cases
  total_cases <- observed_cases / reporting_rate
  total_cases_sum <- colSums(total_cases)
  state_weights <- total_cases_sum / sum(total_cases_sum)
  
  # compute weights only where there are some observed cases, so we can compute
  # row sums in state_weights
  weighted <- sweep(reporting_rate, 2, state_weights, FUN = "*")
  rowSums(weighted)
  
}

# given a greta array, return the posterior mean, 95% credible interval, and
# posterior interquartile range
summarise_samples <- function(ts, draws) {
  draws <- calculate(ts, values = draws)
  draws_mat <- as.matrix(draws)
  list(
    mean = colMeans(draws_mat),
    ci = apply(draws_mat, 2, quantile, c(0.025, 0.975)),
    iqr = apply(draws_mat, 2, quantile, c(0.25, 0.75))
  )
}

plot_timeseries <- function(summary, title = NULL, death_times = c()) {
  
  mean <- summary$mean
  ci <- summary$ci
  iqr <- summary$iqr
  
  plot(mean ~ times, type = "n",
       ylim = c(0, 1),
       xlim = range(times),
       axes = FALSE,
       xlab = "date of symptomatic case report",
       ylab = "probability of detection")
  polygon(x = c(times, rev(times)),
          y = c(ci[1, ], rev(ci[2, ])),
          col = blues9[2], lty = 0)
  polygon(x = c(times, rev(times)),
          y = c(iqr[1, ], rev(iqr[2, ])),
          col = blues9[3], lty = 0)
  lines(mean ~ times, lwd = 4, col = blues9[6])
  axis(2, las = 2)
  
  # subtract 13 days from dates to reflect the date at which symptomatic would
  # have been detected
  first_date <- min(aus_timeseries$date - 13)
  times_axis <- seq(1, n_times, length.out = 5)
  axis(1, at = times_axis, labels = first_date + times_axis - 1,
       cex.axis = 0.8)
  
  title(main = title)
  
  # add dates of recorded deaths in lower rug plot 
  
  if (length(death_times) > 0) {
    rug(death_times, side = 1, lwd = 1)
    text(x = max(times), y = -0.02,
         labels = "deaths",pos = 4,
         xpd = NA, cex = 0.8)
  }
  
}

format_summary <- function (summary) {
  df <- data.frame(
    state = colnames(deaths),
    mean = summary$mean,
    lower_50 = summary$iqr[1, ],
    upper_50 = summary$iqr[2, ],
    lower_95 = summary$ci[1, ],
    upper_95 = summary$ci[2, ]
  )
  rownames(df) <- NULL
  df
}

forest_plot <- function (summary, title, labels = TRUE) {
  
  df <- format_summary(summary)
  
  require(ggplot2)
  
  p <- ggplot(data = df, aes(x = state, y = mean, ymin = lower_95, ymax = upper_95)) +
    geom_linerange(col = blues9[4], size = 1.5) +
    geom_linerange(aes(ymin = lower_50, ymax = upper_50), size = 2.5, col =  blues9[6]) +
    geom_point(col = blues9[9], size = 3) + 
    coord_flip() +
    scale_x_discrete(limits = rev(levels(df$state))) +
    ylim(0, 1) +
    ylab("probability of detection") +
    xlab("") +
    ggtitle(title[[1]], subtitle = title[[2]]) + 
    geom_hline(yintercept = c(0, 1), col = grey(0.5)) +
    theme_minimal()
  
  if(!labels) {
    p <- p + theme(axis.text.y = element_blank())
  }
  
  p
  
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

n_states <- ncol(deaths)
n_times <- nrow(deaths)
times <- seq_len(n_times)


# estimating those parameters directly as with this model, the numbers are so
# large there's essentially no uncertainty on them. Could relax the model a bit
# later (e.g. hierarchical model by state on the proportions), and infer them

other_is_overseas_sigma <- normal(0, 0.5, truncation = c(0, Inf))
other_is_overseas_mean <- normal(0, sqrt(0.5))
other_is_overseas_raw <- normal(0, 1, dim = n_states)
other_is_overseas_z <- other_is_overseas_mean + other_is_overseas_raw * other_is_overseas_sigma
p_other_is_overseas <- iprobit(other_is_overseas_z)

# # check prior on probability (slightly convex is fine, strong 'U' shape not so
# # good)
# hist(calculate(p_other_is_overseas[1], nsim = 10000)[[1]], breaks = 100)

all_is_unknown_local_sigma <- normal(0, 0.5, truncation = c(0, Inf))
all_is_unknown_local_mean <- normal(0, sqrt(0.5))
all_is_unknown_local_raw <- normal(0, 1, dim = n_states)
all_is_unknown_local_z <- all_is_unknown_local_mean + all_is_unknown_local_raw * all_is_unknown_local_sigma
p_all_is_unknown_local <- iprobit(all_is_unknown_local_z)

# get data to inform these, skipping some states for whichthere is no data
state_names <- tibble::tibble(state = colnames(deaths))

overseas_vs_other <- aus_timeseries %>%
  select(overseas, known_local) %>%
  na.omit() %>%
  mutate(other = overseas + known_local) %>%
  summarise(successes = sum(overseas),
            trials = sum(other)) %>%
  right_join(state_names)

not_missing <- which(!is.na(overseas_vs_other$trials))
distribution(overseas_vs_other$successes[not_missing]) <-
  binomial(overseas_vs_other$trials[not_missing],
           p_other_is_overseas[not_missing])

unknown_local_vs_all <- aus_timeseries %>%
  select(unknown_local, cases) %>%
  na.omit() %>%
  summarise(successes = sum(unknown_local),
            trials = sum(cases)) %>%
  right_join(state_names)

not_missing <- which(!is.na(unknown_local_vs_all$trials))
distribution(unknown_local_vs_all$successes[not_missing]) <-
  binomial(unknown_local_vs_all$trials[not_missing],
           p_all_is_unknown_local[not_missing])

# fill in values
unknown_local <- impute_values(unknown_local, cases, p_all_is_unknown_local)
other <- impute_values(other, cases, 1 - p_all_is_unknown_local)
overseas <- impute_values(overseas, other, p_other_is_overseas)
known_local <- impute_values(known_local, other, 1 - p_other_is_overseas)

# now disaggregate these cases using the delay distribution
unknown_local <- cases_known_outcome_matrix(unknown_local)
known_local <- cases_known_outcome_matrix(known_local)
overseas <- cases_known_outcome_matrix(overseas)

# define hierarchical probit-GPs for reporting rates of unknown local, and known
# local cases. Assume reporting rate for overseas-acquired cases is perfect.

unknown_local_gp_lengthscale <- lognormal(4, 0.5)
unknown_local_gp_sigma <- lognormal(-2, 0.5)
unknown_local_sigma <- normal(0, 0.5, truncation = c(0, Inf))
unknown_local_intercepts_raw <- normal(0, 1, dim = c(1, n_states))
unknown_local_intercepts <- unknown_local_intercepts_raw * unknown_local_sigma
unknown_local_gp <- temporal_gp(
  times = times,
  lengthscale = unknown_local_gp_lengthscale,
  sigma = unknown_local_gp_sigma,
  tol = 0
)
unknown_local_z <- kronecker(
  unknown_local_gp,
  unknown_local_intercepts,
  FUN = "+"
)
unknown_local_reporting <- iprobit(unknown_local_z)

# model the difference between these two, to decorrelate them, and so that
# shrinkage is towards no difference
known_diff_gp_lengthscale <- lognormal(4, 0.5)
known_diff_gp_sigma <- lognormal(-2, 0.5)
known_diff_sigma <- normal(0, 0.5, truncation = c(0, Inf))
known_diff_intercepts_raw <- normal(0, 1, dim = n_states)
known_diff_intercepts <- known_diff_intercepts_raw * known_diff_sigma
known_diff_gp <- temporal_gp(
  times = times,
  lengthscale = known_diff_gp_lengthscale,
  sigma = known_diff_gp_sigma,
  tol = 0
)
known_diff_z <- kronecker(
  known_diff_gp,
  t(known_diff_intercepts),
  FUN = "+"
)

# combine them to get probit-reporting rate for known locals
known_local_z <- unknown_local_z + known_diff_z
known_local_reporting <- iprobit(known_local_z)

overseas_reporting <- ones(n_times, n_states)

# # visualise priors:
# nsim <- 10000
# sims <- calculate(unknown_local_reporting[n_times, 1],
#                   nsim = nsim)[[1]]
# plot(sims[1, , 1] ~  times, type = "n", ylim = c(0, 1))
# for(i in seq_len(nsim)) lines(sims[i, , 1] ~ times, lwd = 0.2)

# find zeros in the imputed daily case counts, so we can skip them
known_cases <- (overseas + unknown_local + known_local)
known_cases_vals <- calculate(known_cases,
                              values = list(
                                p_all_is_unknown_local = rep(0.5, n_states),
                                p_other_is_overseas = rep(0.5, n_states)
                              ))[[1]]
some_cases_known <- which(known_cases_vals > 0)

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
  known_diff_gp_lengthscale,
  known_diff_gp_sigma,
  known_diff_sigma
)

n_chains <- 10

draws <- mcmc(
  m,
  chains = n_chains,
  n_samples = 1000,
  one_by_one = TRUE
)

# check convergence before continuing
r_hats <- coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
n_eff <- coda::effectiveSize(draws)
max(r_hats)
min(n_eff)

# compute different ascertainment rate from all sources, and compute a national aggregate
# want:
#  - national-level timeseries by source (including 'all')
#  - latest point estimates by state and source (including 'all')

# compute weights based on timeseries totals not contemporary totals, since the
# timeseries totals have too many zeros and are causing all sorts of numerical
# issues

# reported cases from all sources
all_sources <- overseas + known_local + unknown_local

# true number of cases from each and all sources
overseas_total <- colSums(overseas / overseas_reporting)
known_local_total <- colSums(known_local / known_local_reporting)
unknown_local_total <- colSums(unknown_local / unknown_local_reporting)
all_sources_total <- overseas_total + known_local_total + unknown_local_total

# fraction of true number of cases from each source, within each state
overseas_weight <- overseas_total / all_sources_total
known_local_weight <- known_local_total / all_sources_total
unknown_local_weight <- unknown_local_total / all_sources_total 

# reporting rate for all sources
all_sources_reporting <-
  sweep(overseas_reporting, 2, overseas_weight, FUN = "*") + 
  sweep(known_local_reporting, 2, known_local_weight, FUN = "*") + 
  sweep(unknown_local_reporting, 2, unknown_local_weight, FUN = "*")

# timeseries of national averages for each source
unknown_local_reporting_nat <- national_average(unknown_local_reporting, unknown_local)
known_local_reporting_nat <- national_average(known_local_reporting, known_local)
all_sources_reporting_nat <- national_average(all_sources_reporting, all_sources)

# contemporary estimates of reporting rates in each state, and by each source
unknown_local_reporting_latest <- unknown_local_reporting[n_times, ]
known_local_reporting_latest <- known_local_reporting[n_times, ]
all_sources_reporting_latest <- all_sources_reporting[n_times, ]

# report these in three side-by-side forest plots, using the code below

# predict each national aggregate timeseries for each source
all_sources_nat_sry <- summarise_samples(all_sources_reporting_nat, draws)
known_local_nat_sry <- summarise_samples(known_local_reporting_nat, draws)
unknown_local_nat_sry <- summarise_samples(unknown_local_reporting_nat, draws)

# vector of death times for plotting
death_vec <- rowSums(deaths)
death_times <- times[death_vec > 0]
death_counts <- death_vec[death_times]
death_jitter <- jitter(rep(death_times, death_counts))


# plot timeseries for two sources and overall in forest plots
png("reporting_rate_timeseries_national.png",
    width = 1500,
    height = 500,
    pointsize = 30)
par(mfrow = c(1, 3),
    mar = c(5, 4, 4, 3))

plot_timeseries(all_sources_nat_sry, "all cases", death_jitter)
plot_timeseries(known_local_nat_sry, "local transmission - known source", death_jitter)
plot_timeseries(unknown_local_nat_sry, "local transmission - unknown source", death_jitter)

dev.off()

# posterior summaries of the latest rates by stransmission source for each state
all_sources_latest_sry <- summarise_samples(all_sources_reporting_latest, draws)
known_local_latest_sry <- summarise_samples(known_local_reporting_latest, draws)
unknown_local_latest_sry <- summarise_samples(unknown_local_reporting_latest, draws)

p1 <- forest_plot(all_sources_latest_sry,
                  list("all cases", waiver()))
p2 <- forest_plot(known_local_latest_sry,
                  list("local transmission", "known source"),
                  labels = FALSE)
p3 <- forest_plot(unknown_local_latest_sry,
              list("local transmission", "unknown source"),
              labels = FALSE)

p1 + p2 + p3 + plot_annotation(title = paste("reporting rates as of", max(aus_timeseries$date - 13)))
ggsave("latest_reporting_rates_by_source.png")


# format outputs
all_sources_df <- format_summary(all_sources_latest_sry)
known_local_df <- format_summary(known_local_latest_sry)
unknown_local_df <- format_summary(unknown_local_latest_sry)

df <- rbind(
  cbind(source = "all cases", all_sources_df),
  cbind(source = "known local", known_local_df),
  cbind(source = "unknown local", unknown_local_df)
)

write.csv(df,
          "latest_reporting_rates_by_source.csv",
          row.names = FALSE)
knitr::kable(df, digits = 3)

# get national latest figure

latest_national <- c(all_sources_reporting_nat[n_times], known_local_reporting_nat[n_times], unknown_local_reporting_nat[n_times])
sry <- summarise_samples(latest_national, draws)
nat_df <- data.frame(
  source = c("all cases", "known local", "unknown local"),
  mean = sry$mean,
  lower_50 = sry$iqr[1, ],
  upper_50 = sry$iqr[2, ],
  lower_95 = sry$ci[1, ],
  upper_95 = sry$ci[2, ]
)
rownames(nat_df) <- NULL
write.csv(nat_df,
          "latest_reporting_rates_by_source_national.csv",
          row.names = FALSE)
knitr::kable(nat_df, digits = 3)
