# Time-series estimates of reporting rates for Australian states

# load misc functions
source("load_data.R")
source("modelling_functions.R")
source("output_functions.R")

# get timeseries of cases and deaths, by state and - where available - source
aus_timeseries <- load_aus_timeseries()
  
# get date-by-state matrices for daily case counts by state for each source
date_state_matrices <- get_date_state_matrices(aus_timeseries)

# unpack them
deaths <- date_state_matrices$deaths
unknown_local_cases <- date_state_matrices$unknown_local_cases_imputed
known_local_cases <- date_state_matrices$known_local_cases_imputed
overseas_cases <- date_state_matrices$overseas_cases_imputed
p_all_is_unknown_local <- date_state_matrices$p_all_is_unknown_local
p_other_is_overseas <- date_state_matrices$p_other_is_overseas

# now disaggregate these cases using the delay distribution
unknown_local <- cases_known_outcome_matrix(unknown_local_cases)
known_local <- cases_known_outcome_matrix(known_local_cases)
overseas <- cases_known_outcome_matrix(overseas_cases)

# define hierarchical probit-GPs for reporting rates of unknown local, and known
# local cases. Assume reporting rate for overseas-acquired cases is perfect.
library(greta.gp)

n_states <- ncol(deaths)
n_times <- nrow(deaths)
times <- seq_len(n_times)

unknown_local_gp_lengthscale <- lognormal(4, 0.5)
unknown_local_gp_sigma <- lognormal(-2, 0.5)
unknown_local_intercepts <- hierarchical_normal(dim = n_states, mean_sd = NULL)
unknown_local_gp <- temporal_gp(
  times = times,
  lengthscale = unknown_local_gp_lengthscale,
  sigma = unknown_local_gp_sigma,
  tol = 0
)
unknown_local_z <- kronecker(
  unknown_local_gp,
  t(unknown_local_intercepts),
  FUN = "+"
)
unknown_local_reporting <- iprobit(unknown_local_z)

# model the difference between these two, to decorrelate them, and so that
# shrinkage is towards no difference
known_diff_gp_lengthscale <- lognormal(4, 0.5)
known_diff_gp_sigma <- lognormal(-2, 0.5)
known_diff_intercepts <- hierarchical_normal(dim = n_states, mean_sd = NULL)
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
  known_diff_gp_lengthscale,
  known_diff_gp_sigma
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
