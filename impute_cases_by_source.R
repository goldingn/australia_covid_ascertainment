# Imputing numbers of cases, split by transmission source, in Australian states

# load misc functions
source("load_data.R")
source("modelling_functions.R")
source("output_functions.R")

# get timeseries of cases and deaths, by state and - where available - source
aus_timeseries <- load_aus_timeseries()

# get date-by-state matrices for daily case counts by state for each source (does imputation)
date_state_matrices <- get_date_state_matrices(aus_timeseries)

# extract the things we care about
unknown_local_cases <- date_state_matrices$unknown_local_cases_imputed
known_local_cases <- date_state_matrices$known_local_cases_imputed
overseas_cases <- date_state_matrices$overseas_cases_imputed
p_all_is_unknown_local <- date_state_matrices$p_all_is_unknown_local
p_other_is_overseas <- date_state_matrices$p_other_is_overseas

# fit the model
library(greta)
m <- model(p_all_is_unknown_local, p_other_is_overseas)
draws <- mcmc(m, chains = 10, n_samples = 3000)

# check convergence before continuing
r_hats <- coda::gelman.diag(draws, autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
n_eff <- coda::effectiveSize(draws)
max(r_hats)
min(n_eff)

# # extract posterior samples of the *expected* number of cases
# unknown_local_cases_expected_sims <- get_date_state_sims(unknown_local_cases, draws)
# known_local_cases_expected_sims <- get_date_state_sims(known_local_cases, draws)
# overseas_cases_expected_sims <- get_date_state_sims(overseas_cases, draws)
# 
# # reset the counts that were known (should be the same, but this sets dimnames too)
# unknown_local_cases_expected_sims <- set_known_counts(
#   unknown_local_cases_expected_sims,
#   date_state_matrices$unknown_local_cases_raw
# )
# known_local_cases_expected_sims <- set_known_counts(
#   known_local_cases_expected_sims,
#   date_state_matrices$known_local_cases_raw
# )
# overseas_cases_expected_sims <- set_known_counts(
#   overseas_cases_expected_sims,
#   date_state_matrices$overseas_cases_raw
# )

# extract samples of the integer number of cases, conditioned on the posterior
# for the expectation

# flatten to vectors for multinomial draws
unknown_local_cases_vec <- c(unknown_local_cases)
known_local_cases_vec <- c(known_local_cases)
overseas_cases_vec <- c(overseas_cases)
probs <- cbind(unknown_local_cases_vec, known_local_cases_vec, overseas_cases_vec)
total_cases_vec <- c(date_state_matrices$cases)
cases_discrete <- multinomial(total_cases_vec, probs)

# reshape to matrices
n_times <- nrow(date_state_matrices$deaths)
n_states <- ncol(date_state_matrices$deaths)
unknown_local_discrete <- cases_discrete[, 1]
dim(unknown_local_discrete) <- c(n_times, n_states)
known_local_discrete <- cases_discrete[, 2]
dim(known_local_discrete) <- c(n_times, n_states)
overseas_discrete <- cases_discrete[, 3]
dim(overseas_discrete) <- c(n_times, n_states)

# get simulations
unknown_local_cases_discrete_sims <- get_date_state_sims(
  unknown_local_discrete,
  draws,
  simulate = TRUE
)
known_local_cases_discrete_sims <- get_date_state_sims(
  known_local_discrete,
  draws,
  simulate = TRUE
)
overseas_cases_discrete_sims <- get_date_state_sims(
  overseas_discrete,
  draws,
  simulate = TRUE
)

# reset the counts that were known (remove multinomial variation)
unknown_local_cases_discrete_sims <- set_known_counts(
  unknown_local_cases_discrete_sims,
  date_state_matrices$unknown_local_cases_raw
)
known_local_cases_discrete_sims <- set_known_counts(
  known_local_cases_discrete_sims,
  date_state_matrices$known_local_cases_raw
)
overseas_cases_discrete_sims <- set_known_counts(
  overseas_cases_discrete_sims,
  date_state_matrices$overseas_cases_raw
)

