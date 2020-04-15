# greta modelling functions

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
