# plotting and outputting functions
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