# load and format data on timeseries of cases by source, in each state
load_data <- function () {
  
  require(gert)
  require(dplyr)
  require(tibble)
  require(purrr)
  
  # clone Chris' git repo to download the data

  # create a temporary directory to download data, and delete on exit
  dir.create(tmp <- tempfile())
  on.exit(unlink(tmp))
  gert::git_clone("https://github.com/cmbaker00/covid19data.git", path = tmp)
  
  # find the csv files
  files <- list.files(tmp, pattern = ".csv$", full.names = TRUE)

  # compile information on the datasets
  datasets <- tibble::as_tibble(
    rbind(
      c(state = "Australian Capital Territory",
        file = "Daily confirmed cases in ACT by transmission source.csv",
        cumulative = FALSE),
      c(state = "New South Wales",
        file = "Daily confirmed cases in NSW by transmission source.csv",
        cumulative = FALSE),
      c(state = "Northern Territory",
        file = "Daily confirmed cases in NT by transmission source.csv",
        cumulative = FALSE),
      c(state = "South Australia",
        file = "SA Cumulative view of transmission sources over time.csv",
        cumulative = TRUE),
      c(state = "Tasmania",
        file = "Daily confirmed cases in Tasmania by transmission source.csv",
        cumulative = FALSE),
      c(state = "Victoria",
        file = "Daily confirmed cases in Victoria by transmission source.csv",
        cumulative = FALSE),
      c(state = "Western Australia",
        file = "WA Cumulative view of unknown local and other transmission sources over time.csv",
        cumulative = TRUE)
    )
  )
  
  # load and format them all
  datasets %>%
    mutate(cumulative = as.logical(cumulative)) %>%
    left_join(tibble::tibble(file = basename(files), path = files)) %>%
    select(state, path, cumulative) %>%
    purrr::pmap(load_dataset) %>%
    bind_rows()
  
}

load_dataset <- function (state, path, cumulative) {

  require(lubridate)
  
  vals <- read.csv(path,
                   stringsAsFactors = FALSE)
  names <- names(vals)
  names[1] <- "date"
  names <- gsub("Overseas...Interstate", "Overseas", names)
  names(vals) <- names
  if (!"Other" %in% names) {
    vals$Other <- 0
  }
  
  data <- vals %>%
    mutate(date = paste0(date, "/2020"),
           date = lubridate::dmy(date),
           overseas = Overseas,
           known_local = Known.Local,
           unknown_local = Unknown.Local,
           other = Other) %>%
    arrange(date) %>%
    select(date,
           overseas,
           known_local,
           unknown_local,
           other)    
  
  # use a diff function that preserves first element, then drops it
  if (cumulative) {
    data <- data %>%
      mutate_at(vars(!matches("date")), difference) %>%
      filter(date > min(date))
  }
  
  # add on column for state name
  data$state <- state
  data
}

difference <- function (x) {
  c(0, diff(x))
}

# horrifyingly long working to get ballpark translation to detection rate for local transmissions

# total_rate = reported / all 
# reported = local_rate * local + import
# all = local + import
# 
# total_rate = (local_rate * local + import) / (local + import)
# local_apparent = local_rate * local
# local = local_apparent / local_rate
# 
# total_rate = (local_apparent + import) / (import + local_apparent / local_rate)
# all_apparent = local_apparent + import
# 
# total_rate = all_apparent / (import + local_apparent / local_rate)
# import + local_apparent / local_rate = all_apparent / total_rate
# local_apparent / local_rate = all_apparent / total_rate - import
# 
# local_rate = local_apparent / (all_apparent / total_rate - import)
# import = all_apparent - local_apparent
# 
# local_rate = local_apparent / (all_apparent / total_rate - all_apparent + local_apparent)
# local_rate = local_apparent / (all_apparent / total_rate - total_rate * all_apparent / total_rate + total_rate * local_apparent / total_rate)
# local_rate = local_apparent / ((all_apparent - total_rate * all_apparent + total_rate * local_apparent)  / total_rate)
# local_rate = local_apparent / (all_apparent * (1 - total_rate + total_rate * local_apparent / all_apparent)  / total_rate)
# 
# prop_apparent_local <- local_apparent / all_apparent  
# local_rate = local_apparent / (all_apparent * (1 - total_rate + total_rate * prop_apparent_local)  / total_rate)
# 
# local_rate = total_rate * local_apparent / (all_apparent * (1 - total_rate + total_rate * prop_apparent_local))
# local_rate = total_rate * prop_apparent_local / (1 - total_rate + total_rate * prop_apparent_local)
# 
# local_rate = prop_apparent_local / (prop_apparent_local - 1 + 1 / total_rate)
# local_rate = prop_apparent_local / (prop_apparent_local + (1 - total_rate) / total_rate)


prop_apparent_local <- 1/3
total_rate <- 0.93
(local_rate <- prop_apparent_local / (prop_apparent_local + (1 - total_rate) / total_rate))

local_apparent = 10
all_apparent = 30
import = 20


