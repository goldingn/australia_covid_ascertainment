# load and format public data on cases in each state

load_aus_timeseries <- function() {
  
  # load data from covid19data.org.au (scraped by Chris Baker) with cases by
  # source for some states.
  source_data <- load_cases_by_source_and_state()
  
  # subset the case data to Aus states
  aus <- load_cases_and_deaths_by_state()
  
  # expand the type of count for each state, attach the sources where known, and
  # compute cases with known outcomes
  aus %>%
    tidyr::pivot_wider(names_from = type, values_from = count) %>%
    select(-recovered, cases = confirmed, deaths = death) %>%
    left_join(source_data) %>%
    mutate_at(c("cases", "deaths", "overseas", "known_local", "unknown_local", "other"),
              ~ pmax(., 0))
  
}

# download the latest Johns Hopkins data (it's in the package directly, so need
# to reinstall regularly) and subset/format to Australia
load_cases_and_deaths_by_state <- function() {
  
  require(dplyr)
  
  remotes::install_github("RamiKrispin/coronavirus", upgrade = "never")
  coronavirus::coronavirus %>%
    filter(Country.Region == "Australia") %>%
    transmute(state = Province.State, date, count = cases, type) %>%
    group_by(state)
}


load_cases_by_source_and_state <- function() {
  
  require(processx)
  require(dplyr)
  require(tibble)
  require(purrr)
  
  # use Chris' git repo (included as a submodule) to download the data
  python <- switch(Sys.info()["user"],
                   nick = "/Users/nick/miniconda3/bin/python",
                   "python")

  # scrape data into the submodule
  processx::run(python, "scrape_data.py", wd = "covid19data")
  
  # find the csv files
  files <- list.files("covid19data/", pattern = ".csv$", full.names = TRUE)

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


