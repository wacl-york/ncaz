library(tidyverse)
library(lubridate)
library(openair)
library(KFAS)
source("utils.R")

OUTPUT_DIR <- OUTPUT_DIR_FROM_CRON

# Load states
states <- lapply(SITES, function(x) {
  readRDS(sprintf("%s/states/univariate_%s.rds",
                  OUTPUT_DIR,
                  x$human_readable))
})
models <- lapply(SITES, function(x) {
  readRDS(sprintf("%s/models/univariate_%s.rds",
                  OUTPUT_DIR,
                  x$human_readable))
})

# for each site:
for (site in names(SITES)) {
  
  # get most recent date and data to update with
  this_states <- states[[site]]
  last_update <- this_states[[length(this_states)]]$date
  this_df <- load_data(site) |>
                filter(time > last_update)

  if (nrow(this_df) == 0) next
  
  # update states
  a_last <- this_states[[length(this_states)]]$a
  P_last <- this_states[[length(this_states)]]$P
  detrended <- array(dim=nrow(this_df))
  for (i in 1:nrow(this_df)) {
    updated <- update_univariate(models[[site]], this_df[i, ], a_last, P_last)
    a_last <- updated$a
    P_last <- updated$P
    this_states <- append(this_states, list(list(
      date=this_df$time[i],
      a=a_last,
      P=P_last
    )))
    detrended[i] <- exp(a_last['level', ])
  }
  
  # Save states + results separately
  saveRDS(this_states, sprintf("%s/states/univariate_%s.rds", OUTPUT_DIR, SITES[[site]]$human_readable))
  
  df_to_save <- this_df |>
                  mutate(code = site,
                         detrended = as.numeric(detrended)) |>
                  select(time, code, no2, detrended)
  write_csv(df_to_save, sprintf("%s/data/results.csv", OUTPUT_DIR), append = TRUE)
} 

