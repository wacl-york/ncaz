library(tidyverse)
library(lubridate)
library(openair)
library(KFAS)
source("utils.R")

# get 2023 data from AURN and pre-process it accordingly (daily average correctly handling wind)
df <- importAURN(site = c("NEWC", "NCA3"), year=2023)
df_daily <- df |>
  mutate(time = as_datetime(as_date(date)),
         wind_y = ws * sin(2 * pi * wd / 360),
         wind_x = ws * cos(2 * pi * wd / 360)) |>
  group_by(code, time) |>
  summarise(no2 = mean(no2, na.rm=T),
            air_temp = mean(air_temp,na.rm=T),
            wind_y = mean(wind_y, na.rm=T),
            wind_x = mean(wind_x, na.rm=T),
            ws = mean(ws, na.rm=T)) |>
  ungroup() |>
  # Convert vector wind back to polar
  mutate(
      wd = as.vector(atan2(wind_y, wind_x) * 360 / 2 / pi),
      wd = ifelse(wd < 0, wd + 360, wd),
      ws_vec = (wind_x ^ 2 + wind_y ^ 2) ^ 0.5,
      wind_y_nows = sin(2 * pi * wd / 360),
      wind_x_nows = cos(2 * pi * wd / 360),
      intervention=as.numeric(time >= INTERVENTION_DATE)
)

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
  this_df <- df_daily |>
                filter(code == site,
                       time > last_update)
  
  # update states
  a_last <- this_states[[length(this_states)]]$a
  P_last <- this_states[[length(this_states)]]$P
  for (i in 1:nrow(this_df)) {
    updated <- update_univariate(models[[site]], this_df[i, ], a_last, P_last)
    a_last <- updated$a
    P_last <- updated$P
    this_states <- append(this_states, list(list(
      date=this_df$time[i],
      a=a_last,
      P=P_last
    )))
  }
  
  # Save states + results separately
  saveRDS(this_states, sprintf("%s/states/univariate_%s.rds", OUTPUT_DIR, SITES[[site]]$human_readable))
  
  df_to_save <- this_df |>
                  mutate(code = site,
                         detrended = exp(sapply(this_states, function(x) x$a['level', ]))) |>
                  select(time, code, no2, detrended)
  write_csv(df_to_save, sprintf("%s/data/results.csv", OUTPUT_DIR))
} 

