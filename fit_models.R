# Fits Univariate detrending models for 2 Newcastle AURN sites
# Outputs a saved model object along with the latest states
# TODO: Fit multivariate as well?
library(tidyverse)
library(lubridate)
library(openair)
library(KFAS)
library(bssm)
library(forecast)
source("utils.R")

OUTPUT_DIR <- OUTPUT_DIR_FROM_SHINY

# Will imagine that this started on 15th January to allow for a few days
# of manually updating the filter
START_DATE <- as_date("2023-01-15")
#raw_df <- importAURN(site=names(SITES), year = 2010:2023) |>
#            select(time=date, code, no2, ws, wd, air_temp)
#saveRDS(raw_df, sprintf("%s/data/training_2010_2023-01-24.rds", OUTPUT_DIR))
        
raw_df <- readRDS(sprintf("%s/data/training_2010_2023-01-24.rds", OUTPUT_DIR))

# Average to hourly, including vector averaging of wind direction
df_daily <- raw_df |>
  mutate(time = as_datetime(as_date(time)),
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
      intervention=0,
      sin_year = sin(2 * pi * (yday(time) / 365)),
      cos_year = cos(2 * pi * (yday(time) / 365))
)

earliest_clean_dates <- df_daily %>%
                          mutate(is_clean = complete.cases(.)) |>
                          filter(is_clean) |>
                          group_by(code) |>
                          summarise(earliest_clean = min(time)) |>
                          ungroup()
df_daily <- df_daily |>
  inner_join(earliest_clean_dates, by="code") |>
  filter(time >= earliest_clean) |>
  select(-earliest_clean)

# Firstly get priors for annual cycle
calculate_doy_average <- function(df, pollutant = "no2") {
  df$doy <- yday(df$time)
  df$pollutant <- df[[pollutant]]
  df_smooth <- df |>
    group_by(doy) |>
    summarise(pollutant = mean(pollutant, na.rm=T)) |>
    ungroup() |>
    mutate(smooth = filters::sma(pollutant, 30)) 
  
   # Now interpolate over the new year
  rbind(df_smooth |> mutate(year=1),
        df_smooth |> mutate(year=2)) |>
    mutate(rownum = row_number(),
           interp = zoo::na.spline(smooth, na.rm=F)) |>
    filter(between(rownum, 30*2, (365.25*2)-(30*2))) |>
    distinct(doy, interp) |>
    group_by(doy) |>
    summarise(interp = mean(interp, na.rm=T)) |>
    ungroup() |>
    arrange(doy) |>
    pull(interp)
}

calculate_fft_coefs <- function(x, n_coefs=1) {
  # Calculate FFT coefficients of LOG NO2
  fft_yearly_raw <- fft(x)
  fft_yearly <- fft_yearly_raw[2:(length(fft_yearly_raw)/2 + 1)]
  coefs <- array(0, dim=c(n_coefs*2, 1))
  coefs[seq(1, (n_coefs*2)-1, by=2), 1] <- -Im(fft_yearly[1:n_coefs])
  coefs[seq(2, n_coefs*2, by=2), 1] <- Re(fft_yearly[1:n_coefs])
  coefs / length(x)
}

fit_univariate <- function(df, H=NA, Q=NA) {
  n_yearly <- 1
  n_covars <- 5
  n_covars_total <- n_yearly*2 + n_covars
  # NB: removed days of week as it doesn't seem to be doing what I expect, 
  # i.e. it doesn't have any variance in Z-matrix
  form <- as.formula("~ air_temp + log(ws) + wind_x_nows + wind_y_nows + intervention + sin_year + cos_year")
  
  q_diag <- matrix(0, ncol=n_covars_total, nrow=n_covars_total)
  
  # Obtain prior for yearly Fourier
  yearly_smooth <- calculate_doy_average(df)
  coefs_yearly <- calculate_fft_coefs(log(yearly_smooth))

  a1 <- c(rep(0, n_covars),
          coefs_yearly[1:(n_yearly*2), 1])
  
  mod <- SSModel(log(df$no2) ~ 
                 SSMtrend(1, Q=Q) +                                   
                 SSMregression(form,
                               a1=a1,
                               data=df, Q=q_diag),
                  H=H)
  mod <- rename_states(mod, 
                c(
                  "Temperature", "WindSpeed", "WindX", "WindY", "Intervention",
                  paste0(paste0("yearly_", rep(c("sin", "cos"), n_yearly)), rep(1:n_yearly, each=2)),
                  "level"
                )
  )
  
  n_nas <- sum(is.na(mod$Q)) + sum(is.na(mod$H))
  if (n_nas > 0) {
    mod <- fitSSM(mod, inits=log(rep(.01, n_nas)))$model
  }
  mod
}

#---------------------------------------------
# Test manual update code on a single site
#---------------------------------------------
df_central <- df_daily |> filter(code == 'NEWC')
mod_central <- fit_univariate(df_central |> filter(time < START_DATE ))
filt_central <- KFS(mod_central)

# Sanity check on states
df_central |>
  filter(time < START_DATE) |>
  select(time, no2) |>
  mutate(level = exp(filt_central$alphahat[, 'level'])) |>
  filter(time >= '2015-01-01') |>
  mutate(level = level - min(level)) |>
  pivot_longer(-time) |>
  ggplot(aes(x=time, y=value, colour=name)) +
    geom_line(alpha=0.5)

# Now iterate through each day and add new matrices onto the old
df_to_update <- df_central[which(df_central$time >= START_DATE), ]
results <- update_multiple_dates(df_to_update, 
                                 filt_central, 
                                 list(list(date=START_DATE - days(1),
                                           a=NULL,
                                           P=NULL)))

# Compare what happens if would filter on full dataset, i.e. do we get the same results?
filt_all <- KFS(fit_univariate(df_central, Q=mod_central$Q[8, 8, 1], H=mod_central$H[1, 1, 1]))

# That'll do!
days_since_start <- as.numeric(difftime(max(df_daily$time), START_DATE, units="days"))
manual_filtered <- do.call('rbind', lapply(2:length(results), function(i) {
  as.numeric(t(results[[i]]$a))
}))
full_filtered <- filt_all$a[(filt_all$dims$n+1 - days_since_start):(filt_all$dims$n+1), ]

exp(full_filtered)
exp(manual_filtered)
cbind(manual=exp(manual_filtered[, 8]), kfs=exp(full_filtered[, 8]))

#-------------------------------------
# Fit models for all sites
#-------------------------------------
dfs <- list()
for (site in names(SITES)) {
  # Fit model, run the Kalman Filter, and save the final states
  mod <- fit_univariate(df_daily |> filter(code == site, time < START_DATE ))
  filt <- KFS(mod, filtering = "state", smoothing = "none")
  final_states <- list(list(date=START_DATE - days(1),
                              a=t(t(filt$a[filt$dims$n+1, ])),
                              P=filt$P[, , filt$dims$n+1]))
  
  detrended_ind <- which(colnames(filt$att) == 'level')
  
  # Save dataframe containing the detrended time-series
  dfs[[site]] <- df_daily |> 
                  filter(time < START_DATE, code == site) |> 
                  select(time, code) |>
                  mutate(detrended = as.numeric(filt$att[, detrended_ind]),
                         detrended_var = as.numeric(filt$Ptt[detrended_ind, detrended_ind, ]))
  
  # Save models and states
  saveRDS(filt, sprintf("%s/models/univariate_%s.rds", OUTPUT_DIR, site))
  saveRDS(final_states, sprintf("%s/states/univariate_%s.rds", OUTPUT_DIR, site))
}

# Combine detrended time-series from every site with the measured NO2s
df_to_save <- df_daily |> filter(time < START_DATE) |> select(time, code, no2)
df_detrended <- bind_rows(dfs)
df_to_save <- df_to_save |> inner_join(df_detrended, by=c("time", "code"))
write_csv(df_to_save, sprintf("%s/data/results.csv", OUTPUT_DIR))
