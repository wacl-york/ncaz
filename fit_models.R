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
meta <- importMeta() 
newcastle_sites <- meta |>
          filter(grepl("Newcastle", site)) |>
          pull(code)
        
df <- readRDS(sprintf("%s/data/training_2010_2023-01-19.rds", OUTPUT_DIR))

# Average to hourly, including vector averaging of wind direction
df_daily <- df |>
  mutate(time = as_datetime(as_date(time)),
         wind_y = ws * sin(2 * pi * wd / 360),
         wind_x = ws * cos(2 * pi * wd / 360)) |>
  group_by(site, code, time) |>
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
      intervention=0
)

# Prep data for modelling with yearly terms
yearly_coefs <- as_tibble(
                          fourier(ts(df_daily |> filter(code == 'NEWC') |> pull(time), frequency=365), K=1)) |>
                mutate(time = df_daily |> filter(code == 'NEWC') |> pull(time))
df_daily <- df_daily |>
  inner_join(yearly_coefs, by="time")

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
  # TODO remove days of week!
  # doesn't seem to be doing what I expect, i.e. it doesn't have any variance
  # in Z-matrix
  form_covar <- "~ air_temp + log(ws) + wind_x_nows + wind_y_nows + intervention"
  form_yearly <- paste0(rep(c("S", "C"), n_yearly), rep(1:n_yearly, each=2), "-365")
  form <- as.formula(sprintf("%s + `%s`", 
                  form_covar,
                  paste0(form_yearly, collapse="` + `")
                  ))
  
  q_diag <- matrix(0, ncol=n_covars_total, nrow=n_covars_total)
  
  # Obtain prior for yearly Fourier
  yearly_smooth <- calculate_doy_average(df)
  coefs_yearly <- calculate_fft_coefs(log(yearly_smooth))

  a1 <- c(rep(0, n_covars),
          coefs_yearly[1:(n_yearly*2), 1])
  
  mod <- SSModel(log(df$no2) ~ 
                 SSMtrend(1, Q=Q) +                                   
                 #SSMseasonal(7, Q=0, sea.type="dummy") +
                 SSMregression(form,
                               a1=a1,
                               data=df, Q=q_diag),
                  H=H)
  mod <- rename_states(mod, 
                c(
                  "Temperature", "WindSpeed", "WindX", "WindY", "Intervention",
                  paste0(paste0("yearly_", rep(c("sin", "cos"), n_yearly)), rep(1:n_yearly, each=2)),
                  "level"
                  #'Fri',
                  #'Sat',
                  #'Sun',
                  #'Mon',
                  #'Tues',
                  #'Weds'
                )
  )
  
  n_nas <- sum(is.na(mod$Q)) + sum(is.na(mod$H))
  if (n_nas > 0) {
    mod <- fitSSM(mod, inits=log(rep(.01, n_nas)))$model
  }
  mod
}

# ------------- FIT MODELS FOR EACH SITE!
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

# ------------- Update states
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
filt_all$a[(filt_all$dims$n - 2):(filt_all$dims$n+1), ]
rbind(
  as.numeric(t(results[[2]]$a)),
  as.numeric(t(results[[3]]$a)),
  as.numeric(t(results[[4]]$a)),
  as.numeric(t(results[[5]]$a))
)

# So fit second site, generate states up until 15th January, then manually update states using update_models.R
mod_outer <- fit_univariate(df_daily |> filter(code == 'NCA3', time < START_DATE ))
filt_outer <- KFS(mod_outer)

saveRDS(filt_central, sprintf("%s/models/univariate_central.rds", OUTPUT_DIR))
saveRDS(filt_outer, sprintf("%s/models/univariate_outer.rds", OUTPUT_DIR))

# Save initial states which will pull from model later
initial_states_central <- list(list(date=START_DATE - days(1),
                                    a=t(t(filt_central$a[filt_central$dims$n+1, ])),
                                    P=filt_central$P[, , filt_central$dims$n+1]))
initial_states_outer <- list(list(date=START_DATE - days(1),
                                    a=t(t(filt_outer$a[filt_outer$dims$n+1, ])),
                                    P=filt_outer$P[, , filt_outer$dims$n+1]))
    
saveRDS(initial_states_central, sprintf("%s/states/univariate_central.rds", OUTPUT_DIR))
saveRDS(initial_states_outer, sprintf("%s/states/univariate_outer.rds", OUTPUT_DIR))

# Save data to DB
df_to_save <- df_daily |> filter(time < START_DATE) |> select(time, code, no2)
# Now add new states in
df_detrended_central <- df_daily |> 
                filter(time < START_DATE) |>
                select(time, code) |>
                filter(code == 'NEWC') |>
                mutate(detrended = exp(c(filt_central$att[, 'level'])))
df_detrended_outer <- df_daily |> 
                filter(time < START_DATE) |>
                select(time, code) |>
                filter(code == 'NCA3') |>
                mutate(detrended = exp(c(filt_outer$att[, 'level'])))
df_to_save <- df_to_save |>
  inner_join(df_detrended_central |> rbind(df_detrended_outer), by=c("time", "code")) |>
  mutate(intervention=0, intervention_var=0)
write_csv(df_to_save, sprintf("%s/data/results.csv", OUTPUT_DIR))
