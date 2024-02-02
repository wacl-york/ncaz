OUTPUT_DIR_FROM_SHINY <- "/mnt/shiny/ncaz"
OUTPUT_DIR_FROM_CRON <- "/shared/storage/shiny0/ncaz"
SITES <- list(
  NEWC = list(human_readable = 'Newcastle Central',
              intervention = as_datetime("2023-01-30"),
              stable_date=as_datetime("2011-01-01")),
  NCA3 = list(human_readable = 'Cradlewell Roadside',
              intervention = as_datetime("2023-01-30"),
              stable_date=as_datetime("2011-01-01")),
  SHDG = list(human_readable = 'Devonshire Green',
              intervention = as_datetime("2023-02-27"),
              stable_date=as_datetime("2016-01-01"))
)

load_data <- function(site) {
  tmp_fn <- "tmp.RData"
  base_url <-
    sprintf("https://uk-air.defra.gov.uk/openair/R_data/%s_2024.RData",
            site)
  download.file(base_url, tmp_fn, method = "wget")
  load(tmp_fn)
  # Load hourly data
  df <- get(sprintf("%s_2024", site)) %>%
    as_tibble() %>%
    mutate(date = as_datetime(date)) %>%
    select(site, code, date, no2 = NO2, air_temp = temp, ws, wd)
  file.remove(tmp_fn)
  
  # Average to daily
  df %>%
    mutate(
      time = as_datetime(as_date(date)),
      wind_y = ws * sin(2 * pi * wd / 360),
      wind_x = ws * cos(2 * pi * wd / 360)
    ) %>%
    group_by(code, time) %>%
    summarise(
      no2 = mean(no2, na.rm = T),
      air_temp = mean(air_temp, na.rm = T),
      wind_y = mean(wind_y, na.rm = T),
      wind_x = mean(wind_x, na.rm = T),
      ws = mean(ws, na.rm = T)
    ) %>%
    ungroup() %>%
    # Convert vector wind back to polar
    mutate(
      wd = as.vector(atan2(wind_y, wind_x) * 360 / 2 / pi),
      wd = ifelse(wd < 0, wd + 360, wd),
      ws_vec = (wind_x ^ 2 + wind_y ^ 2) ^ 0.5,
      wind_y_nows = sin(2 * pi * wd / 360),
      wind_x_nows = cos(2 * pi * wd / 360),
      intervention = as.numeric(time >= SITES[[site]][['intervention']])
    )
}

update_univariate <- function(model,
                              newdata,
                              intervention_start=NULL,
                              intervention_end=NULL,
                              alpha = NULL,
                              P = NULL,
                              transform_func=identity) {
  if (is.null(P)) {
    P <- model$Ptt[, , model$dims$n]
  }
  if (is.null(alpha)) {
    alpha <- t(t(model$att[model$dims$n,]))
  }
  
  m <- model$dims$m
  intervention_offset <- as.numeric(!is.null(intervention_start))
  
  # Z Needs to handle newdata!
  # time-update alpha and P
  # NB: not updating alpha since it should be T*alpha + Bu_k, but T is diagonal and u_k is 0
  # Likewise this should technically be T * P * t(T) + Q but the T is redundant as diagonal
  Z <- array(dim = c(1, m))
  # Temp, ws, windx, windy, intervention, yearly-sin, yearly-cos, level
  Z[1, 1] <- newdata$air_temp
  Z[1, 2] <- log(newdata$ws)
  Z[1, 3] <- cos(2 * pi * newdata$wd / 360)
  Z[1, 4] <- sin(2 * pi * newdata$wd / 360)
  if (!is.null(intervention_start)) {
    if (!is.null(intervention_end)) {
      Z[1, 5] <- ifelse(newdata$time >= intervention_start & newdata$time <= intervention_end, 1, 0)
    } else {
      Z[1, 5] <- ifelse(newdata$time >= intervention_start, 1, 0)
    }
  }
  Z[1, 5 + intervention_offset] <- sin(2 * pi * (yday(newdata$time) / 365))
  Z[1, 6 + intervention_offset] <- cos(2 * pi * (yday(newdata$time) / 365))
  Z[1, 7 + intervention_offset] <- 1
  
  if (any(is.na(Z)) || is.na(newdata$no2)) {
    return(
      list(a = alpha, P=P)
    )
  }
  
  P <- P + model$model$Q[, , 1]
  H <- model$model$H[, , 1]
  # Most recent observation
  y <- transform_func(newdata$no2)
  
  # Get: P_t|t-1, alpha_t|t-1, Z_t, R_t
  # Calculate kalman gain as:
  K = P %*% t(Z) %*% (Z %*% P %*% t(Z) + H) ^ -1
  alpha_update = alpha + K %*% (y - Z %*% alpha)
  P_update = (diag(m) - K %*% Z) %*% P
  
  list(a = alpha_update,
       P = P_update)
}

update_multiple_dates <- function(data, model, intervention_date, states) {
  a_new <- states[[length(states)]]$a
  P_new <- states[[length(states)]]$P
  for (i in 1:nrow(data)) {
    updated <- update_univariate(model, data[i,], intervention_date, a_new, P_new)
    a_new <- updated$a
    P_new <- updated$P
    states <- append(states, list(list(
      date = data$time[i],
      a = a_new,
      P = P_new
    )))
  }
  states
}

# Function to calculate the Daily NO2 average to use as a prior for the annual effect
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
    mutate(rownum = row_number(doy),
           interp = zoo::na.spline(smooth, na.rm=F)) |>
    filter(between(rownum, 30*2, (365.25*2)-(30*2))) |>
    distinct(doy, interp) |>
    group_by(doy) |>
    summarise(interp = mean(interp, na.rm=T)) |>
    ungroup() |>
    arrange(doy) |>
    pull(interp)
}

calculate_hod_average <- function(df, pollutant = "no2") {
  df$hod <- hour(df$time)
  df$pollutant <- df[[pollutant]]
  df_smooth <- df |>
    group_by(hod) |>
    summarise(pollutant = mean(pollutant, na.rm=T)) |>
    ungroup()
  
  df_smooth$pollutant
}

# Function to calculate the FFT coefficients of a trend to use
# as a prior
calculate_fft_coefs <- function(x, n_coefs=1) {
  # Calculate FFT coefficients of LOG NO2
  fft_yearly_raw <- fft(x)
  fft_yearly <- fft_yearly_raw[2:(length(fft_yearly_raw)/2 + 1)]
  coefs <- array(0, dim=c(n_coefs*2, 1))
  coefs[seq(1, (n_coefs*2)-1, by=2), 1] <- -Im(fft_yearly[1:n_coefs])
  coefs[seq(2, n_coefs*2, by=2), 1] <- Re(fft_yearly[1:n_coefs])
  coefs / length(x)
}

# Fit the univariate model
fit_univariate <- function(df, 
                           y_col,
                           yearly_prior_col = NULL,
                           slope=FALSE,
                           H=NULL, 
                           a1_level=NULL,
                           P1_level=NULL,
                           P1inf_level=NULL,
                           Q_level=NULL, 
                           a1_regression=NULL,
                           P1_regression=NULL,
                           P1inf_regression=NULL,
                           Q_regression=NULL,
                           transform_func=identity,
                           intervention_date=NULL,
                           intervention_end=NULL,
                           remove_intercept=TRUE,
                           covariate_formula="~ air_temp + log(ws) + wind_x_nows + wind_y_nows",
                           n_day_of_year=1,
                           n_hour_of_day=0
                           ) {
  trend_degree <- 1 + as.integer(slope)
  n_covars <- str_count(covariate_formula, "\\+") + 1
  if (!is.null(intervention_date)) {
    n_covars <- n_covars + 1
  }
  if (!remove_intercept) {
    n_covars <- n_covars + 1
  }
  n_covars_total <- n_day_of_year*2 + n_covars + n_hour_of_day * 2
  # NB: removed days of week as it doesn't seem to be doing what I expect, 
  # i.e. it doesn't have any variance in Z-matrix
  if (!is.null(intervention_date)) {
    covariate_formula <- paste(covariate_formula, "intervention", sep=" + ")
  }
  if (n_day_of_year > 0) {
    form_yearly <- paste(rep(c("sin_year", "cos_year"), n_day_of_year), 
                         rep(seq(n_day_of_year), each=2),
                         collapse="+",
                         sep="")
    covariate_formula <- paste(covariate_formula, form_yearly, sep='+')
  }
  if (n_hour_of_day > 0) {
    form_hour <- paste(rep(c("sin_day", "cos_day"), n_hour_of_day), 
                       rep(seq(n_hour_of_day), each=2),
                       collapse="+",
                       sep="")
    covariate_formula <- paste(covariate_formula, form_hour, sep='+')
  }
  form <- as.formula(covariate_formula)
  
  if (is.null(H)) {
    H <- NA
  }
  
  # Level
  if (is.null(Q_level)) {
    Q_level <- lapply(1:trend_degree, function(x) NA)
  }
  if (is.null(a1_level)) {
    a1_level <- rep(0, trend_degree)
  }
  if (is.null(P1_level)) {
    P1_level <- matrix(0, ncol=trend_degree, nrow=trend_degree)
  }
  if (is.null(P1inf_level)) {
    P1inf_level <- diag(1, trend_degree)
  }
  
  # Regression
  if (is.null(Q_regression)) {
    Q_regression <- matrix(0, ncol=n_covars_total, nrow=n_covars_total)
  }
  if (is.null(a1_regression)) {
    a1_regression <- c(rep(0, n_covars))
    
    if (n_day_of_year > 0) {
      yearly_smooth <- calculate_doy_average(df, yearly_prior_col)
      coefs_yearly <- calculate_fft_coefs(transform_func(yearly_smooth),
                                          n_coefs = n_day_of_year)
      a1_regression <- c(a1_regression, coefs_yearly[, 1])
    }
    
    if (n_hour_of_day > 0) {
      hourly_smooth <- calculate_hod_average(df, yearly_prior_col)
      coefs_hourly <- calculate_fft_coefs(transform_func(hourly_smooth), 
                                          n_coefs = n_hour_of_day)
      a1_regression <- c(a1_regression, coefs_hourly[, 1])
    }
  }
  if (is.null(P1_regression)) {
    P1_regression <- matrix(0, nrow=n_covars_total, ncol=n_covars_total)
  }
  if (is.null(P1inf_regression)) {
    P1inf_regression <- diag(n_covars_total)
  }
  
  if (!is.null(intervention_date)) {
    if (!is.null(intervention_end)) {
      df$intervention <- as.numeric(df$time >= intervention_date & df$time <= intervention_end)
    } else {
      df$intervention <- as.numeric(df$time >= intervention_date)
    }
  }
  
  mod <- SSModel(transform_func(df[[y_col]]) ~ 
                 SSMtrend(trend_degree, 
                          a1 = a1_level,
                          P1 = P1_level,
                          P1inf = P1inf_level,
                          Q=Q_level) +
                 SSMregression(form,
                               a1=a1_regression,
                               P1=P1_regression,
                               P1inf = P1inf_regression,
                               remove.intercept = remove_intercept,
                               data=df, Q=Q_regression),
                  H=H)
  #state_names <- c(
  #  "Temperature", "WindSpeed", "WindX", "WindY"
  #)
  #
  #if (!is.null(intervention_date)) {
  #  state_names <- c(state_names, "intervention")
  #}
  #if (!remove_intercept) {
  #  state_names <- c("intercept", state_names)
  #}
  #if (n_day_of_year > 0) {
  #  state_names <- c(state_names, 
  #                   paste0(paste0("yearly_", rep(c("sin", "cos"), n_day_of_year)), 
  #                          rep(1:n_day_of_year, each=2)))
  #}
  #if (n_hour_of_day > 0) {
  #  state_names <- c(state_names, 
  #                   paste0(paste0("hourly_", rep(c("sin", "cos"), n_hour_of_day)), 
  #                          rep(1:n_hour_of_day, each=2))
  #                  )
  #  
  #}
  #state_names <- c(state_names, 'level')
  #if (slope) {
  #  state_names <- c(state_names, 'slope')
  #}
  #mod <- rename_states(mod, state_names)
  
  n_nas <- sum(is.na(mod$Q)) + sum(is.na(mod$H))
  if (n_nas > 0) {
    if (n_nas == 1) {
      mod <- fitSSM(mod, inits=log(rep(.01, n_nas)), method="Brent",
                    lower=transform_func(0.00001), upper=transform_func(20))$model
    } else {
      mod <- fitSSM(mod, inits=log(rep(.01, n_nas)))$model
    }
  }
  mod
}

manual_filter <- function(model,
                          df,
                          first_date,
                          intervention_start=NULL,
                          intervention_end=NULL,
                          transform_func=identity,
                          all_states=FALSE) {
  rows_to_update <- which(df$time >= as_datetime(first_date))
  
  if (!is.null(intervention_start)) {
    intervention_start <- as_datetime(intervention_start)
  }
  
  states <- list(list())
  a_new <- t(t(model$att[model$dims$n, ]))
  P_new <- model$Ptt[, , model$dims$n]
  for (i in rows_to_update) {
    # The artificial bump I've used to kickstart the intervention
    if (!is.null(intervention_start) && df$time[i] == intervention_start) {
      P_new[5, 5] <- 1e5
    }
    updated <- update_univariate(model,
                                 df[i,], 
                                 intervention_start=intervention_start, 
                                 intervention_end=intervention_end,
                                 alpha=a_new, 
                                 P=P_new,
                                 transform_func=transform_func)
    a_new <- updated$a
    P_new <- updated$P
    states <- append(states, list(list(
      date = df$time[i],
      a = a_new,
      P = P_new
    )))
  }
  
  level_ind <- which(colnames(model$att) == 'level')
  if (is.null(intervention_start)) {
    df_1 <- data.table(
        detrended=as.numeric(model$att[, level_ind]),
        detrended_var=as.numeric(model$Ptt[level_ind, level_ind, ])
      )
    df_2 <- data.table(
        detrended = vapply(states[2:length(states)],  
                           function(x) x$a[level_ind, 1], numeric(1)),
        detrended_var = vapply(states[2:length(states)],  
                           function(x) x$P[level_ind, level_ind], numeric(1))
      )
  } else {
    int_ind <- which(colnames(model$att) == 'intervention')
    df_1 <- data.table(
        detrended=as.numeric(model$att[, level_ind]),
        detrended_var=as.numeric(model$Ptt[level_ind, level_ind, ]),
        intervention=as.numeric(model$att[, int_ind]),
        intervention_var=as.numeric(model$Ptt[int_ind, int_ind, ])
      )
    df_1[, adjusted := detrended + intervention ]
    df_2 <- data.table(
        detrended = vapply(states[2:length(states)],  
                           function(x) x$a[level_ind, 1], numeric(1)),
        detrended_var = vapply(states[2:length(states)],  
                           function(x) x$P[level_ind, level_ind], numeric(1)),
        intervention = vapply(states[2:length(states)], 
                              function(x) x$a[int_ind, 1], numeric(1)),
        intervention_var = vapply(states[2:length(states)],  
                           function(x) x$P[int_ind, int_ind], numeric(1))
      )
    df_2[, adjusted := detrended + intervention ]
  }
  
  states_int_detrended <- rbindlist(
      list(df_1, df_2)
  )
  
  # Either return all states or just the intervention + detrended
  if (all_states) {
    map_dfr(states[2:length(states)], function(x) {
      as_tibble(t(x$a[, 1]))
    })
  } else {
    states_int_detrended 
  }
}
