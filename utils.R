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
    sprintf("https://uk-air.defra.gov.uk/openair/R_data/%s_2023.RData",
            site)
  download.file(base_url, tmp_fn, method = "wget")
  load(tmp_fn)
  # Load hourly data
  df <- get(sprintf("%s_2023", site)) %>%
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
                              intervention_date,
                              alpha = NULL,
                              P = NULL) {
  if (is.null(P)) {
    P <- model$Ptt[, , model$dims$n]
  }
  if (is.null(alpha)) {
    alpha <- t(t(model$att[model$dims$n,]))
  }
  # Z Needs to handle newdata!

  # time-update alpha and P
  # NB: not updating alpha since it should be T*alpha + Bu_k, but T is diagonal and u_k is 0
  # Likewise this should technically be T * P * t(T) + Q but the T is redundant as diagonal
  P <- P + model$model$Q[, , 1]
  
  Z <- array(dim = c(1, 8))
  # Temp, ws, windx, windy, intervention, yearly-sin, yearly-cos, level
  Z[1, 1] <- newdata$air_temp
  Z[1, 2] <- log(newdata$ws)
  Z[1, 3] <- cos(2 * pi * newdata$wd / 360)
  Z[1, 4] <- sin(2 * pi * newdata$wd / 360)
  Z[1, 5] <- ifelse(newdata$time >= intervention_date, 1, 0)
  Z[1, 6] <- sin(2 * pi * (yday(newdata$time) / 365))
  Z[1, 7] <- cos(2 * pi * (yday(newdata$time) / 365))
  Z[1, 8] <- 1
  
  H <- model$model$H[, , 1]
  # Most recent observation
  y <- log(newdata$no2)
  
  # Get: P_t|t-1, alpha_t|t-1, Z_t, R_t
  # Calculate kalman gain as:
  K = P %*% t(Z) %*% (Z %*% P %*% t(Z) + H) ^ -1
  alpha_update = alpha + K %*% (y - Z %*% alpha)
  P_update = (diag(8) - K %*% Z) %*% P
  
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
