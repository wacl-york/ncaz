INTERVENTION_DATE <- as_datetime("2023-01-30")
OUTPUT_DIR_FROM_SHINY <- "/mnt/shiny/ncaz"
OUTPUT_DIR_FROM_CRON <- "/shared/storage/shiny0/ncaz"
SITES <- list(
  NEWC = list(human_readable = 'central',
              intervention = as_datetime("2023-01-30"),
              stable_date=as_datetime("2011-01-01")),
  NCA3 = list(human_readable = 'outer',
              intervention = as_datetime("2023-01-30"),
              stable_date=as_datetime("2011-01-01")),
  SHDG = list(human_readable = 'sheffield',
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
                              alpha = NULL,
                              P = NULL) {
  if (is.null(P)) {
    P <- model$P[, , model$dims$n + 1]
  }
  if (is.null(alpha)) {
    alpha <- t(t(model$a[model$dims$n + 1,]))
  }
  # Z Needs to handle newdata!
  
  Z <- array(dim = c(1, 8))
  # Temp, ws, windx, windy, intervention, yearly-sin, yearly-cos, level
  Z[1, 1] <- newdata$air_temp
  Z[1, 2] <- log(newdata$ws)
  Z[1, 3] <- cos(2 * pi * newdata$wd / 360)
  Z[1, 4] <- sin(2 * pi * newdata$wd / 360)
  Z[1, 5] <- ifelse(newdata$time >= as_date("2023-02-01"), 1, 0)
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
  P_update = (diag(8) - K %*% Z) * P * t(diag(8) - K %*% Z) + K %*% H %*%
    t(K)
  
  list(a = alpha_update,
       P = P_update)
}

update_multiple_dates <- function(data, model, states) {
  a_new <- states[[length(states)]]$a
  P_new <- states[[length(states)]]$P
  for (i in 1:nrow(data)) {
    updated <- update_univariate(model, data[i,], a_new, P_new)
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
