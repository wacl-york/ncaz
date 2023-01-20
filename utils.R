INTERVENTION_DATE <- as_datetime("2023-02-01")
OUTPUT_DIR <- "/mnt/shiny/ncaz"
SITES <- list(
  NEWC=list(
    human_readable='central'
  ),
  NCA3=list(
    human_readable='outer'
  )
)

update_univariate <- function(model, newdata, alpha=NULL, P=NULL) {
  if (is.null(P)) {
    P <- model$P[, , model$dims$n+1]
  }
  if (is.null(alpha)) {
    alpha <- t(t(model$a[model$dims$n+1, ]))
  }
  # Z Needs to handle newdata!
  
  Z <- array(dim=c(1, 8))
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
  K = P %*% t(Z) %*% (Z%*%P%*%t(Z) + H)^-1
  alpha_update = alpha + K %*% (y - Z%*%alpha)
  P_update = (diag(8) - K%*%Z) * P * t(diag(8)-K%*%Z) + K%*%H%*%t(K)
  
  list(
    a=alpha_update,
    P=P_update
  )
}

update_multiple_dates <- function(data, model, states) {
  a_new <- states[[length(states)]]$a
  P_new <- states[[length(states)]]$P
  for (i in 1:nrow(data)) {
    updated <- update_univariate(model, data[i, ], a_new, P_new)
    a_new <- updated$a
    P_new <- updated$P
    states <- append(states, list(list(
      date=data$time[i],
      a=a_new,
      P=P_new
    )))
  }
  states
}