library(shiny)
library(bslib)
library(shinydashboard)
library(shinyWidgets)
library(shinycssloaders)
library(shinyjs)
library(dplyr)
library(data.table)
library(lubridate)
library(plotly)
library(openair)
library(tidytable)
library(KFAS)
library(digest)
source("utils.R")
DISCLAIMER <- "Disclaimer: the estimates shown here are not validated and are still undergoing active research, as such they should not be treated as definitive and should be viewed with caution."

DATA_DIR <- sprintf("%s/data", OUTPUT_DIR_FROM_SHINY)
SITES_META <- readRDS(sprintf("%s/aurn_sites_meta.rds", DATA_DIR))
ALL_SITES <- SITES_META$code
names(ALL_SITES) <- SITES_META$site

ratio_to_pct_difference <- function(x) {
  (x * 100) - 100
}

TRANSFORMS <- list(
  "Absolute"= list(
    func=identity,
    inverse_func=identity,
    units="ppb",
    subtract_func=`-`,
    add_func=`+`,
    mean_func=function(mean, var) mean,
    sd_func=function(mean, var) sqrt(var),
    int_post_func=identity,
    sig_threshold=0
  ),
  "Relative"= list(
    func=log,
    inverse_func=exp,
    units="%",
    subtract_func=`/`,
    add_func=`*`,
    mean_func=function(mean, var) exp(mean + 0.5 * var),
    sd_func=function(mean, var) sqrt(mean**2 * (exp(var)-1)),
    int_post_func=ratio_to_pct_difference,
    sig_threshold=1
  )
)

Q_to_sd <- function(q, transform, transform_type) {
  q_sd <- transform(sqrt(q))
  if (transform_type == 'Relative') {
    q_sd <- (q_sd * 100) - 100
  }
  q_sd
}

sd_to_Q <- function(q_sd, transform, transform_type) {
  if (transform_type == 'Relative') {
    q_sd <- q_sd/100 + 1
  }
  (transform(q_sd))**2
}

dataset_hash <- function(location, transform, intervention, date) {
  str <- if (intervention == 'Yes') {
    paste(location, transform, intervention, date[1], date[2], sep='_')
  } else {
    paste(location, transform, intervention, sep='_')
  }
  digest(str)
}

model_hash <- function(location, transform, intervention, date, Q) {
  str <- if (intervention == 'Yes') {
    paste(location, transform, intervention, date[1], date[2], Q, sep='_')
  } else {
    paste(location, transform, intervention, Q, sep='_')
  }
  digest(str)
}

preprocess_aurn <- function(df) {
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
      sin_year = sin(2 * pi * (yday(time) / 365)),
      cos_year = cos(2 * pi * (yday(time) / 365))
    )
}

load_dataset <- function(site) {
  # Firstly load the training data up until mid January
  fn <- sprintf("%s/sites_training_2015_2023-01-26/%s.rds", DATA_DIR, site)
  df_1 <- readRDS(fn) |> select(-intervention)
  
  # Then load AURN data to fill in the gaps up until today
  years_to_load <- seq(year(max(df_1$time)), year(today()), by=1)
  df_2 <- importAURN(site = site, year=years_to_load) |> as.data.table()
  
  # Process AURN data
  df_2 <- preprocess_aurn(df_2) |>
            filter(time > max(df_1$time))
  
  # Combine
  rbindlist(list(df_1, df_2))
}

estimateQ <- function(df, transform, intervention, date) {
  # HOWEVER, if also have intervention then don't want to use intervention in training!
  # Maybe should split data into training and test
  # So if have intervention then training date is up until the day before
  # If have no intervention then training date is end of 2019
  int_start <- if (intervention == 'Yes') date[1] else NULL
  int_end <- if (intervention == 'Yes') date[2] else NULL
  train_end <- if (intervention == 'Yes') date[1] - days(1) else as_datetime('2019-12-31')
  mod <- fit_univariate(df |> filter(time <= train_end),
                        'no2',
                        yearly_prior_col = 'no2',
                        transform_func=TRANSFORMS[[transform]][['func']],
                        intervention_date=int_start,
                        intervention_end=int_end
                        )
  level_col <- 7 + as.numeric(intervention == 'Yes')
  mod$Q[level_col, level_col, 1]
}
  


ui <- navbarPage(
  theme=bs_theme(bootswatch="flatly", version=5),
  "UK Clean Air Zone monitors",
  
  tabPanel("Newcastle & Gateshead",
           sidebarLayout(
             sidebarPanel(
               useShinydashboard(),
               useShinyjs(),
               tags$head(
                 tags$link(rel="stylesheet", type="text/css", href="styles.css")
               ),
               h4("Newcastle & Gateshead CAZ"),
               HTML(
                 "The Newcastle and Gateshead Clean Air Zone covers a large portion of central Newcastle as well as roads into Gateshead."
               ),
               HTML(
                 "Initially penalties only apply to taxis, buses, coaches, and HGVs from 30th January 2023, with vans and LGVs being affected in July 2023."
               ),
               HTML(
                 "See the official website <a href='https://www.breathe-cleanair.com/'>for further details</a>."
               ),
               br(),
               br(),
               HTML(
                 "In the top panel, the <span style='color:#1F77B4'>blue</span> lines show raw daily average NO2 measurements from 2 recording sites that are part of DEFRA's Automatic Urban and Rural Network (AURN):"
               ),
               HTML(
                 "<a href='https://uk-air.defra.gov.uk/networks/site-info?site_id=NEWC'>Newcastle Centre</a> and "
               ),
               HTML(
                 "<a href='https://uk-air.defra.gov.uk/networks/site-info?site_id=NCA3'>Newcastle Cradlewell Roadside</a>, and display considerable short-term variance."
               ),
               HTML(
                 "The <span style='color:#FF7F04'>orange</span> lines display the detrended NO2, where we would expect to see a decrease
                 if the CAZ is having its desired effect.
                 The <span style='color:#9467BD'>purple</span> lines in the bottom panel displays an estimate of the CAZ's relative effect, with a summary in the bottom 2 boxes."
               ),
               HTML("See the Methodology tab for further details."),
               br(),
               br(),
               HTML(
                 sprintf("<span class='disclaimer'>%s</span>", DISCLAIMER)
               )
             ),
             mainPanel(tabsetPanel(
               type = "tabs",
               tabPanel("Newcastle Central",
                        withSpinner(uiOutput("NEWC"))),
               tabPanel("Cradlewell Roadside",
                        withSpinner(uiOutput("NCA3")))
             ))
           )),
  tabPanel("Sheffield",
           sidebarLayout(
             sidebarPanel(
               h4("Sheffield CAZ"),
               HTML(
                 "The Sheffield CAZ covers the inner ring road and the city centre."
               ),
               HTML(
                 "Initially penalties will only apply to the most polluting HGVs, LGVs, vans, buses, coaches, and taxis from 27th February 2023."
               ),
               HTML(
                 "See the official website <a href='https://www.sheffield.gov.uk/campaigns/clean-air-zone-sheffield'>for further details</a>."
               ),
               br(),
               br(),
               HTML(
                 "In the top panel, the <span style='color:#1F77B4'>blue</span> line shows raw daily average NO2 measurements from the <a href='https://uk-air.defra.gov.uk/networks/site-info?site_id=SHDG'>Devonshire Green</a>
                 monitoring site that is part of DEFRA's Automatic Urban and Rural Network (AURN)."
               ),
               HTML(
                 "The <span style='color:#FF7F04'>orange</span> line display the detrended NO2, where we would expect to see a decrease
                 if the CAZ is having its desired effect.
                 The <span style='color:#9467BD'>purple</span> line in the bottom panel displays an estimate of the CAZ's relative effect, with a summary in the bottom 2 boxes."
               ),
               HTML("See the Methodology tab for further details."),
               br(),
               br(),
               HTML(
                 sprintf("<span class='disclaimer'>%s</span>", DISCLAIMER)
               )
             ),
             mainPanel(tabsetPanel(
               type = "tabs",
               tabPanel("Devonshire Green",
                        withSpinner(uiOutput("SHDG")))
             ))
           )),
  tabPanel(
    "Fit a custom model",
    h3("Fit a custom model"),
    sidebarLayout(
       sidebarPanel(
         p("This page allows you to experiment with the modelling parameters to identify
           their effect on the resulting detrended series. You can generate a detrended model
           from one of the 49 AURN sites that had 75% NO2 availability since 2015."),
         selectInput("custom_site", "AURN location",
                        choices=c("Choose a site"="", ALL_SITES)),
         tags$div(title=paste("Absolute means that the intervention reduction is always the same,",
                              "but a relative effect is proortional to the NO2 background.",
                              "I.e. an absolute effect of 5 ppb is the same if the background is",
                              "10ppb or 50ppb, but a 5% relative effect results in a 0.5ppb",
                              "drop in a 10ppb background, but 2.5ppb when the background is 50ppb"),
           radioButtons("custom_transform", "Model type (hover for details)",
                        c("Absolute", "Relative"), inline=TRUE,
                        selected=character(0))
         ),
         radioButtons("custom_intervention", "Model an intervention?",
                      c("Yes", "No"), inline=TRUE,
                      selected=character(0)),
         shinyjs::hidden(dateRangeInput("custom_date", "Intervention range",
                        start=today() - weeks(1), end=today())),
         fluidRow(
           column(
             width=6,
             uiOutput("custom_q_ui")
           ),
           column(
             width=6,
             shinyjs::disabled(actionButton("estimate_q", "Estimate automatically"))
           )
         ) |> tagAppendAttributes(class = "custom_q_row"),
         shinyjs::disabled(actionButton("fit", "Fit model"))
       ),
       mainPanel(
         withSpinner(uiOutput("custom"))
       )
    )
  ),
  tabPanel(
    "Methodology",
   h3("Motivation"),
   HTML(
     "To see if a CAZ has resulted in a reduction in harmful air pollution, we can't just look at the measured pollution 
     before and after the CAZ was introduced, as NO2 concentrations are the result of complex physical relationships 
     between both natural and man-made systems and as such vary considerably day-to-day."
   ),
   HTML(
     "Notably, average NO2 levels vary throughout the year (peaking in Winter), and are further impacted by 
      local meteorological conditions."
   ),
   HTML(
     "As such, a model is used to <em>detrend</em> the raw measurements so that any changes are due to external factors,
     such as policy interventions."
   ),
   HTML("In particular, we require a statistical model that can handle time-series data with external covariates, can 
        track underlying changes, and provides a full probabilistic formulation of every parameter."),
   h3("State-space modelling"),
  HTML(
      "A State-space model (SSM), which is sometimes referred to as a Kalman-Filter, is a statistical model of a 
      time-series that defines the observed measurements as"
    ),
    HTML(
      "a linear function of one or more unobserved <emph>states</emph> (observation equation), along with a fully specified model of the 
      state's dynamics (state equation).
      The general forms of the SSM equations are as follows."
    ),
    br(),
    br(),
    withMathJax(
      HTML(
        "Observation equation: $$y_t = \\alpha_t + \\beta X_t + \\epsilon_t$$"
      )
    ),
    HTML("State equation: $$\\alpha_t = \\alpha_{t-1} + \\eta_t$$"),
    HTML("Where: "),
    HTML("<ul>"),
    HTML("<li>\\(y_t\\) = observed data of the outcome</li>"),
    HTML("<li>\\(\\alpha_t\\) = the underlying trend: part 1 of the state</li>"),
    HTML(
      "<li>\\(\\beta\\) = regression coefficients for the external covariates: part 2 of the state</li>"
    ),
    HTML("<li>\\(X_t\\) = observed data of the external covariates</li>"),
    HTML("<li>\\(\\epsilon_t\\) = error unexplained by the model, assumed to follow a constant variance \\(\\Sigma_\\epsilon\\)</li>"),
    HTML(
      "<li>\\(\\eta_t\\) = how much the trend can update each time-step on a random walk, assumed to follow a constant variance \\(\\Sigma_\\eta\\). 
      NB: \\(\\beta\\) can be included here so that the regression coefficients are time-varying, but this is not used in this model</li>"
    ),
    HTML("</ul>"),
    h3("Detrending NO2"),
    HTML("In particular for the NO2 detrending, these parameters represent:"),
    HTML("<ul>"),
    HTML("<li>\\(y_t\\) = measured NO2</li>"),
    HTML("<li>\\(\\alpha_t\\) = underlying detrended NO2</li>"),
    HTML("<li>\\(X_t\\) = meteorological measurements and temporal variables</li>"),
    HTML(
      "<li>\\(\\beta\\) = coefficients for meteorological and temporal factors</li>"
    ),
    HTML("<li>\\(\\epsilon_t\\) = error unexplained by the model</li>"),
    HTML(
      "<li>\\(\\eta_t\\) = how much the NO2 detrended series can update each day (on a random walk)</li>"
    ),
    HTML("</ul>"),
    HTML("\\(\\Sigma_\\epsilon\\) and \\(\\Sigma_\\eta\\) can be manually specified if they are known from a mechanistic
          knowledge of the dynamics, but in this case they are automatically estimated using Maximum Likelihood Estimation.
          "),
    br(),
    HTML("Furthermore, a log transform is applied to the outcome in order to help stabilize the variance (NO2 exhibits 
         right skew), and to enforce positivity."),
    h3("Estimating the CAZ's impact"),
    HTML("To quantify the impact of the CAZ on reducing NO2 concentrations in a more formal manner than simply observing a decrease in the detrended series,
         an intervention variable is added to the \\(X_t\\) covariate matrix, equal to a 1 when the CAZ is in effect, and a 0 beforehand.
         The resulting associated coefficient in \\(\\beta\\) quantifies the effect of the NO2, which is in relative/percentage terms since a log transform
         is used on the outcome.
         It is this value that is displayed as the <span style='color:#9467BD'>purple</span> time-series, along with its associated 95% confidence interval."),
    h3("Daily update"),
    HTML("Every day at 6am GMT, the model is updated with the previous day's average concentration, providing a new estimate of both the detrended series and the intervention effect.
    A SSM can be updated in two ways: <em>filter</em> or <em>smoother</em>. The filter only uses the most recent data available for each estimate, i.e. the detrend estimate
    for the 15th January only uses data that was available on the 15th January, even if it's now 30th January. A smoother by contrast uses all data, so on the 30th January
    rather than just updating the detrend estimate for the 29th, the smoother would create state estimates for every day. The smoother produces, well, smoother results,
    and is less susceptible to quick changes in trend. However it feels slightly misleading to portray it as an <em>online</em> method unlike the filter. In high resolution
    time-series the filter is preferred for online situations since it is far less computationally demanding (since it just updates one timepoint rather than all); since
    this is just daily data the smoother could be used here as there wouldn't be any time-limitations, but it feels more in keeping with the stated purpose of a 'real-time'
    dashboard to use the filter.
    "),
    br(),
    br()
  ),
  tabPanel(
    "About",
    h3("Modelling the effect of UK Clean Air Zones on reducing NO2 concentrations"),
    HTML("<p>This website contains live dashboards to monitor the impact of Clean Air Zones (CAZ) on NO2 concentrations in multiple cities in the UK.
      It has been developed by researchers at the <a href='https://www.york.ac.uk/chemistry/research/wacl/'>Wolfson Atmospheric Chemistry Laboratories</a> 
      at the University of York as part of a research project into techniques for identifying local changes in NO2 emissions arising from policy changes.
      The dashboard was setup prior to the introduction of CAZs in two UK cities: Newcastle (30th January 2023) and Sheffield (27th February 2023),
      providing a real-time online estimate of the CAZ's effectiveness.
      This provides a more realistic estimate of how much information can be gleamed in real-time, rather than a post-hoc study once a significant duration has passed with the benefit of hindsight.
      However, it also means there is additional uncertainty in the estimates, both due to the limited information contained in the data, but also due teething issues being identified in the methodology itself.
      As such, the following disclaimer is provided on each dashboard as a reminder that this is a work-in-progress and the estimates should not be used as-is without thoroughly understanding the limitations of the 
      approach.
      </p>"),
    HTML(sprintf("<span class='disclaimer'>%s</span>", DISCLAIMER)),
    HTML("<p>
    A modelling based approach is used to extract changes in the underlying NO2 concentrations, free of confounding factors such as the local meteorology and seasonal factors.
       The resulting detrended series is used to identify the changes since the intervention took place. See the Methodology tab for full details of the modelling.
         </p>"),
    h3("Contact"),
    HTML("Please send me an email <code>stuart.lacy</code> at <code>york.ac.uk</code> if you have any questions or would like to discuss this work."),
    br(),
    HTML("The source code for both the state-space modelling and the Shiny web-app are publicly available on GitHub: 
    <a href='https://github.com/wacl-york/ncaz'>https://github.com/wacl-york/ncaz</a>
    ")
  ),
)

plot_detrended <- function(df, x_start, x_end) {
  plot_ly(df, x =  ~ time) %>%
    add_lines(y =  ~ no2,
              name = "Measured",
              color = I("#1F77B4")) %>%
    add_lines(
      y =  ~ detrended_abs,
      name = "Detrended",
      legendgroup='Detrended',
      color = I("#FF7F0E")
    ) %>%
    add_ribbons(
      ymin = ~ detrended_lower,
      ymax = ~ detrended_upper,
      name = "Detrended +/- 2sds",
      legendgroup='Detrended',
      line = list(color = 'rgba(255, 127, 4, 0)'),
      fillcolor = 'rgba(255, 127, 4, 0.2)',
      showlegend = FALSE,
      hoverinfo = "none",
      connectgaps=TRUE
    ) %>%
    add_lines(y =  ~ bau,
              name = "Business-as-usual",
              legendgroup='BAU',
              color = I("#2CA02C")) %>%
    add_ribbons(
      ymin = ~ bau_lower,
      ymax = ~ bau_upper,
      name = "BAU +/- 2sds",
      legendgroup='BAU',
      line = list(color = 'rgba(44, 160, 44, 0)'),
      fillcolor = 'rgba(44, 160, 44, 0.2)',
      showlegend = FALSE,
      hoverinfo = "none",
      connectgaps=TRUE
    ) %>%
    config(
      displaylogo=FALSE
    ) %>%
    layout(
      xaxis = list(range = c(x_start, x_end),
                   title = ""),
      yaxis = list(title = "NO2 (ppb)"),
      legend = list(
        orientation = "h",
        xanchor = "center",
        x = 0.5
      )
    )
}

plot_intervention <- function(df, x_start, x_end) {
  plot_ly(df, x =  ~ time) %>%
    add_lines(
      y =  ~ intervention_mean_pct,
      name = "Intervention effect",
      color = I("#9467BD"),
      showlegend = FALSE
    ) %>%
    add_ribbons(
      ymin = ~ intervention_lower_pct,
      ymax = ~ intervention_upper_pct,
      name = "Intervention effect +/- 2sds",
      line = list(color = 'rgba(148, 103, 189, 0)'),
      fillcolor = 'rgba(148, 103, 189, 0.2)',
      showlegend = FALSE,
      hoverinfo = "none",
      connectgaps=TRUE
    ) %>%
    config(
      displaylogo=FALSE
    ) %>%
    layout(
      xaxis = list(range = c(x_start, x_end),
                   title = ""),
      yaxis = list(title = "% change in NO2")
    )
}

generate_tab <- function(df, site) {
  this_df <- df[code == site & !is.na(detrended_abs)]
  p1 <- plot_detrended(this_df, x_start=today() - months(1), today())
  p2 <- plot_intervention(this_df, x_start=today() - months(1), today())
  
  # Time-series of detrended NO2 + intervention effect
  obj1 <-
    subplot(p1,
            p2,
            nrows = 2,
            shareX = TRUE,
            titleY = TRUE) %>%
    layout(hovermode = "x unified")
  
  # Cards displaying the current intervention effect and the date when it went statistically significant
  curr_effect <- this_df %>% slice_max(time, n = 1)
  obj2 <-
    valueBox(
      sprintf("%.0f%%", curr_effect$intervention_mean_pct),
      subtitle = sprintf("Intervention effect as of %s", curr_effect$time),
      icon = icon("bolt-lightning"),
      color = "green",
      width = 6
    )
  first_sig <- this_df %>%
    filter(intervention_upper < 1) %>%
    slice_min(time, n = 1)
  if (nrow(first_sig) == 0)
    first_sig <- list(time = "Not yet")
  obj3 <- valueBox(
    sprintf("%s", first_sig$time),
    subtitle = "When intervention was first statistically significant",
    icon = icon("calendar"),
    color = "green",
    width = 6
  )
  
  tagList(fluidRow(obj1), br(), br(), fluidRow(obj2, obj3))
}

generate_custom_tab <- function(df, transform, int_start, int_end) {
  if (is.null(int_start)) {
    df <- df |> mutate(intervention = 0, intervention_var = 0)
    x_start <- today() - months(1)
    x_end <- today()
    int_start <- min(df$time)
    int_end <- max(df$time)
  } else {
    x_start <- int_start
    x_end <- int_end
  }
  
  # Preprocess!
  df <- df |>
          mutate(
            # Convert from raw mean and var to mean and SD
            detrended = transform$mean_func(detrended, detrended_var),
            intervention = transform$mean_func(intervention, intervention_var),
            detrended_sd = transform$sd_func(detrended, detrended_var),
            intervention_sd = transform$sd_func(intervention, intervention_var),
            # Calculate BAU and adjust detrended for intervention for plot
            bau=ifelse(time >= int_start & time <= int_end, transform$subtract_func(no2, intervention), NA),
            # CIs
            detrended_lower = detrended- 2 * detrended_sd,
            intervention_lower = intervention- 2 * intervention_sd,
            bau_lower = transform$add_func(bau, intervention_lower),
            detrended_upper = detrended + 2 * detrended_sd,
            intervention_upper = intervention + 2 * intervention_sd,
            bau_upper = transform$add_func(bau, intervention_upper),
            # Get intervention into human-readable format
            intervention_mean_pct = transform$int_post_func(intervention),
            intervention_lower_pct = transform$int_post_func(intervention_lower),
            intervention_upper_pct = transform$int_post_func(intervention_upper)
          ) |>
          rename(detrended_abs = detrended)  # Rename to be compatible with legacy code
  
  # Rescale detrending
  #min_once_settled <- df[ time >= (min(time) + years(1)), min(detrended, na.rm=T)]
  #df[, detrended := detrended - min_once_settled ]
  
  # Currently have time, detrended, intervention, adjusted
  p1 <- plot_detrended(df |> filter(time > (min(time) + years(1))), x_start, x_end)
  p2 <- plot_intervention(df |> filter(time > (min(time) + years(1))), x_start, x_end)
  
  # Time-series of detrended NO2 + intervention effect
  obj1 <-
    subplot(p1,
            p2,
            nrows = 2,
            shareX = TRUE,
            titleY = TRUE) %>%
    layout(hovermode = "x unified")
  
  # Cards displaying the current intervention effect and the date when it went statistically significant
  curr_effect <- df %>% slice_max(time, n = 1)
  obj2 <-
    valueBox(
      sprintf("%.0f%s", curr_effect$intervention_mean_pct, transform$units),
      subtitle = sprintf("Intervention effect as of %s", curr_effect$time),
      icon = icon("bolt-lightning"),
      color = "green",
      width = 6
    )
  first_sig <- df %>%
    filter(intervention_upper < transform$sig_threshold) %>%
    slice_min(time, n = 1)
  if (nrow(first_sig) == 0)
    first_sig <- list(time = "Not yet")
  obj3 <- valueBox(
    sprintf("%s", first_sig$time),
    subtitle = "When intervention was first statistically significant",
    icon = icon("calendar"),
    color = "green",
    width = 6
  )
  
  tagList(fluidRow(obj1), br(), br(), fluidRow(obj2, obj3))
}

server <- function(input, output) {
  df <-
    fread(sprintf("%s/results.csv", DATA_DIR))
  sites_dt <- rbindlist(lapply(SITES, as.data.table), idcol = "code")[, .(code, stable_date, intervention_date=intervention)]
  df <- sites_dt[df, on=.(code)]
  df <- df[ time >= stable_date]
  # When have 0 measurements, don't plot detrended or intervention. The whole point of SSM is that
  # we still have these latent estimates even when missing observations, but it's more intuitive
  # for a general audience to not plot
  df[ is.na(no2), `:=` (detrended=NA, detrended_var=NA, intervention=NA, intervention_var=NA)]
  
  # Convert log(mean) and log(var) into mean and sd
  df[, c("detrended_abs", "intervention_abs") := .(exp(detrended + 0.5 * detrended_var),
                                                   exp(intervention + 0.5 * intervention_var))]
  df[, c("detrended_sd", "intervention_sd") := .(sqrt(detrended_abs ** 2 * (exp(detrended_var) - 1)),
                                                 sqrt(intervention_abs ** 2 * (exp(intervention_var) - 1)))]
  
  # Create Business as Usual and relative intervention (%) columns
  df[, c("bau", "intervention_mean_pct") := .(
    ifelse(time >= intervention_date, no2 / intervention_abs, NA),
    (intervention_abs * 100) - 100
  )]

  # Rescale detrending
  df[, detrended_abs := detrended_abs - min(detrended_abs, na.rm = T)]
  
  # Create CIs
  df[, c("detrended_lower",
         "intervention_lower",
         "bau_lower",
         "intervention_lower_pct") := .(
           detrended_abs - 2 * detrended_sd,
           intervention_abs - 2 * intervention_sd,
           bau * (intervention_abs - 2 * intervention_sd),
           (intervention_abs - 2 * intervention_sd) * 100 - 100
         )]
  df[, c("detrended_upper",
         "intervention_upper",
         "bau_upper",
         "intervention_upper_pct") := .(
           detrended_abs + 2 * detrended_sd,
           intervention_abs + 2 * intervention_sd,
           bau * (intervention_abs + 2 * intervention_sd),
           (intervention_abs + 2 * intervention_sd) * 100 - 100
         )]
  
  output$NEWC <- renderUI({
    generate_tab(df, "NEWC")
  })
  output$NCA3 <- renderUI({
    generate_tab(df, "NCA3")
  })
  output$SHDG <- renderUI({
    generate_tab(df, "SHDG")
  })
  
  ######## CUSTOM MODEL
  observeEvent(input$custom_intervention, {
    if (input$custom_intervention == 'Yes') {
      shinyjs::show("custom_date")
    } else {
      shinyjs::hide("custom_date")
    }
  }, ignoreInit = TRUE)
  
  output$custom_q_ui <- renderUI({
    numericInput("custom_q", label="Day-to-day NO2 standard deviation", value = NA, min = 0, max=100, step=0.5)
  })
  
  observeEvent(input$custom_transform, {
    lab <- "Day-to-day NO2 standard deviation"
    updateNumericInput(inputId="custom_q",
                       label=sprintf("%s (%s)", lab, TRANSFORMS[[input$custom_transform]][['units']]),
                       value=character(0))
  }, ignoreNULL = TRUE)
  
  # Only enable the automatic identification of Q when have filled in all the values
  observeEvent({
    input$custom_site
    input$custom_transform
    input$custom_intervention
    input$custom_date
  }, {
    if (input$custom_site != '' && !is.null(input$custom_transform) && !is.null(input$custom_intervention)) {
     shinyjs::enable("estimate_q")
    }
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  # Only enable Fit model button when all fields have been populated
  observeEvent({
    input$custom_site
    input$custom_transform
    input$custom_intervention
    input$custom_date
    input$custom_q
  }, {
    if (input$custom_site != '' && !is.null(input$custom_transform) && !is.null(input$custom_intervention) && !is.na(input$custom_q)) {
     shinyjs::enable("fit")
    }
    
  }, ignoreNULL = TRUE, ignoreInit = TRUE)
  
  estimated_Qs <- list()
  datasets <- list()
  models <- list()
  current_filter <- reactiveVal()
  
  observeEvent(input$estimate_q, {
    hash <- dataset_hash(input$custom_site, input$custom_transform, input$custom_intervention, input$custom_date)
    if (!hash %in% names(estimated_Qs)) {
      if (! input$custom_site %in% names(datasets)) {
        withProgress(message = 'Downloading data from AURN...', value= 0.33, {
          datasets[[input$custom_site]] <<- load_dataset(input$custom_site)
        })
      }
      # Estimate Q, but NB: Q is used as a variance in KFAS, and also if we have a log transform we want to get back to pcts
      withProgress(message = 'Estimating variance...', value= 0.66, {
        raw_q <- estimateQ(datasets[[input$custom_site]], input$custom_transform, input$custom_intervention, input$custom_date)
      })
      estimated_Qs[[hash]] <<- Q_to_sd(raw_q, TRANSFORMS[[input$custom_transform]][['inverse_func']], input$custom_transform)
    }
      
    q <- estimated_Qs[[hash]]
    updateNumericInput(inputId="custom_q",
                       value=q)
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  observeEvent(input$fit, {
    hash <- model_hash(input$custom_site, input$custom_transform, input$custom_intervention, input$custom_date, input$custom_q)
    if (!hash %in% names(models)) {
      
      if (! input$custom_site %in% names(datasets)) {
        datasets[[input$custom_site]] <<- load_dataset(input$custom_site)
      }
      df <- datasets[[input$custom_site]]
      
      # Estimate Q, but NB: Q is used as a variance in KFAS, and also if we have a log transform we want to get back to pcts
      int_start <- if (input$custom_intervention == 'Yes') input$custom_date[1] else NULL
      int_end <- if (input$custom_intervention == 'Yes') input$custom_date[2] else NULL
      train_end <- if (input$custom_intervention == 'Yes') input$custom_date[1] - days(1) else as_datetime('2019-12-31')
      # Fit H and model. If have already estimated Q then in theory can already know H but ah well
      mod <- fit_univariate(df |> filter(time <= train_end),
                            'no2',
                            yearly_prior_col = 'no2',
                            transform_func=TRANSFORMS[[input$custom_transform]][['func']],
                            intervention_date=int_start,
                            intervention_end=int_end,
                            Q_level=sd_to_Q(input$custom_q, TRANSFORMS[[input$custom_transform]][['func']], input$custom_transform)
                            )
      # Now filter
      filt_1 <- KFS(mod, filtering="state")
      
      filt_2 <- manual_filter(filt_1,
                    df,
                    first_date=train_end + days(1),
                    intervention_start=int_start,
                    intervention_end = int_end,
                    transform_func = TRANSFORMS[[input$custom_transform]][['func']]
                    ) |>
      mutate(time=df$time, no2=df$no2)
      models[[hash]] <<- filt_2
    }
    
    # Update reactive value with latest filter
    current_filter(models[[hash]])
      
  }, ignoreInit = TRUE, ignoreNULL = TRUE)
  
  output$custom <- renderUI({
    if (!is.null(current_filter())) {
      intervention <- isolate(input$custom_intervention)
      date <- isolate(input$custom_date)
      int_start <- if (intervention == 'Yes') date[1] else NULL
      int_end <- if (intervention == 'Yes') date[2] else NULL
      generate_custom_tab(current_filter(), 
                          TRANSFORMS[[isolate(input$custom_transform)]],
                          int_start,
                          int_end)
    }
  })
  
}

shinyApp(ui = ui, server = server)
