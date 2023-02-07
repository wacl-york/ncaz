library(shiny)
library(bslib)
library(shinydashboard)
library(shinyWidgets)
library(shinycssloaders)
library(dplyr)
library(data.table)
library(lubridate)
library(plotly)
source("utils.R")
DISCLAIMER <- "Disclaimer: the estimates shown here are not validated and are still undergoing active research, as such they should not be treated as definitive and should be viewed with caution."

ui <- navbarPage(
  theme=bs_theme(bootswatch="flatly", version=5),
  "UK Clean Air Zone monitors",
  
  tabPanel("Newcastle & Gateshead",
           sidebarLayout(
             sidebarPanel(
               useShinydashboard(),
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
               useShinydashboard(),
               tags$head(
                 tags$link(rel="stylesheet", type="text/css", href="styles.css")
               ),
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

plot_detrended <- function(df) {
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
      hoverinfo = "none"
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
      hoverinfo = "none"
    ) %>%
    layout(
      xaxis = list(range = c(today() - months(3), today()),
                   title = ""),
      yaxis = list(title = "NO2 (ppb)"),
      legend = list(
        orientation = "h",
        xanchor = "center",
        x = 0.5
      )
    )
}

plot_intervention <- function(df) {
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
      hoverinfo = "none"
    ) %>%
    layout(
      xaxis = list(range = c(today() - months(3), today()),
                   title = ""),
      yaxis = list(title = "% change in NO2")
    )
}

generate_tab <- function(df, site) {
  this_df <- df[code == site & !is.na(detrended_abs)]
  p1 <- plot_detrended(this_df)
  p2 <- plot_intervention(this_df)
  
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

server <- function(input, output) {
  df <-
    fread(sprintf("%s/data/results.csv", OUTPUT_DIR_FROM_SHINY))
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
  df[, c("bau", "detrended_abs", "intervention_mean_pct") := .(
    ifelse(time >= intervention_date, no2 / intervention_abs, NA),
    ifelse(time >= intervention_date, detrended_abs * intervention_abs, detrended_abs),
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
  
  
}

shinyApp(ui = ui, server = server)
