library(shiny)
library(shinydashboard)
library(shinyWidgets)
library(shinycssloaders)
library(tidyverse)
library(lubridate)
library(plotly)
source("utils.R")

ui <- navbarPage(
  "Newcastle and Gateshead CAZ",
  useShinydashboard(),
  tabPanel("Live monitor",
           sidebarLayout(
             sidebarPanel(
               h3("Measuring the impact of the Newcastle and Gateshead Clean Air Zone"),
               h4("Newcastle & Gateshead CAZ"),
               HTML(
                 "A Clean Air Zone (CAZ) is being introduced to Newcastle and Gateshead in an attempt to improve air quality."
               ),
               HTML(
                 "The zone covers a large portion of central Newcastle as well as roads into Gateshead."
               ),
               HTML(
                 "Initially penalties only apply to taxis, buses, coaches, and HGVs from 30th January 2023, with vans and LGVs being affected in July 2023."
               ),
               HTML(
                 "See the official website <a href='https://www.breathe-cleanair.com/'>for further details</a>."
               ),
               br(),
               h4("Detrending NO2"),
               HTML(
                 "This page attempts to quantify the impact of the CAZ on reducing emissions of NO2."
               ),
               HTML(
                 "Unfortunately, this isn't as straightforward as simply measuring NO2 directly, as NO2 concentrations are the result of complex physical relationships between both natural and man-made systems."
               ),
               HTML(
                 "Notably, average NO2 levels vary throughout the year, and are further impacted by local meteorological conditions."
               ),
               HTML(
                 "The <span style='color:#1F77B4'>blue</span> lines show raw daily average NO2 measurements from 2 recording sites that are part of DEFRA's Automatic Urban and Rural Network (AURN):"
               ),
               HTML(
                 "<a href='https://uk-air.defra.gov.uk/networks/site-info?site_id=NEWC'>Newcastle Centre</a> and "
               ),
               HTML(
                 "<a href='https://uk-air.defra.gov.uk/networks/site-info?site_id=NCA3'>Newcastle Cradlewell Roadside</a>, and display considerable short-term variance."
               ),
               
               br(),
               h4("Statistical modelling"),
               HTML(
                 "To quantify the impact of the CAZ on reducing NO2 concentrations, a statistical model has been developed to remove the effects of seasonality and weather on measured NO2."
               ),
               HTML(
                 "The resulting values are shown in the <span style='color:#FF7F04'>orange</span> lines, which are much smoother."
               ),
               HTML(
                 "If the CAZ does result in reduced NO2 emissions, it should be evident here."
               ),
               HTML("See the Methodology tab for further details."),
               br(),
               br(),
               HTML(
                 "<span style='color:#FF0000'>Disclaimer: this is a non-validated piece of research work that should not be referenced in its current state.</span>"
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
  tabPanel(
    "Methodology",
    HTML(
      "This page contains a description of the state-space model (SSM) that has been used for the meteorological and seasonal detrending."
    ),
    HTML(
      "A SSM (sometimes referred to as a Kalman-Filter) is a statistical model of a time-series that defines the observed measurements as"
    ),
    HTML(
      "a linear function of one or more unobserved <emph>states</emph>, along with a fully specified model of the state's dynamics."
    ),
    br(),
    br(),
    withMathJax(
      HTML(
        "Observation equation: $$y_t = \\alpha_t + \\lambda_t + \\beta X_t + \\epsilon_t$$"
      )
    ),
    HTML("State equation: $$\\alpha_t = \\alpha_{t-1} + \\eta_t$$"),
    HTML("Where: "),
    HTML("<ul>"),
    HTML("<li>\\(y_t\\) = measured NO2</li>"),
    HTML("<li>\\(\\alpha_t\\) = underlying NO2 detrended (state)</li>"),
    HTML("<li>\\(\\lambda_t\\) = seasonal factors (state)</li>"),
    HTML("<li>\\(X_t\\) = meteorological measurements</li>"),
    HTML(
      "<li>\\(\\beta\\) = coefficients for meteorological factors (state)</li>"
    ),
    HTML("<li>\\(\\epsilon_t\\) = variance unexplained by the model</li>"),
    HTML(
      "<li>\\(\\eta_t\\) = how much the NO2 detrend can update each time-step (on a random walk)</li>"
    ),
    HTML("</ul>")
  )
)

plot_detrended <- function(df) {
  plot_ly(df, x =  ~ time) %>%
    add_lines(y =  ~ no2,
              name = "Measured",
              color = I("#1F77B4")) %>%
    add_lines(
      y =  ~ detrended_abs,
      name = "Detrended",
      color = I("#FF7F0E")
    ) %>%
    add_ribbons(
      ymin = ~ detrended_lower,
      ymax = ~ detrended_upper,
      name = "Detrended +/- 2sds",
      line = list(color = 'rgba(255, 127, 4, 0)'),
      fillcolor = 'rgba(255, 127, 4, 0.2)',
      showlegend = FALSE,
      hoverinfo = "none"
    ) %>%
    add_lines(y =  ~ bau,
              name = "Business-as-usual",
              color = I("#2CA02C")) %>%
    add_ribbons(
      ymin = ~ bau_lower,
      ymax = ~ bau_upper,
      name = "BAU +/- 2sds",
      line = list(color = 'rgba(44, 160, 44, 0)'),
      fillcolor = 'rgba(44, 160, 44, 0.2)',
      showlegend = FALSE,
      hoverinfo = "none"
    ) %>%
    layout(
      xaxis = list(
        range = c(today() - months(3), today()),
        title = ""
      ),
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
      color = I("#9467BD")
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
      xaxis = list(
        range = c(today() - months(3), today()),
        title = ""
      ),
      yaxis = list(title = "% change in NO2 due to CAZ")
    )
}

generate_tab <- function(df, site) {
  this_df <- df |> filter(code == site, !is.na(detrended_abs))
  p1 <- plot_detrended(this_df)
  p2 <- plot_intervention(this_df)
  
  # Time-series of detrended NO2 + intervention effect
  obj1 <- subplot(p1, p2, nrows = 2, shareX = TRUE, titleY = TRUE) %>%
             layout(hovermode = "x unified")
  
  # Cards displaying the current intervention effect and the date when it went statistically significant
  curr_effect <- this_df %>% slice_max(time, n=1)
  obj2 <- valueBox(sprintf("%.0f%%", curr_effect$intervention_mean_pct),
                   subtitle=sprintf("Intervention effect as of %s", curr_effect$time),
                   icon=icon("bolt-lightning"),
                   color = "green",
                   width=6)
  first_sig <- this_df %>%
                filter(intervention_upper < 1) |>
                slice_min(time, n=1)
  if (nrow(first_sig) == 0) first_sig <- list(time = "Not yet")
  obj3 <- valueBox(sprintf("%s", first_sig$time),
                   subtitle="When intervention was first statistically significant",
                   icon=icon("calendar"),
                   color = "green",
                   width=6)
  
  tagList(fluidRow(obj1), br(), br(), fluidRow(obj2, obj3))
}

server <- function(input, output) {
  
  
  df <-
    read_csv(sprintf("%s/data/results.csv", OUTPUT_DIR_FROM_SHINY))
  df <- df %>%
    filter(time >= as_date("2011-01-01")) %>%
    mutate(
      detrended_abs = exp(detrended + 0.5 * detrended_var),
      detrended_sd = sqrt(detrended_abs ** 2 * (exp(detrended_var) - 1)),
      intervention_abs = exp(intervention + 0.5 * intervention_var),
      intervention_sd = sqrt(intervention_abs ** 2 * (exp(intervention_var) - 1)),
      detrended_abs = detrended_abs - min(detrended_abs, na.rm = T),
      detrended_upper = detrended_abs + 2 * detrended_sd,
      detrended_lower = detrended_abs - 2 * detrended_sd,
      intervention_upper = intervention_abs + 2 * intervention_sd,
      intervention_lower = intervention_abs - 2 * intervention_sd,
      bau = ifelse(time >= INTERVENTION_DATE, no2 * intervention_abs, NA),
      bau_lower = bau * (intervention_abs - 2 * intervention_sd),
      bau_upper = bau * (intervention_abs + 2 * intervention_sd),
      intervention_mean_pct = 100 - (intervention_abs * 100),
      intervention_lower_pct = 100 - ((intervention_abs - 2 * intervention_sd) * 100),
      intervention_upper_pct = 100 - ((intervention_abs + 2 * intervention_sd) * 100)
    )
  
  output$NEWC <- renderUI({
    generate_tab(df, "NEWC")
  })
  output$NCA3 <- renderUI({
    generate_tab(df, "NCA3")
  })
  
  
}

shinyApp(ui = ui, server = server)
