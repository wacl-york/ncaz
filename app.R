library(shiny)
library(tidyverse)
library(lubridate)
library(plotly)
source("utils.R")

ui <- fluidPage(

    titlePanel("Modelling the impact of the Newcastle and Gateshead Clean Air Zone"),
    mainPanel(
       plotlyOutput("distPlot")
    )
)

server <- function(input, output) {
  
  df <- read_csv(sprintf("%s/data/results.csv", OUTPUT_DIR_FROM_SHINY))
  df <- df |>
    filter(time >= as_date("2011-01-01")) |>
    mutate(
      intervention_upper = intervention + 2 * sqrt(intervention_var),
      intervention_lower = intervention - 2 * sqrt(intervention_var),
      detrended = detrended - min(detrended),
      bau = ifelse(time >= INTERVENTION_DATE, detrended * intervention, NA),
      bau_lower = no2 * detrended * intervention_lower,
      bau_upper = no2 * detrended * intervention_upper
      )
      

    output$distPlot <- renderPlotly({
      p1 <- plot_ly(df |> filter(code == 'NEWC'), x=~time) |>
            add_lines(y=~no2, 
                      name="Measured",
                      color=I("#1F77B4")
                      ) |>
            add_lines(y=~detrended, 
                      name="Detrended",
                      color=I("#FF7F04")
                      ) |>
            add_lines(y=~bau, 
                      name="Business-as-usual", 
                      color=I("#2CA02C")
                      ) |>
            layout(xaxis=list(range=c(as_datetime("2022-12-01"), max(df$time))),
                   yaxis=list(title="NO2 (ppb)")
            ) |>
            add_annotations(
              text="Newcastle Centre",
              x=0.5,
              y=1,
              yref="paper",
              xref="paper",
              xanchor="middle",
              yanchor="top",
              showarrow=FALSE,
              font=list(size=15)
            )
      
      p2 <- plot_ly(df |> filter(code == 'NCA3'), x=~time) |>
            add_lines(y=~no2, 
                      name="Measured",
                      color=I("#1F77B4"),
                      showlegend=FALSE,
                      ) |>
            add_lines(y=~detrended, 
                      color=I("#FF7F04"),
                      name="Detrended",
                      showlegend=FALSE
                      ) |>
            add_lines(y=~bau, 
                      color=I("#2CA02C"),
                      name="Business-as-usual",
                      showlegend=FALSE
            )|>
            layout(xaxis=list(range=c(as_datetime("2022-12-01"), max(df$time))),
                   yaxis=list(title="NO2 (ppb)")
            ) |>
            add_annotations(
              text="Newcastle Cradlewell Roadside",
              x=0.5,
              y=1,
              yref="paper",
              xref="paper",
              xanchor="middle",
              yanchor="top",
              showarrow=FALSE,
              font=list(size=15)
            )
      
      
       subplot(p1, p2, shareX = TRUE, nrows=2,
                    titleY=TRUE) |>
          layout(hovermode="x unified",
                 xaxis=list(title=""))
    })
}

shinyApp(ui = ui, server = server)
