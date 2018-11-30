shinyServer(function(input, output, session) {

  get_df <- reactive({

    df <- q %>%
      filter(
        mag   >= input$slider_mag[1],
        mag   <= input$slider_mag[2],
        depth >= input$slider_depth[1],
        depth <= input$slider_depth[2])
    df$var <- df[[input$select_var]]
    df
  })

  output$hist_var <- renderPlotly({

    df     <- get_df()
    n_bins <- min(c(30, length(unique(df$var))))
    lab_x  <- c("mag"="Magnitude (richter)", "depth"="Depth (m)")[input$select_var]

    g <- ggplot2::ggplot(df, aes(x=var)) +
      geom_histogram(bins=n_bins) +
      xlab(lab_x) + ylab("Count")

    p <- plotly::ggplotly(g)
    p$elementId <- NULL # https://github.com/ropensci/plotly/issues/985
    p
  })

  output$map <- renderLeaflet({

    leaflet(abnj) %>%
      addProviderTiles(providers$Stamen.TonerLite) %>%
      addPolygons() %>%
      fitBounds(-150, -60, 150, 60)

  })
})
