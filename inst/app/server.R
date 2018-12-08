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

    #input = list(select_var = tbl_gfw$name[1])
    tif <- tbl_gfw %>% filter(name == input$select_var) %>% pull(tif)
    r <- raster(tif)

    rng <- range(values(r), na.rm=T)
    pal <- colorNumeric(
      "inferno", rng, na.color = "transparent", reverse=T)

    leaflet() %>%
      addProviderTiles(providers$Esri.OceanBasemap, group="ESRI Ocean Basemap") %>%
      addProviderTiles(providers$Stamen.TonerLite, group="Stamen TonerLite") %>%
      addRasterImage(r, colors = pal, opacity = 0.8, group="Global Fishing Watch") %>%
      addPolygons(data=abnj, group = "High Seas Area") %>%
      addLayersControl(
        baseGroups = c("Stamen TonerLite", "ESRI Ocean Basemap"),
        overlayGroups = c("High Seas Area", "Global Fishing Watch"),
        options = layersControlOptions(collapsed = T)) %>%
      addLegend(
        group = "Global Fishing Watch", pal = pal,
        position = "bottomright",
        values = rng,
        title = input$select_var) %>%
      fitBounds(-150, -60, 150, 60)

  })
})
