shinyServer(function(input, output, session) {

  output$map <- renderLeaflet({
    # input = list(sel_lyr = lyr_choices[[1]][[1]], slider_opacity=0.7)

    data(p_abnj_s05)
    data(p_ppow_s05)

    leaflet() %>%
      addProviderTiles(providers$Esri.OceanBasemap, group="ESRI Ocean Basemap") %>%
      addProviderTiles(providers$Stamen.TonerLite, group="Stamen TonerLite") %>%
      addPolygons(
        data=p_eez_s05, label = ~GeoName, group = "EEZ",
        color = "gray", weight=2) %>%
      # addPolygons(
      #   data=p_abnj_s05, group = "High Seas",
      #   color = "black", weight=4) %>%
      # addPolygons(
      #   data=p_ppow_s05, group = "Pelagic Provinces",
      #   color = "blue", weight=3) %>%
      addPolygons(
        data=p_ihor_s05, group = "IHO 7 Seas",
        color = "blue", weight=3) %>%
      addMouseCoordinates() %>%
      addLayersControl(
        baseGroups = c("Stamen TonerLite", "ESRI Ocean Basemap"),
        #overlayGroups = c("EEZ","High Seas","Pelagic Provinces", "IHO 7 Seas"),
        overlayGroups = c("EEZ","IHO 7 Seas"),
        options = layersControlOptions(collapsed = T)) %>%
      hideGroup("EEZ") %>%
      hideGroup("High Seas") %>%
      hideGroup("Pelagic Provinces") %>%
      hideGroup("IHO 7 Seas") %>%
      fitBounds(-150, -60, 150, 60)
  })

  observeEvent(input$btn_report, {
    #browser()
    # scenario <- input$sel_lyr %>%
    #   str_replace("\\.", "_") %>%
    #   str_replace_all("\\.", "-")
    scenario <- input$sel_lyr
    url <- glue("https://ecoquants.com/bbnj/articles/{scenario}.html")

    showModal(modalDialog(
      size = "l",
      title = input$sel_lyr,
      tagList(
        a(href=url, scenario),
        p(),
        tags$iframe(
          src=url,
          height="100%", width="100%", margin="0", `min-height`="700px", frameborder="0"))
    ))
  })

  # observe layer selection
  observe({

    #browser()

    r <- raster(lyrs_mer, input$sel_lyr)
    opacity <- input$slider_opacity

    rng <- range(values(r), na.rm=T)
    pal <- colorNumeric(
      "Spectral", rng, na.color = "transparent", reverse=T)

    leafletProxy("map") %>%
      removeImage("lyr") %>%
      addRasterImage(
        r, colors = pal, opacity = opacity, project = F,
        layerId="lyr", group="Layer") %>%
      clearControls() %>%
      addLegend(
        group = "Layer", pal = pal,
        position = "bottomright",
        values = rng,
        title = input$sel_lyr %>%
          str_replace(fixed("."), ", ") %>%
          str_replace("_", "<br>")) %>%
      addLayersControl(
        baseGroups = c("Stamen TonerLite", "ESRI Ocean Basemap"),
        #overlayGroups = c("EEZ","High Seas","Pelagic Provinces", "Layer"),
        overlayGroups = c("EEZ","IHO 7 Seas"),
        options = layersControlOptions(collapsed = T))
  })


  ## observe map clicks
  observeEvent(input$map_click, {

    #browser()
    k <- input$map_click[c("lng","lat")]

    r <- raster(lyrs_mer, input$sel_lyr)

    xy <- sp::SpatialPoints(
      as.matrix(as.numeric(k)) %>% t(),
      proj4string = CRS("+init=epsg:4326")) %>%
      spTransform(CRS("+init=epsg:3857"))

    cellid <- cellFromXY(r, xy)
    v <- r[cellid]

    am_link <- glue("https://www.aquamaps.org/SpeciesList.php?xlat={k$lat}&xlong={k$lng}")

    leafletProxy("map") %>%
      removePopup("lyr_val")

    if (!is.na(v)){
      popup <- glue("
                    lon: {format(k$lng, digits=5)}<br>
                    lat: {format(k$lat, digits=5)}<br>
                    value: <strong>{v}</strong><br>
                    species: <a href='{am_link}'>AquaMaps</a>")

      leafletProxy("map") %>%
        addPopups(k$lng, k$lat, popup, layerId = "lyr_val")
    }


  })

})
