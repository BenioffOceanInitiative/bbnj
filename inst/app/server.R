shinyServer(function(input, output, session) {

  # selected layer ---
  rx_prjres <- reactiveVal("")

  rx_lyr <- reactive({

    #req(input$sel_type)
    req(input$sel_lyr)

    sel_type <- lyr_choices %>%
      filter(label == input$sel_lyr) %>%
      pull(type)
    #browser()
    prjres <- lyr_choices %>%
      filter(label == input$sel_lyr) %>%
      pull(prjres)
    rx_prjres(prjres)

    #browser()
    message(glue(
      "input$sel_lyr: {input$sel_lyr}
       sel_type: {sel_type}
       input$sel_type: {input$sel_type}"))
    if (sel_type == "output_scenario"){
      stopifnot(input$sel_lyr %in% scenarios$label)
      o_lyr <- scenarios %>%
        filter(label == input$sel_lyr) %>%
        pull(gcs_shp) %>%
        read_sf()
    } else {
      #browser()
      stopifnot(input$sel_lyr %in% names(s_features))
      o_lyr <- raster(s_features, input$sel_lyr)
    }
    #s_lyrs <- get(glue("s_layers{prjres}"))
    o_lyr
  })

  # map fxn by projection ----
  get_map_mer <- function(){
    leaflet() %>%
      addProviderTiles(providers$Esri.OceanBasemap, group="ESRI Ocean Basemap") %>%
      addProviderTiles(providers$Stamen.TonerLite, group="Stamen TonerLite")
  }
  get_map_gcs <- function(){

    leaflet(
      #elementId = "map_env",
      options = leafletOptions(
        crs                = leafletCRS(crsClass = "L.CRS.EPSG4326"),
        minZoom            = 1,
        worldCopyJump      = T,
        attributionControl = F)) %>%
      addTiles(
        "//tile.gbif.org/4326/omt/{z}/{x}/{y}@1x.png?style=gbif-geyser",
        group = "Basemap (from GBIF: Geyser)")
  }
  get_map_mol <- function(){

    leaflet(
      options =
        leafletOptions(
          crs=leafletCRS(crsClass="L.Proj.CRS", code=glue('EPSG:{projections_lst$mol$epsg}'),
                         proj4def= projections_lst$mol$proj,
                         resolutions = c(65536, 32768, 16384, 8192, 4096, 2048)))) %>%
      addGraticule(style= list(color= '#999', weight= 0.5, opacity= 1)) %>%
      addGraticule(sphere = TRUE, style= list(color= '#777', weight= 1, opacity= 0.25)) %>%
      addPolygons(
        data=p_countries, group = 'land', weight = 1, color = '#4D4D4D') # gplots::col2hex('gray30'): '#4D4D4D'

  }
  add_map_ctl <- function(map){
    if (is.null(map$x$options$crs$code)){
      map %>%
        addLayersControl(
          baseGroups = c("Stamen TonerLite", "ESRI Ocean Basemap"),
          overlayGroups = c("EEZ","IHO 7 Seas","Feature","Scenario"),
          options = layersControlOptions(collapsed = T))
    } else {
      map %>%
        addLayersControl(
          overlayGroups = c("EEZ","IHO 7 Seas","Feature","Scenario"),
          options = layersControlOptions(collapsed = T))
    }
  }
  add_map_plyctl <- function(map){
    map %>%
      addPolygons(
        data=p_eez_s05, label = ~GeoName, group = "EEZ",
        color = "gray", weight=2) %>%
      addPolygons(
        data=p_ihor_s05, group = "IHO 7 Seas",
        color = "blue", weight=3) %>%
      add_map_ctl() %>%
      #addMouseCoordinates() %>%
      hideGroup("EEZ") %>%
      hideGroup("IHO 7 Seas") %>%
      fitBounds(-150, -60, 150, 60)
  }

  # initial map ----
  output$map <- renderLeaflet({
    # input = list(sel_lyr = lyr_choices[[1]][[1]], slider_opacity=0.7)
    prjres <- rx_prjres()

    #message(glue("output$map: pre P <- projections_tbl"))
    P <- projections_tbl %>%
      filter(prjres == !!prjres)

    message(glue("output$map: post P; pre map; P$prj={P$prj}"))
    # missing "_mer36km"
    #browser()
    map <- switch(P$prj,
      "mer" =  get_map_mer(),
      "gcs" =  get_map_gcs(),
      "mol" =  get_map_mol())

    #message(glue("output$map: post map; pre add_map_lyrsctl"))
    map %>%
      add_map_plyctl()

    # p_sol <- p_r %>%
    #   rename(v = s00a.bio.30pct.gl.mer_sol) %>%
    #   filter(!is.na(v), v > 0) %>%
    #   st_union() %>%
    #   st_collection_extract("POLYGON")
    #
    # get_map_mer() %>%
    #   add_map_plylyrs()
    #
    # get_map_mol() %>%
    #   add_map_plylyrs() %>%
    #   addPolygons(
    #     data = p_sol,
    #     color="green", stroke=F, fillOpacity = 0.8)


  })

  # observe sel_type ----
  observe({
    updateSelectInput(
      session, "sel_lyr", "layer",
      choices = lyr_choices %>%
        filter(type == input$sel_type) %>%
        pull(label))
  })

  # ui_lyr ----
  # output$ui_lyr = renderUI({
  #   #req(input$sel_type)
  #
  #   message("ui_lyr pre")
  #
  #   sel_type <- isolate(
  #     ifelse(is.null(input$sel_type), "output_scenario", input$sel_type))
  #
  #   tags <- tagList(
  #     selectInput(
  #       "sel_type", label = "type",
  #       choices = list(
  #         "input" = list(
  #           "features, general"      = "input_general",
  #           "features by taxa, now"  = "input_taxa_now",
  #           "features by taxa, 2100" = "input_taxa_2100"),
  #         "output" = list(
  #           # TODO: categorize outputs
  #           "scenarios" = sel_type)),
  #       selected = "output_scenario"),
  #
  #     selectInput(
  #       "sel_lyr", label = "layer",
  #       choices = lyr_choices %>%
  #         filter(type == sel_type) %>%
  #         pull(label)))
  #
  #   message("ui_lyr post")
  #
  #   tags
  #   })

  # observe btn_report ----
  observeEvent(input$btn_report, {
    scenario <- input$sel_lyr

    session$sendCustomMessage(
      type = 'modalmessage',
      message = list(
        title = scenario,
        src = glue("scenarios/{scenario}.html")))
  })

  # observe sel_lyr and map ----
  observe({

    #r <- raster(lyrs_prj, input$sel_lyr)
    #tif = "~/github/bbnj/inst/app/www/scenarios/s00a.bio.30pct.gl.mol50km_sol.tif"
    #r <- raster(tif)
    opacity <- input$slider_opacity

    message("observe sel_lyr and map")

    o_lyr <- rx_lyr()

    message("  post o_lyr")

    if ("RasterLayer" %in% class(o_lyr)){
      r <- o_lyr

      rng <- range(values(r), na.rm=T)
      pal <- colorNumeric(
        "Spectral", rng, na.color = "transparent", reverse=T)

      message(glue("pre leafletProxy()"))

      leafletProxy("map") %>%
        #removeImage("lyr_feature") %>%
        clearGroup("Feature") %>%
        #removeShape("lyr") %>%
        addRasterImage(
          r, colors = pal, opacity = opacity, project = F,
          layerId="lyr_feature", group="Feature") %>%
        clearControls() %>%
        addLegend(
          group = "Feature", pal = pal,
          position = "bottomright",
          values = rng,
          title = input$sel_lyr %>%
            str_replace(fixed("."), ", ") %>%
            str_replace("_", "<br>")) %>%
        add_map_ctl()
    } else {
      # presumably a polygon feature class
      p <- o_lyr

      leafletProxy("map") %>%
        clearGroup("Scenario") %>%
        addPolygons(
          data = p,
          color="green", stroke=F, fillOpacity = opacity,
          group="Scenario")
    }

    message(glue("post leafletProxy() opacity:{opacity}"))

  })


  # observe map_click ----
  observeEvent(input$map_click, {

    #browser()
    k <- input$map_click[c("lng","lat")]
    xy_gcs <- st_as_sf(as_tibble(k), coords=c("lng","lat"), crs = 4326, remove=F)
    xy_mer <- st_transform(xy_gcs, 3857)
    am_link <- glue("https://www.aquamaps.org/SpeciesList.php?xlat={k$lat}&xlong={k$lng}")

    o_lyr <- rx_lyr()
    if ("RasterLayer" %in% class(o_lyr)){
      v <- extract(o_lyr, as_Spatial(xy_mer))
    } else {
      v <- st_intersects(xy_gcs, o_lyr, sparse=T) %>% apply(1, any)
    }

    popup <- glue("
      lon: {format(k$lng, digits=5)}<br>
      lat: {format(k$lat, digits=5)}<br>
      value: <strong>{v}</strong><br>
      species: <a href='{am_link}'>AquaMaps</a>")

    leafletProxy("map") %>%
      removePopup("lyr_val") %>%
      addPopups(k$lng, k$lat, popup, layerId = "lyr_val")

  })

})
