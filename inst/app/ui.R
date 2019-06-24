dashboardPage(

  dashboardHeader(title = "BBNJ"),

  dashboardSidebar(

    sidebarMenu(
      menuItem("map", tabName = "tab_map", icon = icon("map")),
      menuItem("histogram", tabName = "tab_hist", icon = icon("bar-chart"))),

    selectInput(
      "sel_type", label = "type",
      choices = list(
        "input" = list(
          "features, general"      = "input_general",
          "features by taxa, now"  = "input_taxa_now",
          "features by taxa, 2100" = "input_taxa_2100"),
        "output" = list(
          # TODO: categorize outputs
          "scenarios" = "output_scenario")),
      selected = "output_scenario"),

    selectInput(
      "sel_lyr", label = "layer",
      choices = lyr_choices$label),
      # choices = lyr_choices %>%
      #   filter(type == "output_scenario") %>%
      #   pull(label)),

    conditionalPanel(
      condition = "input.sel_lyr.substring(0,1) == 's'",
      actionButton("btn_report", "Scenario report")),

    sliderInput(
      "slider_opacity", "opacity", 0, 1, 0.7, step=0.1)),

  dashboardBody(
    tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
        tags$script(src = "modal-handler.js")),
    # tags$style(type = "text/css", "#map {height: calc(100vh - 80px) !important;}"),
    # tags$style(type = "text/css", "#hist_var {height: calc(100vh - 80px) !important;}"),

    bs_modal(
      "modal", "title", size = "large",
      HTML('<iframe data-src="modal.html" height="100%" width="100%" frameborder="0"></iframe>')),

    tabItems(
      tabItem("tab_map",
              leafletOutput("map")),
      tabItem("tab_hist",
              "todo"
              #plotlyOutput("hist_var")
              ))))
