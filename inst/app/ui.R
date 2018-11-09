dashboardPage(

  dashboardHeader(title = "BBNJ"),

  dashboardSidebar(

    sidebarMenu(
      menuItem("Map", tabName = "tab_map", icon = icon("map")),
      menuItem("Histogram", tabName = "tab_hist", icon = icon("bar-chart"))),

    sliderInput(
      "slider_mag", label = "Magnitude",
      min = min(quakes$mag), max = max(quakes$mag), step = 0.2,
      value = c(min(quakes$mag), max = max(quakes$mag))),

    sliderInput(
      "slider_depth", label = "Depth",
      min = min(quakes$depth), max = max(quakes$depth), step = 50,
      value = c(min(quakes$depth), max = max(quakes$depth))),

    selectInput(
      "select_var", label = "Variable to plot",
      choices = c(
        "Magnitude" = "mag",
        "Depth" = "depth"))),

  dashboardBody(
    tags$style(type = "text/css", "#map {height: calc(100vh - 80px) !important;}"),
    tags$style(type = "text/css", "#hist_var {height: calc(100vh - 80px) !important;}"),
    tabItems(
      tabItem("tab_map",
              leafletOutput("map")),
      tabItem("tab_hist",
              plotlyOutput("hist_var")))))
