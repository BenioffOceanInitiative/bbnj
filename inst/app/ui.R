dashboardPage(

  dashboardHeader(title = "BBNJ"),

  dashboardSidebar(

    sidebarMenu(
      menuItem("Map", tabName = "tab_map", icon = icon("map")),
      menuItem("Histogram", tabName = "tab_hist", icon = icon("bar-chart"))),

    selectInput(
      "sel_lyr", label = "Layer",
      choices = lyr_choices),

    sliderInput(
      "slider_opacity", "Opacity", 0, 1, 0.7, step=0.1)),

  dashboardBody(
    tags$style(type = "text/css", "#map {height: calc(100vh - 80px) !important;}"),
    tags$style(type = "text/css", "#hist_var {height: calc(100vh - 80px) !important;}"),
    tabItems(
      tabItem("tab_map",
              leafletOutput("map")),
      tabItem("tab_hist",
              plotlyOutput("hist_var")))))
