library(shiny)
library(listviewer)
library(yaml)

# put some data in environment so it will show up
data(mtcars)

ui <- shinyUI(
  fluidPage(
    jsoneditOutput( "jsed" )
  )
)

server <- function(input,output){
  output$jsed <- renderJsonedit({
    yml <- ("~/Gdrive Ecoquants/projects/bbnj/data/bbnj_config_v01.yml")

    jsonedit(
      #as.list( .GlobalEnv )
      #yaml.load_file(system.file("htmlwidgets/jsonedit.yaml",package="listviewer"))
      yaml.load_file(yml)
      ,"change" = htmlwidgets::JS('function(){
        console.log( event.currentTarget.parentNode.editor.get() )
      }')
    )

  })
}

runApp( list( ui = ui, server = server ) )
