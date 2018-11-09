library(shiny)
library(tidyverse)
library(leaflet)
library(glue)
library(shinydashboard)
library(plotly)

# load your data
q <- quakes

add <- function(a, b){
  a + b
}