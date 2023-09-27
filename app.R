#Libraries
library(shiny)
library(shinydashboard)
library(DT)
library(Biostrings)
library(shinyFiles)
library(reticulate)
library(dplyr)
library(shinyjs)

source("ui.R")
source("server.R")


shinyApp(ui, server)