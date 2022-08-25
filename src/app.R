library(shiny)
library(shinydashboard)

source('ui/ui.R', local=TRUE)
source('server/server.R')

shinyApp(ui = ui, server = server)
