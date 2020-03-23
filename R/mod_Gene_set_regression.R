#' Gene_set_regression UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_Gene_set_regression_ui <- function(id){
  ns <- NS(id)
  tagList(
    sidebarPanel(
      fileInput("file1", "Choose an input file",
                multiple = F),

      uiOutput(ns("Significance_threshold")
      ),
      
    ),
    mainPanel(plotOutput(ns('PvalPlot')))
  )
}
    
#' Gene_set_regression Server Function
#'
#' @noRd 
mod_Gene_set_regression_server <- function(input, output, session){
  ns <- session$ns
 
}
    
## To be copied in the UI
# mod_Gene_set_regression_ui("Gene_set_regression_ui_1")
    
## To be copied in the server
# callModule(mod_Gene_set_regression_server, "Gene_set_regression_ui_1")
 
