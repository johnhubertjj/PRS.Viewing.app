#' Combined_PRS UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_Combined_PRS_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' Combined_PRS Server Function
#'
#' @noRd 
mod_Combined_PRS_server <- function(input, output, session){
  ns <- session$ns
 
}
    
## To be copied in the UI
# mod_Combined_PRS_ui("Combined_PRS_ui_1")
    
## To be copied in the server
# callModule(mod_Combined_PRS_server, "Combined_PRS_ui_1")
 
