#' Decile UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList 
mod_Decile_ui <- function(id){
  ns <- NS(id)
  tagList(
 
  )
}
    
#' Decile Server Function
#'
#' @noRd 
mod_Decile_server <- function(input, output, session){
  ns <- session$ns
 
}
    
## To be copied in the UI
# mod_Decile_ui("Decile_ui_1")
    
## To be copied in the server
# callModule(mod_Decile_server, "Decile_ui_1")
 
