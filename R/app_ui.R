#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    
    # List the first level UI elements here 

      fluidPage(
        
                  theme = shinythemes::shinytheme("united"),
                  #theme = "theme1.css",
                  title = "PRSent your data clearly",
                  #titlePanel(title=img(src='www/PRSent_logo.png')),
                  titlePanel(title = "PRSent your data clearly"),
                    tabPanel("Gene-set Analysis Viewer",
                      mod_Gene_set_regression_ui("Gene_set_regression_ui_1")
                    )
                  
      )
    )
  
}

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' @import data.table
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
 
  
  tags$head(
    favicon(ext = 'png'),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'PRSent'
    ),
    
    tags$link(rel = "shortcut icon", type = "image/png", href = "www/PRSent_logo.png"),
    tags$link(rel = "shortcut icon", type = "image/png", href = "www/GitHub-Mark.png"),
    tags$link(rel="stylesheet", type="text/css", href="www/test.css"),
    
    tags$title("PRSent")
    

    #shinyalert::useShinyalert()
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}

