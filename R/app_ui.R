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
    fluidPage(
    # List the first level UI elements here 
      navbarPage(
  tags$script(HTML("var header = $('.navbar > .container-fluid');
 header.append('<div style=\"float:right\"><h3>This is R</h3></div>');"
      )),
      tags$script(HTML("var header = $('.navbar > .container-fluid');
header.append('<div style=\"float:right\"><ahref=\"URL\"><img src=\"www/PRSent_logo.png\" alt=\"alt\" style=\"float:right;width:33px;height:41px;padding-top:10px;\"> </a></div>');
    console.log(header)")
      ),
      title = tags$div(img(src="www/PRSent_logo.png", height = '40px', width = '40px'), "something"),
                  fluid = T, 
                  theme = shinythemes::shinytheme("sandstone"),
                  #titlePanel(title=img(src='www/PRSent_logo.png')),
                  tabPanel("Gene-set Analysis Viewer",
                    mod_Gene_set_regression_ui("Gene_set_regression_ui_1")
                  )
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
      app_title = 'PRS.viewing.app'
    ),
    
    tags$link(rel = "shortcut icon", type = "image/png", href = "www/PRSent_logo.png"),
    tags$title("Browser tab title")
    

    #shinyalert::useShinyalert()
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}

