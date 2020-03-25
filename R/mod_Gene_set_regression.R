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
      fileInput(ns("file1"), "Choose an input file",
                multiple = F),

      uiOutput(ns("Significance_threshold")
      ),
      checkboxGroupInput(ns("Gene_regions"), label = "Length of Gene regions:",
                           choices = c("Genome-wide","Gene-set")),
      uiOutput(ns("DSM")),
      uiOutput(ns("geneset")),
      
    ),
    mainPanel(plotOutput(ns('PvalPlot')))
  )
}
    
#' Gene_set_regression Server Function
#'
#' @noRd 
mod_Gene_set_regression_server <- function(input, output, session){
  
  ns <- session$ns
  
  My_data <- reactive({ req(input$file1)
    cat(file=stderr(), "accessing_file input system")
    
    # insert alterations functions that cleans up string names
    # source("alterations_of_gene_sets_script.R")
    
    
    ## Read in data
    Full_data <- data.table::fread(input$file1$datapath)
    
    data.table::setnames(Full_data, old = c("Set","Threshold"), new = c("Genesets", "Significance_thresholds"))
    data.table::Full_data[, Gene_regions := Genesets]
    
    #First runthrough of PRS positions
    Positions_of_PRS_in_table <- Calculate_positions_of_genome_wide_PRS(Full_data,"Genesets")
    
    #Genome_wide_PRS <- which(Full_data$Gene_regions == "Base")
    #Gene_set_PRS <- which(Full_data$Gene_regions != "Base")
    
    data.table::Full_data[,Gene_regions := "NA"]
    data.table::Full_data[Positions_of_PRS_in_table$Genome_wide_PRS, Gene_regions := "Genome-wide"]
    data.table::Full_data[Positions_of_PRS_in_table$Gene_set_PRS, Gene_regions := "Gene-set"]
    
    
    ## Create arguments to shiny app
    Gene.sets.input <- unique(Full_data$Genesets)
    significance_threshold.input <- unique(Full_data$Significance_thresholds)
    DSM.input <- "Everything"
    
    
    data.table::Full_data[Positions_of_PRS_in_table$Genome_wide_PRS, Type := "Whole_genome"]
    data.table::Full_data[!Positions_of_PRS_in_table$Genome_wide_PRS, Type:= "Pathway"]
    
    
    # Set all P values to work within the app and add annotations
    if(any(Full_data$P == 0) == T){
      data.table::Full_data[,P_altered := P]
      Full_data[P_altered == 0, P_altered := 1e-300]
      Full_data$logp <- -log10(Full_data$P_altered)
    }else{
      Full_data$logp <- -log10(Full_data$P)
    }
    
    Full_data[,estimate := scale(Coefficient)]
    Full_data[,SE := scale(Standard.Error)]
    
    Full_data$SE_higher <- Full_data$estimate + Full_data$SE
    Full_data$SE_lower <- Full_data$estimate - Full_data$SE
    Full_data$r2_dir <- 100 * (as.numeric(Full_data$R2) *
                                 (sign(as.numeric(Full_data$estimate))))
    Full_data$p_value_text <- paste("p =", scientific(Full_data$P, digits = 2), sep = " ")
    
    if(any(Full_data$p_value_text == "p = 0.0e+00") == TRUE){
      Full_data[p_value_text == "p = 0.0e+00", p_value_text := "p < 1e-300"]
    }
    
    
    # Run a duplication check on the gene-set PRS names
    Full_data <- Duplication_of_gene_sets_check(Data_table = Full_data, Genome_wide_positions = Positions_of_PRS_in_table$Genome_wide_PRS, Significance_thresholds_name = "Significance_thresholds", gene_set_values = Full_data$Genesets)
    
    # Re-calculate the positions of whole_genome_PRS
    #Second runthrough of PRS positions
    Positions_of_PRS_in_table <- Calculate_positions_of_genome_wide_PRS(Full_data,"Genesets")
    
    
    # Change names in alterations 
    alterations <- Pathway_cleanup(Full_data$Genesets, Positions_of_PRS_in_table$Genome_wide_PRS)   
    Full_data[, alterations := alterations]
    
    #Add DSM argument script
    Full_data[, samples.i. := "Everything"]
    
    Full_data[Positions_of_PRS_in_table$Gene_set_PRS, score := paste0(Gene_regions,"_SCORE_",Genesets,"_",Significance_thresholds)]
    Full_data[Positions_of_PRS_in_table$Genome_wide_PRS, score := paste0("All.genome_SCORE_whole_genome_",Significance_thresholds)]
    Full_data[,.id := "PRS Results"]
    
    Full_data <- list(Full_data = Full_data ,Gene.sets.input = Gene.sets.input, significance_threshold.input = significance_threshold.input, DSM.input = DSM.input)
    Full_data
    
  })
  
  
  output$Significance_threshold <- renderUI({ 
    significance_threshold.input <- as.numeric(My_data()$significance_threshold.input)
    checkboxGroupInput(ns("Significance_threshold"), label = "PRS P Value Threshold:",
                       choices = significance_threshold.input[order(significance_threshold.input)], selected = significance_threshold.input)
  })
  
  output$DSM <- renderUI({ 
    DSM.input <- My_data()$DSM.input
    selectInput(ns("DSM"), label = "DSM type:",
                choices = DSM.input)
  })
  
  output$geneset <- renderUI({ 
    Gene.sets.input <- My_data()$Gene.sets.input
    checkboxGroupInput(ns("geneset"), label = "Geneset PRS to include:",
                       choices = Gene.sets.input, selected = Gene.sets.input)
  })
  
  
  
  output$PvalPlot <- renderPlot({
    
    My_data()
    
    
    
    # Plot the resulting table for comparisons
    p <- ggplot(part_2(), aes(x=score, y=logp, fill = Type, group=Significance_thresholds))
    
    p <- p +
      geom_point(aes(colour = Type))
    
    
    p <- p + scale_x_discrete(labels=levels(part_2()$alterations))
    
    p <- p + facet_grid(. ~ as.double(Significance_thresholds),scales = "free_x", space = "free_x") +
      theme(strip.text.x = element_text(size = 10))
    
    p <- p + geom_hline(aes(yintercept=1.30103), colour = "red", linetype= "solid", alpha = 0.25)
    p <- p + scale_fill_brewer(palette = "Paired")
    p <- p + theme(axis.text.x = element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
    p <- p + ggtitle(part_2()$.id[1])
    #p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5,vjust = 0.5,angle = 90))
    p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5))
    p <- p + ylab(label = expression(-'log'[10]*'(p)'))
    p <- p + xlab(label = "Polygenic risk score")
    p
})
}
    
## To be copied in the UI
# mod_Gene_set_regression_ui("Gene_set_regression_ui_1")
    
## To be copied in the server
# callModule(mod_Gene_set_regression_server, "Gene_set_regression_ui_1")
 
