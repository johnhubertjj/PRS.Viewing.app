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
    sidebarLayout(
    sidebarPanel(
      fileInput(ns("file1"), "Choose an input file",
                multiple = F),

      uiOutput(ns("Significance_threshold")
      ),
      checkboxGroupInput(ns("Gene_regions"), label = "Length of Gene regions:",
                           choices = c("Genome-wide","Gene-set")),
      uiOutput(ns("DSM")),
      uiOutput(ns("geneset")),
      sliderInput('plotHeight', 'Bar which does nothing, use if bored', 
                  min = 100, max = 2000, value = 1000)
      
    ),
    mainPanel(
                tabsetPanel(id = "tabs",
                            
                            tabPanel("Plots",
                                     plotOutput(ns('PvalPlot')),
                                     plotOutput(ns('Beta_plot')),
                                     plotOutput(ns('R2_plot'))),
                            tabPanel("Table", dataTableOutput(ns('summary_table'))),
                            tabPanel("Input variables",
                                     rclipboard::rclipboardSetup(),
                                     
                                     
                                     textAreaInput("text_2",label = "Full message will appear here:",width = "500px",height = "100px", resize = "both",
                                                   placeholder = "Twitter handles will appear here at the end of your message depending on the options selected to the left (eg: Pint of Science is Great! @virustinkerer)")
                                     ,
                                     # UI ouputs for the copy-to-clipboard buttons
                                     uiOutput(ns("clip")))
                                    )
                            )
              )
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
    
    ## Read in data
    Full_data <- data.table::fread(input$file1$datapath)
    Full_data <- data.table(Full_data)
    data.table::setnames(Full_data, old = c("Set","Threshold"), new = c("Genesets", "Significance_thresholds"))
    Full_data[, Gene_regions := Genesets]
    
    # First runthrough of PRS positions 
    # From utils_helpers.R
    Positions_of_PRS_in_table <- Calculate_positions_of_genome_wide_PRS(Full_data,"Genesets")
    
    
    Full_data <- data.table(Full_data)
    #Genome_wide_PRS <- which(Full_data$Gene_regions == "Base")
    #Gene_set_PRS <- which(Full_data$Gene_regions != "Base")
    
    Full_data[,Gene_regions := "NA"]
    Full_data[Positions_of_PRS_in_table$Genome_wide_PRS, Gene_regions := "Genome-wide"]
    Full_data[Positions_of_PRS_in_table$Gene_set_PRS, Gene_regions := "Gene-set"]
    
    
    ## Create arguments to shiny app
    Gene.sets.input <- unique(Full_data$Genesets)
    significance_threshold.input <- unique(Full_data$Significance_thresholds)
    DSM.input <- "Everything"
    
    
    Full_data[Positions_of_PRS_in_table$Genome_wide_PRS, Type := "Whole_genome"]
    Full_data[!Positions_of_PRS_in_table$Genome_wide_PRS, Type:= "Pathway"]
    
    
    # Set all P values to work within the app and add annotations
    if(any(Full_data$P == 0) == T){
      Full_data[,P_altered := P]
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
    Full_data$p_value_text <- paste("p =", scales::scientific(Full_data$P, digits = 2), sep = " ")
    
    if(any(Full_data$p_value_text == "p = 0.0e+00") == TRUE){
      Full_data[p_value_text == "p = 0.0e+00", p_value_text := "p < 1e-300"]
    }
    
    
    # Run a duplication check on the gene-set PRS names
    # From utils_helpers.R
    Full_data <- Duplication_of_gene_sets_check(Data_table = Full_data, Genome_wide_positions = Positions_of_PRS_in_table$Genome_wide_PRS, Significance_thresholds_name = "Significance_thresholds", gene_set_values = Full_data$Genesets)
    Full_data <- data.table(Full_data)
    
    # Re-calculate the positions of whole_genome_PRS
    #Second runthrough of PRS positions
    Positions_of_PRS_in_table <- Calculate_positions_of_genome_wide_PRS(Full_data,"Genesets")
    Full_data <- data.table(Full_data)
    
    # Change names in alterations 
    # From utils_helpers.R
    alterations <- Pathway_cleanup(Full_data$Genesets, Positions_of_PRS_in_table$Genome_wide_PRS)   
    Full_data <- data.table(Full_data)
    
    Full_data[, alterations := alterations]
    
    #Add DSM argument script
    Full_data[, samples.i. := "Everything"]
    
    Full_data[Positions_of_PRS_in_table$Gene_set_PRS, score := paste0(Gene_regions,"_SCORE_",Genesets,"_",Significance_thresholds)]
    Full_data[Positions_of_PRS_in_table$Genome_wide_PRS, score := paste0("All.genome_SCORE_whole_genome_",Significance_thresholds)]
    Full_data[,.id := "PRS Results"]
    
    Full_data <- list(Full_data = Full_data ,Gene.sets.input = Gene.sets.input, significance_threshold.input = significance_threshold.input, DSM.input = DSM.input)
    Full_data
    
  })
  
  part_2 <- reactive({ 
    
    #Gene_regions <- c("Genome-wide", "Gene-set")
    #Significance_thresholds <- c(0.5,1)
    #Genesets <- "GO_CARDIAC_DEVELOPMENT.bed"
    
    Current_table <- My_data()$Full_data %>%
      
      dplyr::filter(samples.i. == input$DSM,
             Gene_regions %in% Gene_region_debounce(),
             Significance_thresholds %in% sigthreshold_debounce(),
             Genesets %in% gene_set_debounce())
    
    Sample_analysis_2 <- as.data.table(Current_table)
    Sample_analysis_2$score <- factor(Sample_analysis_2$score, levels = Sample_analysis_2$score[order(Sample_analysis_2$score, Sample_analysis_2$Type)])
    Sample_analysis_2$alterations <- factor(Sample_analysis_2$alterations, levels = unique(Sample_analysis_2$alterations[order(Sample_analysis_2$score, Sample_analysis_2$Type)]))
    Sample_analysis_2
    
    ## Format DF to DT and apply fixes to the number of decimal points, format "g" = change to nn.dde-dd only if required
    
    ## Okay, so faceting was not meant to have differing axis lables, but in order to place the axis in the right order, I need to specify just one threshold and repeat across all facets
    ## I've used a short-cut here, the line 108 sorts the alterations column by the score and type and then only selects the unique labels for these columns so that the structure is "repeated" across all thresholds
    ## despite not knowing how many thresholds are in the analysis...i've saved a few lines of code and thought here.
    
    
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
  
  # debouncing algorithm -> MUST GO HERE because of the rendering UI scripts above and reactivity scripts have not been processed yet
  sigthreshold_debounce <- reactive({ input$Significance_threshold }) %>% debounce(1000)
  gene_set_debounce <- reactive({ input$geneset }) %>% debounce(1000)
  Gene_region_debounce <- reactive({ input$Gene_regions }) %>% debounce(1000)
  
  output$PvalPlot <- renderPlot({
    
    My_data()
    
    if (is.null(sigthreshold_debounce())) {
      return(NULL)
    }    
    if (is.null(gene_set_debounce())) {
      return(NULL)
    }    
    if (is.null(Gene_region_debounce())) {
      return(NULL)
    }  
    
    
    # Plot the resulting table for comparisons
    p <- ggplot2::ggplot(part_2(), ggplot2::aes(x=score, y=logp, fill = Type, group=Significance_thresholds))
    
    p <- p +
      ggplot2::geom_point(ggplot2::aes(colour = Type))
    
    
    p <- p + ggplot2::scale_x_discrete(labels=levels(part_2()$alterations))
    
    p <- p + ggplot2::facet_grid(. ~ as.double(Significance_thresholds),scales = "free_x", space = "free_x") +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = 10))
    
    p <- p + ggplot2::geom_hline(ggplot2::aes(yintercept=1.30103), colour = "red", linetype= "solid", alpha = 0.25)
    p <- p + ggplot2::scale_fill_brewer(palette = "Paired")
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
    p <- p + ggplot2::ggtitle(part_2()$.id[1])
    #p <- p + theme(plot.title = element_text( face = "bold",hjust = 0.5,vjust = 0.5,angle = 90))
    p <- p + ggplot2::theme(plot.title = ggplot2::element_text( face = "bold",hjust = 0.5))
    p <- p + ggplot2::ylab(label = expression(-'log'[10]*'(p)'))
    p <- p + ggplot2::xlab(label = "Polygenic risk score")
    p
})
  output$Beta_plot <- renderPlot({
    
    My_data()
    if (is.null(sigthreshold_debounce())) {
      return(NULL)
    }    
    if (is.null(gene_set_debounce())) {
      return(NULL)
    }    
    if (is.null(Gene_region_debounce())) {
      return(NULL)
    }    
    
    
    # Put in the code below above, removing all of the excess alterations work to create the pdf plots...
    
    p <- ggplot2::ggplot(part_2(), ggplot2::aes(x=score, y=estimate, fill = Type, group=Significance_thresholds))
    
    p <- p +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = SE_higher, ymax = SE_lower), position = "dodge", width = 0.25) +
      ggplot2::geom_point(ggplot2::aes(colour = Type))
    
    p <- p + ggplot2::scale_x_discrete(labels= levels(part_2()$alterations))
    p <- p + ggplot2::facet_grid(. ~ as.double(Significance_thresholds), scales = "free_x", space = "free_x") +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = 10))
    p <- p + ggplot2::geom_hline(ggplot2::aes(yintercept=0), colour = "red", linetype= "solid", alpha = 0.25)
    p <- p + ggplot2::scale_fill_brewer(palette = "Paired")
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,size = 10,hjust = 1,vjust = 0.5))
    p <- p + ggplot2::ggtitle(part_2()$.id[1])
    p <- p + ggplot2::theme(plot.title = ggplot2::element_text( face = "bold",hjust = 0.5))
    p <- p + ggplot2::ylab(label = "BETA")
    p <- p + ggplot2::xlab(label = "Polygenic risk score")
    p
    
    # ggplotly(p) %>% 
    # layout(height = input$plotHeight, autosize=TRUE)
    
    # Possible improvements:
    # Implement in switch from whole genome to gene-sets
    # Implement data-table of the raw results
    # Implement output file of the plots
    # Colour rows for significant values
    # Incorporate into its own app
    # 
    
  })
  output$R2_plot <- renderPlot({
    
    My_data()
    
    if (is.null(sigthreshold_debounce())) {
      return(NULL)
    }    
    if (is.null(gene_set_debounce())) {
      return(NULL)
    }    
    if (is.null(Gene_region_debounce())) {
      return(NULL)
    }    
    
    p <- ggplot2::ggplot(part_2(), ggplot2::aes(x=score, y=r2_dir, fill = Type, group=Significance_thresholds))
    p <- p +
      ggplot2::geom_bar(stat = "identity", ggplot2::aes(colour = Type), position = "dodge") +
      ggplot2::geom_text(data=subset(part_2(), P < 0.05),
                         ggplot2::aes(x=score,y=r2_dir,label=p_value_text, hjust=ifelse(sign(r2_dir)>0, 0, 0)), angle = 90, position = ggplot2::position_dodge(width = 1), size = 2.9)
    
    #Problem with labels with a workaround
    # I use the score column in the format of factors and reference each relevant dataset for ggplot.
    # However this relies on having 0.05 and 0.5 in the value name.
    # scale_x_discrete accepts functions, but I also need to convert SCORE_0.05 and Score_0.5 into a "Whole_genome_PRS" which is almost impossible to write"
    # However as the labels function accepts key:value pairs, I wrote a vector in R that maps the original names of the pathways to "human readable" format using names function in R
    # This should work for most instances
    
    p <- p + ggplot2::scale_x_discrete(labels= levels(part_2()$alterations))
    p <- p + ggplot2::scale_y_continuous(expand = ggplot2::expand_scale(mult = c(0.2,.6)))
    p <- p + ggplot2::facet_grid(. ~ as.double(Significance_thresholds), scales = "free_x", space = "free_x") +
      ggplot2::theme(strip.text.x = ggplot2::element_text(size = 10))
    p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, size = 10, hjust = 1,vjust = 0.5))
    p <- p + ggplot2::ggtitle(part_2()$.id[1])
    p <- p + ggplot2::theme(plot.title = ggplot2::element_text( face = "bold",hjust = 0.5))
    p <- p + ggplot2::ylab(label = "R2_dir (%)")
    p <- p + ggplot2::xlab(label = "Polygenic risk score")
    p
    
    #ggplotly(p) %>% 
    #  layout(height = input$plotHeight, autosize=TRUE)
    
    # Possible improvements:
    # Implement in switch from whole genome to gene-sets
    # Implement data-table of the raw results
    # Implement output file of the plots
    # Colour rows for significant values
    # Incorporate into its own app
    # 
    
  })
  output$summary_table <- renderDataTable({
    
    My_data()
    #browser()
    # These are required in case no tick boxes are selected
    if (is.null(input$Significance_threshold)) {
      return(NULL)
    }    
    if (is.null(input$geneset)) {
      return(NULL)
    }    
    if (is.null(input$Gene_regions)) {
      return(NULL)
    }    
    
    # Select columns you wish to output
    cols <- c("estimate", "SE","R2","P", "Num_SNP")
    
    ## Limit data table to input arguments and pipe to limiting columns and ordering based on significance
    sample_analysis <- My_data()$Full_data %>%
      dplyr::filter(samples.i. == input$DSM,
             Gene_regions %in% Gene_region_debounce(),
             Significance_thresholds %in% sigthreshold_debounce(),
             Genesets %in% gene_set_debounce()
      )  %>%
      dplyr::select(c(Genesets,Gene_regions,Significance_thresholds,estimate,SE,P,R2,Num_SNP)) %>%
      dplyr::arrange(P)
    
    ## Format DF to DT and apply fixes to the number of decimal points, format "g" = change to nn.dde-dd only if required
    Sample_analysis_2 <- data.table::as.data.table(sample_analysis)
    Sample_analysis_2[, (cols) := lapply(.SD, formatC, digits = 3, format = "g"), .SDcols = cols]
    
    ## leave datatable function in for "prettyfying" the final result    
    DT::datatable(data = Sample_analysis_2,
              options = list(pageLength = 10),
              rownames = F)
    
    # Possible improvements:
    # colour rows for significant values
    # incorporate into its own app
    # 
  })
  observe({
    
    
    if (!is.null(input$Significance_threshold)) {
      Sig_print <- paste(input$Significance_threshold, collapse = ",")
      
    } else{
      Sig_print <- "placeholder"
    }  
    
    if (!is.null(input$geneset)) {
      geneset_print <- paste("\"",input$geneset,"\"", collapse = ",", sep = "")
    }else{
      geneset_print <- "placeholder"
    }    
    
    if (!is.null(input$Gene_regions)) {
      Generegion_print <- paste("\"",input$Gene_regions,"\"", collapse = ",", sep="")
    } else{
      Generegion_print <- "placeholder"
    }   
    
    if (!is.null(input$DSM)) {
      DSM_print <- paste0("\"",input$DSM,"\"")
    }else{
      DSM_print <- "placeholder"
    }
    
    if (!is.null(input$file1$datapath)){
      data_path_print <- paste(input$file1$datapath)
    }else{
      data_path_print <- "path_not_found"
    }
    
    output_text <- paste0( "input <- list() \n",
                           "\n input$Significance_threshold <- c(", Sig_print,")",
                           "\n input$geneset <- c(", geneset_print,")", 
                           "\n input$Gene_regions <- c(", Generegion_print,")", 
                           "\n input$DSM <- ", DSM_print, 
                           "\n data_print_path <- ", data_path_print)
    
    
    
    updateTextInput(session,"text_2",value = output_text)
    
    #number_of_characters <- paste(text_output_speaker_2, collapse = " ")
    #number_of_characters <- nchar(input$text_2)
    #output$length_text_left <- renderText(280 - number_of_characters)
  })
  
  output$clip <- renderUI({
    rclipboard::rclipButton("clipbtn", "Copy to Clipboard", input$text_2, icon("clipboard"))
  })
  
  
  
  
}
    
## To be copied in the UI
# mod_Gene_set_regression_ui("Gene_set_regression_ui_1")
    
## To be copied in the server
# callModule(mod_Gene_set_regression_server, "Gene_set_regression_ui_1")
 
