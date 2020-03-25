# Clean up names of PRS scripts

#Data_table <- Full_data
#Genome_wide_positions <- whole_genome_plot_all_positions
#gene_set_column_name <- "Genesets"
#Significance_thresholds_name <- "Significance_thresholds"
#gene_set_values <- Full_data$Genesets

# Duplications of the gene-set in the table
Duplication_of_gene_sets_check <- function(Data_table,Genome_wide_positions, Significance_thresholds_name, gene_set_values){
  
  ## Identify the number of gene sets and the significance thresholds input into PRSAVE
  number_of_unique_significance_thresholds <- length(unique(Data_table[[Significance_thresholds_name]]))
  number_of_PRS_gene_set <- length(unique(gene_set_values[-Genome_wide_positions]))
  
  ## Calculate whether values to compare against
  expected_number_of_gene_set_names <- number_of_unique_significance_thresholds * number_of_PRS_gene_set
  actual_number_of_gene_set_names <- length(gene_set_values[-Genome_wide_positions])
  
  ## If statement checking that there are any duplications of gene-set names
  if (expected_number_of_gene_set_names != actual_number_of_gene_set_names){
    
    ## Dplyr solution to figuring out whether the rows within the data frame are duplicated or not
    Data_table_new <- Data_table %>%
      dplyr::distinct()
  }
  Data_table_new
}

Pathway_cleanup <- function(alterations, Genome_wide_positions){
  alterations[-Genome_wide_positions] <- tolower(alterations[-Genome_wide_positions])
  alterations[Genome_wide_positions] <- "Whole Genome PRS ALL"
  alterations
}

Calculate_positions_of_genome_wide_PRS <- function(Data_table, Gene_regions_column){
  Type_of_PRS_positions <- list()
  
  Type_of_PRS_positions[["Genome_wide_PRS"]] <- which(Data_table[[gene_set_column_name]] == "Base")
  Type_of_PRS_positions[["Gene_set_PRS"]] <- which(Data_table[[gene_set_column_name]] != "Base")
  
  Type_of_PRS_positions
}