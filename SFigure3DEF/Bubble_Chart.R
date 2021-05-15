# upload data from "SupplementaryData5_DifferentlyExpressedProteins"
funImport <- function(sheetNo){
  
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/ProteomicsManuscript_03.07.2019/NewFigures&Data_2020/5-SupplementaryData")
  library(readxl)
  
  dat <- read_excel("SupplementaryData6_DifferentlyExpressedProteins.xlsx", sheet = sheetNo)
  dat <- data.frame(dat)
  
  for (i in 1:nrow(dat)){
    dat$Mean_Pvalue_FDR[i] <- mean(as.numeric(dat[i, 5:6]), na.rm = TRUE)
  }
  
  dat <- dat[, c(3, 4, 7)]
  colnames(dat)  <- c("Gene.symbol", "logFC", "adj.P.Val")
  
  return(dat)
  
}

datEM <- funImport(1)
datHE <- funImport(2)
datHM <- funImport(3)

#----------------------------------------------------------------------------------------#

# pathfindR

#----------------------------------------------------------------------------------------#
FunPathfindR <- function(dat){
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/pathfindR/BubbleChart_2020.04.07")
  suppressPackageStartupMessages(library(pathfindR))
  
  # PIN == Protein-Protein-Interaction
  # Biogrid", "STRING", "GeneMania", "IntAct", "KEGG", "mmu_STRING" or custom PIN 
  dat_processed <- input_processing(input = dat, 
                                    p_val_threshold = 0.05, 
                                    pin_name_path  = "Biogrid",
                                    convert2alias = TRUE)
  
  # Number of genes provided in input: 715
  # Number of genes in input after p-value filtering: 715
  
  # Found interactions for all genes in the PIN
  # Final number of genes in input: 715
  
  # using "KEGG" as our gene sets for enrichment
  # The available gene sets in pathfindR are "KEGG", "Reactome", "BioCarta", 
  # "GO-All", "GO-BP", "GO-CC" and "GO-MF", or "Custom"
  GeneSet_list <- fetch_gene_set(gene_sets = "KEGG",
                                 min_gset_size = 10,
                                 max_gset_size = 300)
  GeneSet_gsets <- GeneSet_list[[1]]
  GeneSet_descriptions <- GeneSet_list[[2]]
  
  n_iter <- 100 ## number of iterations
  combined_res <- NULL ## to store the result of each iteration
  
  for (i in 1:n_iter) {
    
    ###### Active Subnetwork Search
    snws_file <- paste0("active_snws_", i) # Name of output file
    active_snws <- active_snw_search(input_for_search = dat_processed, 
                                     pin_name_path = "KEGG", 
                                     snws_file = snws_file,
                                     score_quan_thr = 0.8, # you may tweak these arguments for optimal filtering of subnetworks
                                     sig_gene_thr = 0.02, # you may tweak these arguments for optimal filtering of subnetworks
                                     search_method = "GR")
    # The chosen PIN must be one of: "Biogrid", "STRING", "GeneMania", "IntAct", "KEGG", "mmu_STRING", "/path/to/custom/SIF"
    ###### Enrichment Analyses
    current_res <- enrichment_analyses(snws = active_snws,
                                       sig_genes_vec = dat_processed$GENE,
                                       pin_name_path = "KEGG", 
                                       genes_by_term = GeneSet_gsets,
                                       term_descriptions = GeneSet_descriptions,
                                       adj_method = "bonferroni",
                                       enrichment_threshold = 0.05,
                                       list_active_snw_genes = TRUE) # listing the non-input active snw genes in output
    
    ###### Combine results via `rbind`
    combined_res <- rbind(combined_res, current_res)
    
  }

  ###### Summarize Combined Enrichment Results
  summarized_df <- summarize_enrichment_results(combined_res, 
                                                list_active_snw_genes = TRUE)
  ###### Annotate Affected Genes Involved in Each Enriched Term
  final_res <- annotate_term_genes(result_df = summarized_df, 
                                   input_processed = dat_processed, 
                                   genes_by_term = GeneSet_gsets)
  # "ID", "Term_Description", "Fold_Enrichment", "occurrence", "lowest_p", 
  # "highest_p", "non_Signif_Snw_Genes", "Up_regulated", "Down_regulated" 
  
  #----------------------------------------------------------------------------------#
  
  # interaction chart
  
  #----------------------------------------------------------------------------------#
  
 # visualize_terms(result_df = final_res,
  #                hsa_KEGG = FALSE, # boolean to indicate whether human KEGG gene sets were used for enrichment analysis or not
   #               pin_name_path = "KEGG")
  #----------------------------------------------------------------------------------#
  
  # Bubble chart
  
  #----------------------------------------------------------------------------------#
  source("enrichment_chart.R")
  enrichment_chart(final_res, top_terms = 15)
  
  library(xlsx)
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/pathfindR/BubbleChart_2020.04.07")
  
  file <- paste(format(Sys.time(), "%F %H-%M-%S"), 
                "KEGGannotation_08.04.2020.xlsx", 
                sep = " ")
  
  write.xlsx(final_res, file, 
             row.names = FALSE, 
             sheetName = "KEGG")
}

FunPathfindR(datEM)
FunPathfindR(datHE)
FunPathfindR(datHM)

# Sys.setenv(TZ = "Greenwich Mean Time")


