setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/IPA/ImportFile")
library(readxl)
em <- read_excel("Proteomics_LFQ&SILAC_EM.xlsx", sheet = 1)
em <- data.frame(em) # 883
em <- em[which(em$p.value < 0.05), ] # 728

he <- read_excel("Proteomics_LFQ&SILAC_HE.xlsx", sheet = 1)
he <- data.frame(he) # 672
he <- he[which(he$p.value < 0.05), ] # 548

table(duplicated(em$Gene.Name)) # 13 duplicates
table(duplicated(he$Gene.Name)) # 13 duplicates

em <- em[-which(duplicated(em$Gene.Name)), ] # 715
he <- he[-which(duplicated(he$Gene.Name..GN.)), ] # 535

#----------------------------------------------------------------------------#
# double check the LFQ FDR < 0.05
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/pathfindR/Data/LFQ")
lfq <- read_excel("5 Dundee_LFQ_valid_Log2_imputated_p.value_t.stat_q.value_mean_22.11.2018.xlsx", sheet = 1)
lfq <- data.frame(lfq)

lfq_em <- lfq[which(lfq$Gene.Name..GN. %in% em$Gene.Name), ] # 722
lfq_em <- lfq_em[-which(duplicated(lfq_em$Gene.Name..GN.)), ] # 715

lfq_he <- lfq[which(lfq$Gene.Name..GN. %in% he$Gene.Name..GN.), ] # 537
lfq_he[which(duplicated(lfq_he$Gene.Name..GN.)), ] # LMNA, ANXA6
grep("LMNA", lfq_he$Gene.Name..GN.) # 49, 276
grep("ANXA6", lfq_he$Gene.Name..GN.) # 5, 522
lfq_he <- lfq_he[c(-49, -522), ] # 535

table(lfq_em$EM_significant) # 715
table(lfq_he$HE_significant) # 535

#---------------------------------------------------------------------#
# use LFQ data to do the "pathfindR"
lfq_em$"logFC" <- -lfq_em$Mean_Log2.E.M. # D492/D492M -> D492M/D492

path_em <- lfq_em[, c(3, 31, 15)]
colnames(path_em) <- c("Gene.symbol", "logFC", "adj.P.Val")

path_he <- lfq_he[, c(3, 24, 21)]
colnames(path_he) <- c("Gene.symbol", "logFC", "adj.P.Val")

#--------------------------------------------------------------------------#
# for heatmap
rownames(lfq_em) <- lfq_em$Gene.Name..GN.
rownames(lfq_he) <- lfq_he$Gene.Name..GN.

exp_mat_em <- lfq_em[, c(10:12, 4:6, 7:9)]
exp_mat_he <- lfq_he[, c(10:12, 4:6, 7:9)]

exp_mat_em <- as.matrix(exp_mat_em)
exp_mat_he <- as.matrix(exp_mat_he)

#--------------------------------------------------------------------------#
# pathfindR, do EM and HE one at a time, note EM is actually ME
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/pathfindR/GOannotation/EM&HE&HM_03022019/Version_2020.01")
suppressPackageStartupMessages(library(pathfindR))

# PIN == Protein-Protein-Interaction
# pin_path <- return_pin_path(pin_name_path = "Biogrid") 
# Biogrid", "STRING", "GeneMania", "IntAct", "KEGG", "mmu_STRING" or custom PIN 
em_processed <- input_processing(input = path_em, 
                                 p_val_threshold = 0.05, 
                                 pin_name_path  = "Biogrid",
                                 convert2alias = TRUE)

# Number of genes provided in input: 715
# Number of genes in input after p-value filtering: 715
# pathfindR cannot handle p values < 1e-13. These were changed to 1e-13

# Found interactions for all genes in the PIN
# Final number of genes in input: 715

# using "BioCarta" as our gene sets for enrichment
# The available gene sets in pathfindR are "KEGG", "Reactome", "BioCarta", 
# "GO-All", "GO-BP", "GO-CC" and "GO-MF", or "Custom"
biocarta_list <- fetch_gene_set(gene_sets = "BioCarta",
                                min_gset_size = 10,
                                max_gset_size = 300)
biocarta_gsets <- biocarta_list[[1]]
biocarta_descriptions <- biocarta_list[[2]]

n_iter <- 100 ## number of iterations
combined_res <- NULL ## to store the result of each iteration

for (i in 1:n_iter) {
  
  ###### Active Subnetwork Search
  snws_file <- paste0("active_snws_", i) # Name of output file
  active_snws <- active_snw_search(input_for_search = em_processed, 
                                   pin_name_path = "KEGG", 
                                   snws_file = snws_file,
                                   score_quan_thr = 0.8, # you may tweak these arguments for optimal filtering of subnetworks
                                   sig_gene_thr = 0.02, # you may tweak these arguments for optimal filtering of subnetworks
                                   search_method = "GR")
  # The chosen PIN must be one of: "Biogrid", "STRING", "GeneMania", "IntAct", "KEGG", "mmu_STRING", "/path/to/custom/SIF"
  ###### Enrichment Analyses
  current_res <- enrichment_analyses(snws = active_snws,
                                     sig_genes_vec = em_processed$GENE,
                                     pin_name_path = "KEGG", 
                                     genes_by_term = biocarta_gsets,
                                     term_descriptions = biocarta_descriptions,
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
                                 input_processed = em_processed, 
                                 genes_by_term = biocarta_gsets)
# "ID", "Term_Description", "Fold_Enrichment", "occurrence", "lowest_p", 
# "highest_p", "non_Signif_Snw_Genes", "Up_regulated", "Down_regulated" 

# KEGG pathways
visualize_hsa_KEGG(summarized_df$ID, em_processed)

visualize_terms(result_df = final_res,
                hsa_KEGG = FALSE, # boolean to indicate whether human KEGG gene sets were used for enrichment analysis or not
                pin_name_path = "Reactome")

# Bubble chart
source("./FunpathfindR/enrichment_chart.R")
enrichment_chart(final_res[which(final_res$Fold_Enrichment >= 1), ], top_terms = 15)

library(xlsx)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/pathfindR/GOannotation/EM&HE&HM_03022019Version_2020.01")
write.xlsx(final_res, "KEGGannotation_HE_08.01.2020.xlsx", 
           row.names = FALSE, 
           sheetName = "KEGG_HE")
# enrichment was not run cause it takes too long
# later, wanted to run GO-BP and GO-CC separately, see if the enrichment can take shorter time
# but didnt finish, just did GO-BP not GO-CC (on 04.02.2019)

