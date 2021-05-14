# http://www.sthda.com/english/wiki/survival-analysis-basics
# http://www.sthda.com/english/wiki/survminer-r-package-survival-data-analysis-and-visualization
# https://www.cbioportal.org/datasets (Breast Invasive Carcinoma (TCGA, Provisional))
# DEseq
# FPKM, Counts
# https://towardsdatascience.com/kaplan-meier-curves-c5768e349479
#-----------------------------------------------------------------------------------------#

                                # install all the packages

#-----------------------------------------------------------------------------------------#
#devtools::install_github('siggitrausti/siggitRausti',force=T)
#library(BiocManager)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager") # update to new version

#BiocManager::install("AnnotationHub", version = "3.8")
#BiocManager::install("clusterProfiler", version = "3.8")
#BiocManager::install("org.Hs.eg.db", version = "3.8")
#BiocManager::install("topGO", version = "3.8")
#BiocManager::install("DOSE", version = "3.8")

#install.packages("factoextra")
#install.packages("NbClust")
#install.packages("survminer")
#install.packages("randomForestSRC")
#install.packages("ggRandomForests")

#-------------------------------------------------------------------------------------------#

                                      # Survival plotting

#-------------------------------------------------------------------------------------------#
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/SurvivalCurveGeneration_17.09.2019")
library(dendextend)
library(colorspace)
library(factoextra)
library(NbClust)
library(gplots)
library(tidyr)
library(dplyr)
library(survminer) # useful
library(survival) # useful
library(BiocManager)
# library(clusterProfiler)
# library(org.Hs.eg.db)
# library(topGO)
# library(DOSE)
# library(ReactomePA)
library(ggforce)
library(igraph)
library(siggitRausti) # useful
# library(DESeq2)
library(genefilter)
library(randomForestSRC)
library(ggRandomForests)
library(dplyr)
library(randomForest)
# library(genefu)

source('./Rfiles/SurvivalSignatureTest/findBestSurvival2_qiong.R')
source('./Rfiles/SurvivalSignatureTest/segregate_dataset_qiong.R')
# source('./Rfiles/SurvivalSignatureTest/sample_vector.R') # useless
# source('./Rfiles/SurvivalSignatureTest/enrichment_test.R') # useless
# source('./Rfiles/SurvivalSignatureTest/test_match_order.R') # useless
# source("./Rfiles/SurvivalSignatureTest/ColBrew.R") # useless

# col5 <- ColBrew('JCO')

#------------------------------------------------------------------------------------------#

                                          # import data

#------------------------------------------------------------------------------------------#
# First, use the actual list of interest:
# subtype, survival ("sub_vec", "age", "patient", "status", "survival", "cluster", "ER")
dd_subtype <- readRDS(file = "Subtyped_TCGA_data.rds") # 1066 * 7

# all genes for all cases (20247 * 1066)
# NON-z-scored TCGA data (but normalized)
vsd_tcga <- readRDS(file = 'TCGA_data_whole_vst_transformed_170519.rds') #DEseq package

# clinical data with living/dead status
clinical_data <- read.table('data_bcr_clinical_data_patient.txt',header=T,sep='\t',fill=T) # 1072  * 115
clin_data <- clinical_data[,c("PATIENT_ID", 
                              "OS_STATUS", 
                              "OS_MONTHS",
                              "DFS_STATUS", 
                              "DFS_MONTHS",
                              "RACE")] # 1072 * 6
clin_data$PATIENT_ID <- gsub("-",".",clin_data$PATIENT_ID)
clin_data$PATIENT_ID <- paste0(clin_data$PATIENT_ID,'.01')
rm(clinical_data)

#-------------------------------------------------------------------------------------------#

                              # choose breast cancer subtypes

#-------------------------------------------------------------------------------------------#
# valid_patients <- dd_subtype$patient # all subtypes
# valid_patients <- dd_subtype$patient[which(dd_subtype$sub_vec == 'Basal-like')] # 198
# valid_patients <- dd_subtype$patient[which(dd_subtype$sub_vec == 'HER2-enriched')] # 91, for testing of other subtypes with methodology
# valid_patients <- dd_subtype$patient[which(dd_subtype$sub_vec == 'Luminal A')] # 395
# valid_patients <- dd_subtype$patient[which(dd_subtype$sub_vec == 'Normal-like')] # 15
valid_patients <- dd_subtype$patient[which(dd_subtype$sub_vec == 'Luminal B')] # 367

# select cases from vsd_tcga based on cases with specific subtypes
id_pats <- which(colnames(vsd_tcga) %in% valid_patients)
basal_cor_data <- vsd_tcga[, id_pats]
rm(valid_patients)
rm(id_pats)

#-------------------------------------------------------------------------------------------#

                                 # READ IN GENES OF INTEREST

#-------------------------------------------------------------------------------------------#
# basal_lethal <- data.frame("GFPT2")
# basal_lethal <- data.frame(c("GALE", "PGM2L1", "UGDH", "GFPT2")) # c("GALE", "PGM2L1", "UGDH", "GFPT2")
# basal_lethal <- data.frame(c("CYP1B1","PGM3", "PGM2L1", "UGDH", "ADPGK", "NQO1", "DGKA", "GLUL"))
# basal_lethal <- read.delim("HE.txt", header = FALSE) # 535 proteins with significant differences
# basal_lethal <- read.delim("EM.txt", header = FALSE) # 715 proteins with significant differences
# basal_lethal <- read.delim("MetabolicTargetsHE.txt", header = FALSE) # 12 Metabolic enzymes for HE
# basal_lethal <- read.delim("MetabolicTargetsEM.txt", header = FALSE) # 22 Metabolic enzymes for EM
# basal_lethal <- data.frame(c("GFPT2","UPP1", "GSTM1", "SULT1E1", "PXDN"))

# basal_lethal <- read.delim("MetabolicTargetsFromProteomicDataAnalysis.txt", header = FALSE) 
# basal_lethal <- data.frame(unique(basal_lethal[, 1])) # 29, and 5 out 29 genes are not in the gene expression dataset
# 5 genes are ACSL1, ACAA2, ITPK1, GSTM1 and PXDN

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/ProteomicsManuscript_03.07.2019/NewFigures&Data_2020/5-SupplementaryData")
library(readxl)
datH <- read_excel("SupplementaryData3_SignatureProteins.xlsx", sheet = 1)
basal_lethal <- data.frame(datH$`Gene Name`)

# basal_lethal <- unique(basal_lethal)
colnames(basal_lethal) <- 'GeneName'

#------------------------------------------------------------------------------------------#

                                         # analysis

#------------------------------------------------------------------------------------------#
# OS_STATUS (overall survival), change to 0/1
id_living <- which(clin_data$OS_STATUS == "LIVING") # 912
id_dead <- which(clin_data$OS_STATUS == "DECEASED") # 151

status <- rep(NA,dim(clin_data)[1]) # 1072
status[id_living] <-  0
status[id_dead] <-  1

rm(id_living)
rm(id_dead)

# DFS_STATUS (disease-free survival), change to 0/1
id_living2 <-  which(clin_data$DFS_STATUS == "DiseaseFree") # 859
id_dead2 <-  which(clin_data$DFS_STATUS == "Recurred/Progressed") # 113

dfs_status <- rep(NA,dim(clin_data)[1]) # 1072
dfs_status[id_living2] <-  0
dfs_status[id_dead2] <-  1

rm(id_living2)
rm(id_dead2)

# which rows in clin_data are in basal_cor_data (1066/1072)
index_vector <- c()
for (i in 1:dim(basal_cor_data)[1]){
  idx <- which(as.character(clin_data$PATIENT_ID) == colnames(basal_cor_data)[i])
  if (length(idx) != 0){
    index_vector <- c(index_vector,idx)
  }
}
idx_clind <- index_vector # 1066
rm(index_vector)

# patient ID, OS status (0/1), OS months, DFS status (0/1), DFS months
months <- as.integer(as.character(clin_data$OS_MONTHS[idx_clind])) # NAs introduced by coercion 
dfs_months <- as.integer(as.character(clin_data$DFS_MONTHS[idx_clind])) # NAs introduced by coercion
cor_data2 <- data.frame(cbind(clin_data$PATIENT_ID[idx_clind],
                              status[idx_clind],
                              months,
                              dfs_status[idx_clind],
                              dfs_months))

# change colnames and change from "factor" to "numeric"
cor_data3 <- cor_data2
rm(months)
rm(dfs_months)
rm(cor_data2)
colnames(cor_data3) <- c('patient','status','survival','dfs_status','dfs_survival')
cor_data3$survival <- as.numeric(as.character(cor_data3$survival))
cor_data3$status <- as.numeric(as.character(cor_data3$status))
cor_data3$dfs_survival <- as.numeric(as.character(cor_data3$dfs_survival))
cor_data3$dfs_status <- as.numeric(as.character(cor_data3$dfs_status))

# change dataset names for the "gene expression data" and the "survival data"
basal_data <- vsd_tcga[,which(colnames(vsd_tcga) %in% colnames(basal_cor_data))]
rm(basal_cor_data) # dataset "basal_data" is the same as dataset "basal_cor_data"
basal_patients <- cor_data3
rm(cor_data3)

#---------------------------------------------------------------------------------------------#

                                           # p values

#---------------------------------------------------------------------------------------------#
# targeted genes and p values
genes_present_basal <- c()
p_val_basal <- c()
for (i in 1:length(basal_lethal$GeneName)){
  if (length(which(basal_lethal$GeneName[i] %in% rownames(basal_data))) == 0){
    next
  } else {
    genes_present_basal <- rbind(genes_present_basal,as.character(basal_lethal$GeneName[i]))
    p_val_basal <- rbind(p_val_basal,
                         segregate_dataset(basal_data, # gene expression data
                                           basal_lethal$GeneName[i], # gene list (one by one)
                                           0.20, # 0.10
                                           basal_patients, # survival data
                                           c(), # "print_plot" has to be "FALSE" because of c().
                                           print_plot = FALSE)) 
  } # function "segregate_dataset" for calculating p values (print_plot = FALSE) and plotting (print_plot = TRUE & there are effects)
} # one p value for one gene

# FDR
# p_val_basal_adj <- p.adjust(p_val_basal,'fdr',length(p_val_basal)) # NOT used later
lethal_genes_sign_basal <- data.frame(genes_present_basal[which(p_val_basal < 0.05)]) # actually, no genes are correlated after
# correction for multiple testing. So... just use the non-corrected ones. 
colnames(lethal_genes_sign_basal) <- 'lethal_basal'

rm(genes_present_basal)
rm(p_val_basal)
print(lethal_genes_sign_basal)

#----------------------------------------------------------------------------------------#

                                   # KM plotting 

#----------------------------------------------------------------------------------------#
# find the effects ("Beneficial", "Harmful or "No effect") - function "findBestSurvival2"
basal_signature <- basal_lethal
colnames(basal_signature) <- 'basal'
effect_type <- rep(0, dim(basal_signature)[1])

# function "findBestSurvival2" is in package "siggitRausti", p values (9 replicates) & effect
# the reason for having 9 p values is because of seq(0.3, 0.7, by = 0.05).
for (i in 1:dim(basal_signature)[1]){
  effect_type[i] <- findBestSurvival2(basal_data, # gene expression data
                                      basal_signature[i, ], # gene list (one by one)
                                      seq(0.3, 0.7, by = 0.05),
                                      basal_patients, # survival data
                                      print_plot = F, # TRUE to plot the KM plot for each gene in each loop
                                      effect_type = T) # if TRUE, then "print_plot" has to be "FALSE"
} 

basal_sign_final <- data.frame(cbind(basal_signature,
                                     effect_type))
colnames(basal_sign_final) <- c('genes','effect')

print(basal_sign_final)

# Try segregating:
basal_fin <- segregate_dataset(basal_data, # gene expression data
                               basal_sign_final[, 1], # gene list
                               0.20,
                               basal_patients, # survival data
                               basal_sign_final[, 2],
                               print_plot = T) # "print_plot" can be "TRUE" because of the effect.
basal_fin

#----------------------------------------------------------------------------------------------#

                      # export data with "harmful" or "beneficial" effects

#----------------------------------------------------------------------------------------------#
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/SurvivalCurveGeneration_17.09.2019/Outputs")
library(xlsx)
write.xlsx(basal_sign_final,
           file = "Effect_Basal_Sign_Final.xlsx", 
           sheetName = "Effects", 
           row.names = FALSE)



#----------------------------------------------------------------------------------------#

                              # lethal genes for each case

#----------------------------------------------------------------------------------------#
cor_data_fin <- data.frame(cbind(basal_patients$status,
                                 basal_patients$survival))
rownames(cor_data_fin) <- basal_patients$patient
colnames(cor_data_fin) <- c('status','survival')

# HERE I NEED TO ONLY HAVE THE GENES WITH SIGNIFICANT CORRELATION WITH SURVIVAL. EXTRACT THOSE BEFOREHAND. 
for (i in 1:length(lethal_genes_sign_basal$lethal_basal)){
  hehe <- prepareSurvivalDataTCGA(basal_data,
                                  data.frame(lethal_genes_sign_basal$lethal_basal[i]),
                                  0.5,
                                  0.99) # 1066 of 0s or 1s for each "lethal gene", z-score
  cor_data_fin <- cbind(cor_data_fin,
                        as.integer(hehe))
  colnames(cor_data_fin)[ncol(cor_data_fin)] <- as.character(lethal_genes_sign_basal$lethal_basal[i])
} # function "prepareSurvivalDataTCGA" is in package "siggitRausti"

# Take out the NANvalues (1059/1066)
cor_data_fin <- cor_data_fin[-c(which(is.na(cor_data_fin$status))),] 
