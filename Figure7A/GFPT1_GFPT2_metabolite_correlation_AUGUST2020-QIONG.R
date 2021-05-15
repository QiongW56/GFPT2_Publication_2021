# A script for Ottar and Qiong. 

# Goal: Check what metabolites in Ortmayr et al. correlate with GFPT1 and GFPT2 expression within
# the NCI-60 Human tumor cell line panel. 
# library(BiocManager)
# BiocManager::install('rcellminer')
library(rcellminer)
library(readxl)

#----------------------------------------------------------------------#

# Read in metabolite data:

#----------------------------------------------------------------------#
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/NCI60/Datasets")
met_means <- read_excel('41467_2019_9695_MOESM2_ESM.xlsx','Met_alpha','A1:BG2185')
met_sds <- read_excel('41467_2019_9695_MOESM2_ESM.xlsx','Met_alpha_SE','A1:BG2185')
colnames(met_means) <- unlist(as.character(met_means[3,]))
colnames(met_sds) <- unlist(as.character(met_sds[3,]))
met_means <- met_means[-c(1:3),]
met_sds <- met_sds[-c(1:3),]

#----------------------------------------------------------------------#

# Make adjustments to tissue names so that they will match the ones from the NCI60:

#----------------------------------------------------------------------#
expr_dat <- getAllFeatureData(rcellminerData::molData)[["exp"]]
expr_dat <- data.frame(expr_dat)
colnames(met_means) <- gsub('Breast_','BR.',colnames(met_means))
colnames(met_means) <- gsub('CNS_','CNS.',colnames(met_means))
colnames(met_means) <- gsub('Colon_','CO.',colnames(met_means))
colnames(met_means) <- gsub('Lung_','LC.',colnames(met_means))
colnames(met_means) <- gsub('Melanoma_','ME.',colnames(met_means))
colnames(met_means) <- gsub('Ovarian_','OV.',colnames(met_means))
colnames(met_means) <- gsub('Prostate_','PR.',colnames(met_means))
colnames(met_means) <- gsub('Renal_','RE.',colnames(met_means))
colnames(met_means) <- gsub('\\.','',colnames(met_means))
colnames(expr_dat) <- gsub('\\.','',colnames(expr_dat))

#----------------------------------------------------------------------#

# Find if GFPT1 and/or GFPT2 are in the gene expression data:

#----------------------------------------------------------------------#
grep('GFPT',rownames(expr_dat))

# They are both there. Extract only them:
gfpt_dat <- expr_dat[grep('GFPT',rownames(expr_dat)),]

# Create a dataframe of common cell lines from A) Gene expression and B) Metabolite abundance:
metabolite_dataframe <- data.frame(matrix(NA,nrow = 54,ncol = nrow(met_means)))
rownames(metabolite_dataframe) <- colnames(met_means)[6:ncol(met_means)]
for (i in 1:nrow(met_means)){
  metabolite_dataframe[,i] <- as.numeric(as.character(unlist(met_means[i,6:ncol(met_means)])))
  colnames(metabolite_dataframe)[i] <- make.names(met_means$Annotation[i])
}

# Find the index of the cell lines of the gene expression within the metabolite dataframe:
id_vec <- c()
for (i in 1:nrow(metabolite_dataframe)){
  if (length(which(colnames(gfpt_dat) %in% rownames(metabolite_dataframe)[i])) != 0){
    id_temp <- which(colnames(gfpt_dat) %in% rownames(metabolite_dataframe)[i])
    id_vec <- c(id_vec,id_temp)
  } else {
    id_vec <- c(id_vec,NA)
  }
}

# Transpose the gfpt_dat:
gfpt_dat2 <- t(gfpt_dat)

# Combine data:
total_data <- data.frame(cbind(gfpt_dat2[id_vec,]),metabolite_dataframe)
total_data <- total_data[!is.na(total_data[,1]),]

#----------------------------------------------------------------------#

# Find which columns (from columns 3 and above) correlate with GFPT1:

#----------------------------------------------------------------------#
signif_cor <- c()
max_cor <- 0.5
for (i in 3:ncol(total_data)){
  cor_temp_p <- cor.test(total_data[,1],total_data[,i],method = 'spearman',exact = F)[3] # the p-value from the correlation
  cor_temp_val <- abs(cor.test(total_data[,1],total_data[,i],method = 'spearman',exact=F)$estimate)
  if (cor_temp_p < 0.05){
    signif_cor <- c(signif_cor,i)
  }
  if (cor_temp_val > max_cor){
    max_cor <- cor_temp_val
  }
}
            
# Metabolites that correlate with GFPT1 (p < 0.05):
mets_GFPT1 <- colnames(total_data)[signif_cor]
mets_GFPT1

#----------------------------------------------------------------------#

# Find which columns (from columns 3 and above) correlate with GFPT2:

#----------------------------------------------------------------------#
signif_cor2 <- c()
max_cor2 <- 0.5
max_cor2_id <- NA
for (i in 3:ncol(total_data)){
  cor_temp_p <- cor.test(total_data[,2],total_data[,i],method = 'spearman')[3]
  cor_temp_val <- abs(cor.test(total_data[,2],total_data[,i],method = 'spearman')$estimate)
  if (cor_temp_p < 0.05){
    signif_cor2 <- c(signif_cor2,i)
  }
  if (cor_temp_val > max_cor2){
    max_cor2 <- cor_temp_val
    max_cor2_id <- i
  }
}

# Metabolites that correlate with GFPT1 (p < 0.05):
mets_GFPT2 <- colnames(total_data)[signif_cor2]
mets_GFPT2

#----------------------------------------------------------------------#

# Find if glutathione is a part of the correlated metabolites:

#----------------------------------------------------------------------#
mets_GFPT1[grep('Glutathione',mets_GFPT1,ignore.case = T)]
mets_GFPT2[grep('Glutathione',mets_GFPT2,ignore.case = T)]

# It correlates significantly with GFPT2. Now plot it:
plot(total_data[,2],total_data[,grep('Glutathione',colnames(total_data),ignore.case = F)],
     xlab = 'GFPT2 levels',ylab = 'Glutathione levels')
cor.test(total_data[,2],total_data[,grep('Glutathione',colnames(total_data),ignore.case = F)],method = 'spearman',exact = F)

# Now plot GFPT1 (for comparison):
plot(total_data[,1],total_data[,grep('Glutathione',colnames(total_data),ignore.case = F)],
     xlab = 'GFPT1 levels',ylab = 'Glutathione levels')
cor.test(total_data[,1],total_data[,grep('Glutathione',colnames(total_data),ignore.case = F)],method = 'spearman',exact = F)

#----------------------------------------------------------------------#

# Plot two images (using ggplot) and combine into a panel:

#----------------------------------------------------------------------#
library(ggplot2)
# library(siggitRausti)
library(ggpubr)

total_data2 <- total_data[,c(2,which(colnames(total_data) %in% 'Glutathione'))]
total_data2$GFPT2_grouped <- ifelse(total_data2$GFPT2 < 0,'GFPT2_Low','GFPT2_High')
total_data2$GFPT2_grouped <- as.factor(total_data2$GFPT2_grouped)
dat <- total_data2

# plotting
# function to change the space between keys in legends
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.6, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}
GeomBar$draw_key = draw_key_polygon3

yna <- paste("Reduced",colnames(dat)[2])
p <- ggplot(dat, 
            aes(x = GFPT2_grouped, 
                y = Glutathione, 
                fill = GFPT2_grouped)) +
  stat_boxplot(geom = 'errorbar', 
               position = position_dodge(1),
               width = 0.5,
               size = 2) +
  geom_boxplot(position = position_dodge(1), 
               lwd = 3, 
               fatten = 1) +
  #             outlier.colour = "red", outlier.shape = 8, outlier.size = 4) +
  scale_fill_manual(values = c("orange1", "#FAE48BFF")) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               position = position_dodge(1),
               binwidth = 0,
               dotsize = 10) +
  scale_y_continuous(breaks = seq(0, max(dat[, 2]), 25),
                     labels = scales::number_format(accuracy = 1, 
                                                    decimal.mark = ".")) +
 # scale_x_discrete(breaks = c("GFPT2_High","GFPT2_Low"),
 #                   labels = c(expression(GFPT2^High), expression(GFPT2^Low))) +
  theme_classic() +
  labs(title = NULL, # title = paste0("Correlation between ", colnames(dat)[1], " and ", colnames(dat)[2]), 
       x = NULL, 
       y = yna) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 28), 
        axis.line = element_line(colour = "black", size = 3),
        legend.position = "right",
        legend.key.size = unit(1.5, "cm"),
        legend.text = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        axis.title.y.left = element_text(size = 30, face = "bold", 
                                         margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x.bottom = element_text(size = 30, face = "bold",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text = element_text(size = 25, face = "bold", color = "black"),
        axis.ticks.y.left = element_line(colour = "black", size = 2))
p

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/NCI60/Figures")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), 
              colnames(dat)[1],
              "NCI60_Glutathione.Box.plot.tiff", 
              sep = " ")
ggsave(filename = file, units = "in", dpi = 72, width = 10, height = 6)

print(p)

#----------------------------------------------------------------------#

# Sigi's plot, changed to my style

#----------------------------------------------------------------------#
# Boxplot:
my_comparisons <- list(c("GFPT2_High", "GFPT2_Low"))
p <- ggplot(total_data2, aes(x = GFPT2_grouped, y = Glutathione, color = GFPT2_grouped)) +
  #geom_point(size = 3) +
  geom_boxplot(size = 1) +
  stat_compare_means(comparisons = my_comparisons,label = 'p.signif',size=12) +
  scale_x_discrete(
    breaks=c("GFPT2_High","GFPT2_Low"),
    labels=c(expression(GFPT2^High), expression(GFPT2^Low))
  )
p <- plot_look_for_paper(p,y_text = 'Glutathione level',x_text = '',legend_check = F,rotate_check = T) + 
  scale_y_continuous(breaks = seq(0, 175, by=50), limits=c(0,175))
p <-p + theme(panel.grid.major.y = element_line(color='grey'))
p

#ggsave('GFPT2_Glutathione_AUGUST2020.png',width = 6, height = 6, units = "in")

#----------------------------------------------------------------------#

# Now for GFPT1:

#----------------------------------------------------------------------#
total_data3 <- total_data[,c(1,which(colnames(total_data) %in% 'Glutathione'))]
total_data3$GFPT1_grouped <- ifelse(total_data3$GFPT1 < 0,'GFPT1_Low','GFPT1_High')
total_data3$GFPT1_grouped <- as.factor(total_data3$GFPT1_grouped)
dat <- total_data3

# plotting
# function to change the space between keys in legends
draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.6, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}
GeomBar$draw_key = draw_key_polygon3

yna <- paste("Reduced",colnames(dat)[2])

p <- ggplot(dat, 
            aes(x = GFPT1_grouped, 
                y = Glutathione, 
                fill = GFPT1_grouped)) +
  stat_boxplot(geom = 'errorbar', 
               position = position_dodge(1),
               width = 0.5,
               size = 2) +
  geom_boxplot(position = position_dodge(1), 
               lwd = 3, 
               fatten = 1) +
  #             outlier.colour = "red", outlier.shape = 8, outlier.size = 4) +
  scale_fill_manual(values = c("orange1", "#FAE48BFF")) +
  geom_dotplot(binaxis = 'y', 
               stackdir = 'center',
               position = position_dodge(1),
               binwidth = 0,
               dotsize = 10) +
  scale_y_continuous(breaks = seq(0, max(dat[, 2]), 25),
                     labels = scales::number_format(accuracy = 1, 
                                                    decimal.mark = ".")) +
  # scale_x_discrete(breaks = c("GFPT2_High","GFPT2_Low"),
  #                   labels = c(expression(GFPT2^High), expression(GFPT2^Low))) +
  theme_classic() +
  labs(title = NULL, # title = paste0("Correlation between ", colnames(dat)[1], " and ", colnames(dat)[2]), 
       x = NULL, 
       y = yna) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 28), 
        axis.line = element_line(colour = "black", size = 3),
        legend.position = "right",
        legend.key.size = unit(1.5, "cm"),
        legend.text = element_text(size = 20, face = "bold"),
        legend.title = element_text(size = 20, face = "bold"),
        axis.title.y.left = element_text(size = 30, face = "bold", 
                                         margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.title.x.bottom = element_text(size = 30, face = "bold",
                                           margin = margin(t = 0, r = 0, b = 0, l = 0)),
        axis.text = element_text(size = 25, face = "bold", color = "black"),
        axis.ticks.y.left = element_line(colour = "black", size = 2))
p

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/NCI60/Figures")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), 
              colnames(dat)[1],
              "NCI60_Glutathione.Box.plot.tiff", 
              sep = " ")
ggsave(filename = file, units = "in", dpi = 72, width = 10, height = 6)

print(p)
#----------------------------------------------------------------------#

# Sigi's plot, changed to my style

#----------------------------------------------------------------------#
# Boxplot:
my_comparisons <- list(c("GFPT1_High", "GFPT1_Low"))
p2 <- ggplot(total_data3, aes(x = GFPT1_grouped, y = Glutathione, color = GFPT1_grouped)) +
  #geom_point(size = 3) +
  geom_boxplot(size = 1) + 
  stat_compare_means(comparisons = my_comparisons,label = 'p.signif',size=7) +
  scale_x_discrete(
    breaks=c("GFPT1_High","GFPT1_Low"),
    labels=c(expression(GFPT1^High), expression(GFPT1^Low))
  )
p2 <- plot_look_for_paper(p2,y_text = 'Glutathione level',x_text = '',legend_check = F,rotate_check = T) + 
  scale_y_continuous(breaks = seq(0, 175, by=50), limits=c(0,175))
p2 <-p2 + theme(panel.grid.major.y = element_line(color='grey'))
p2

#ggsave('GFPT1_Glutathione_AUGUST2020.png',width = 6, height = 6, units = "in")





############################################################################
# Part 2: Check whether GFPT2 and PTEN are coexpressed in cancer cell lines:
id_PTEN <- 12657
id_GFPT2 <- grep('GFPT2',rownames(expr_dat),)
id_another <- which(rownames(expr_dat) == 'AKT2')

co_mat <- data.frame(t(expr_dat[c(id_PTEN,id_GFPT2,id_another),]))

# Test correlation of PTEN and GFPT2:
cor.test(co_mat[,1],co_mat[,2],method = 'spearman',exact=F)
plot(co_mat[,1],co_mat[,2])

# Test correlation of PTEN and any other random gene:
cor.test(co_mat[,1],co_mat[,3],method = 'spearman',exact=F)
plot(co_mat[,1],co_mat[,3])

# Test correlation of GFPT2 and any other random gene:
cor.test(co_mat[,2],co_mat[,3],method = 'spearman',exact=F)
plot(co_mat[,2],co_mat[,3])


# Check same thing in TCGA data:
vsd_tcga <- readRDS(file = 'TCGA_data_whole_vst_transformed_170519.rds') # NON-z-sscored TCGA data (but normalized)

id_PTEN <- which(rownames(vsd_tcga) == 'GFPT1')
id_GFPT2 <- grep('GFPT2',rownames(vsd_tcga),)
id_another <- which(rownames(vsd_tcga) == 'PTENP1')


co_mat2 <- data.frame(t(vsd_tcga[c(id_PTEN,id_GFPT2,id_another),]))

# Test correlation of PTEN and GFPT2:
cor.test(co_mat2[,1],co_mat2[,2],method = 'spearman',exact=F)
plot(co_mat2[,1],co_mat2[,2])

# Test correlation of PTEN and any other random gene:
cor.test(co_mat2[,1],co_mat2[,3],method = 'spearman',exact=F)
plot(co_mat2[,1],co_mat2[,3])

# Test correlation of GFPT2 and any other random gene:
cor.test(co_mat2[,2],co_mat2[,3],method = 'spearman',exact=F)
plot(co_mat2[,2],co_mat2[,3])


# Try reading in raw TCGA data before they were variance stabilized:
raw_rna <- read.table('data_RNA_Seq_v2_expression_median.txt',header=T)

raw_rna$Hugo_Symbol[grep('PTEN',raw_rna$Hugo_Symbol)]

