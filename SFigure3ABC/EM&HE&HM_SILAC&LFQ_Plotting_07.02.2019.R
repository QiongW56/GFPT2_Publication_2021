setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/Perseus/EM&HE&HM_06.02.2019_LFQ&SILAC_significant/Perseus_1.6.12.0")
# function to upload data
fun_dat <- function(filename){
  filename <- paste0(filename, ".txt")
  dat <- read.delim(filename)
  dat <- dat[-1, ] # different from output from Perseus 1.6.2.3, that was -1:-2
  dat <- dat[which(dat[, 4] == "YES"),]
  
  for (i in 1:4){
    dat[, i] <- as.character.factor(dat[, i])
  }
  
  for (i in 5:11){
    dat[, i] <- as.numeric(levels(dat[, i]))[dat[, i]]
  }
  
  dat <- dat[, c(2, 9)]
  
  colnames(dat) <- c("GO.annotation", "Enrichment.factor")
  
  return(dat)
}

# The enrichment factor was calculated as (Intersection size/Category size)/(Selection size/Total size) 
# https://archiv.ub.uni-heidelberg.de/volltextserver/22835/1/Dissertation%20Ieva%20Didrihsone%20PDF%20A.pdf

# imput data
# function to reorder the rows
fun_order <- function(dat1){
  dat1 <- dat1[order(dat1$Enrichment.factor), ]
  dat1$GO.annotation <- factor(dat1$GO.annotation, levels = dat1$GO.annotation)
  return(dat1)
}

me_bp <- fun_dat(filename = "PerseusOutPutTable_EM_BP_07.04.2019") # 76
me_bp <- me_bp[which(me_bp$Enrichment.factor > 2), ] # 13
me_bp <- fun_order(me_bp)

me_cc <- fun_dat(filename = "PerseusOutPutTable_EM_CC_07.04.2019") # 20
me_cc <- me_cc[which(me_cc$Enrichment.factor > 1.5), ] # 10
me_cc <- fun_order(me_cc)

he_bp <- fun_dat(filename = "PerseusOutPutTable_HE_BP_07.04.2019") # 67
he_bp <- he_bp[which(he_bp$Enrichment.factor > 2), ] # 17
he_bp <- fun_order(he_bp)

he_cc <- fun_dat(filename = "PerseusOutPutTable_HE_CC_07.04.2019") # 18
he_cc <- he_cc[which(he_cc$Enrichment.factor > 1.25), ] # 15
he_cc <- fun_order(he_cc)

hm_bp <- fun_dat(filename = "PerseusOutPutTable_HM_BP_07.04.2019") # 57
hm_bp <- hm_bp[which(hm_bp$`Enrichment.factor` > 2), ] # 16
hm_bp <- fun_order(hm_bp)

hm_cc <- fun_dat(filename = "PerseusOutPutTable_HM_CC_07.04.2019") # 24
hm_cc <- hm_cc[which(hm_cc$`Enrichment.factor` > 2), ] # 10
hm_cc <- fun_order(hm_cc)

#--------------------------------------------------------------------#
# plotting
library(ggplot2)
library(grid)
library(gridExtra)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/Perseus/EM&HE&HM_06.02.2019_LFQ&SILAC_significant")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "BarPlot_GOAnnotation_Perseus.tiff", sep = " ")
tiff(filename = file, units = "in", res = 300, width = 12.5, height = 12.5)

p1 <- ggplot(data = me_bp, aes(x = GO.annotation, y = Enrichment.factor)) +
  geom_bar(stat="identity", width = 0.75, fill = "orange1") + 
  coord_flip() + 
  ggtitle('D492M vs. D492') +
  theme_classic() +
  xlab("GO Annotation (BP)") +
  scale_y_continuous(breaks = seq(0, 4, 2),
                     labels = scales::number_format(accuracy = 0.1)) +
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "lines"),
        plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 1))

p2 <- ggplot(data = me_cc, aes(x = GO.annotation, y = Enrichment.factor)) +
  geom_bar(stat="identity", width = 0.75, fill = "#FAE48BFF") + 
  coord_flip() + 
  ggtitle('D492M vs. D492') +
  theme_classic() +
  xlab("GO Annotation (CC)") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "lines"),
        plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 1))

p3 <- ggplot(data = he_bp, aes(x = GO.annotation, y = Enrichment.factor)) +
  geom_bar(stat="identity", width = 0.75, fill = "red") + 
  coord_flip() + 
  ggtitle('D492HER2 vs. D492') +
  theme_classic() +
  xlab("GO Annotation (BP)") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "lines"),
        plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 1))

p4 <- ggplot(data = he_cc, aes(x = GO.annotation, y = Enrichment.factor)) +
  geom_bar(stat="identity", width = 0.75, fill = "pink") + 
  coord_flip() + 
  ggtitle('D492HER2 vs. D492') +
  theme_classic() +
  xlab("GO Annotation (CC)") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "lines"),
        plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 1))

p5 <- ggplot(data = hm_bp, aes(x = GO.annotation, y = Enrichment.factor)) +
  geom_bar(stat="identity", width = 0.75, fill = "steelblue") + 
  coord_flip() + 
  ggtitle('D492HER2 vs. D492M') +
  theme_classic() +
  xlab("GO Annotation (BP)") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "lines"),
        plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 1))

p6 <- ggplot(data = hm_cc, aes(x = GO.annotation, y = Enrichment.factor)) +
  geom_bar(stat="identity", width = 0.75, fill = "skyblue1") + 
  coord_flip() + 
  ggtitle('D492HER2 vs. D492M') +
  theme_classic() +
  xlab("GO Annotation (CC)") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "lines"),
        plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 10, face = "bold"),
        axis.title.y = element_text(size = 12, face = "bold"),
        axis.line = element_line(size = 1))

# add "A", "B", "C"... 
p5 <- arrangeGrob(p5, top = textGrob("A.", x = unit(0, "npc")
                                     , y   = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize = 16, fontfamily="Arial", fontface = "bold")))

p6 <- arrangeGrob(p6, top = textGrob("B.", x = unit(0, "npc")
                                     , y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize = 16, fontfamily="Arial", fontface = "bold")))

p3 <- arrangeGrob(p3, top = textGrob("C.", x = unit(0, "npc")
                                     , y  = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize = 16, fontfamily="Arial", fontface = "bold")))

p4 <- arrangeGrob(p4, top = textGrob("D.", x = unit(0, "npc")
                                     , y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize = 16, fontfamily="Arial", fontface = "bold")))

p1 <- arrangeGrob(p1, top = textGrob("E.", x = unit(0, "npc")
                                     , y  = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize = 16, fontfamily="Arial", fontface = "bold")))

p2 <- arrangeGrob(p2, top = textGrob("F.", x = unit(0, "npc")
                                     , y = unit(1, "npc"), just=c("left","top"),
                                     gp=gpar(col="black", fontsize = 16, fontfamily="Arial", fontface = "bold")))

p <- grid.arrange(p5, p6, p3, p4, p1, p2, ncol = 2) # plot HM first, EM last

dev.off()
