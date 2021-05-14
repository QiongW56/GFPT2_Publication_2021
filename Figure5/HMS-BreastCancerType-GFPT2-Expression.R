setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/HMS_Dataset_20348")
library(readxl)
datHMS <- read_excel("13058_2017_855_MOESM2_ESM_SE Smith2017_HMSdataset_final.xlsx", 
                     sheet = "Plotting")
datHMS <- data.frame(datHMS)

# Z-score 
cal_z_score <- function(x){
  (x - mean(x)) / sd(x)
}

datHMS[, ncol(datHMS)] <- cal_z_score(datHMS[, ncol(datHMS)])

# to replace the "Basal B" with "Mesenchymal" for three cell lines
mesCL <- c("BT549","MDAMB231", "MDAMB157", "HS578T")

for (i in 1:nrow(datHMS)){
  if (datHMS$Cell.line[i] %in% mesCL){
    datHMS$Subgroup[i] <- "Claudin-low"
  }
}

datHMS[, 3][grep("Luminal", datHMS[, 3])] <- "Luminal"

library(ggplot2)
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

# start plotting
p <- ggplot(datHMS, aes(x = reorder(Cell.line, GFPT2), y = GFPT2, fill = Subgroup)) + # linetype = 1, 2...
  geom_bar(stat = "identity", position = position_dodge(), color = "black", 
           width = 0.75, size = 4) +
  coord_flip() + 
  scale_fill_manual(values = c("red", "steelblue", "orange1", "#FAE48BFF", "grey100")) + 
  scale_y_continuous(breaks = seq(-2, 7, 1),
                     labels = scales::number_format(accuracy = 1,
                                                    decimal.mark = ".")) +
  theme_classic() +
  labs(title = NULL, 
       y = "RNA_Expression (GFPT2)", 
       x = "Cell_Lines") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 120), 
        axis.line = element_line(colour = "black", size = 4),
        legend.position = "right",
        legend.key.size = unit(2.5, "cm"),
        legend.text = element_text(size = 96, face = "bold"),
        legend.title = element_text(size = 96, face = "bold"),
        axis.title.y.left = element_text(size = 108, face = "bold", 
                                         margin = margin(t = 0, r = 50, b = 0, l = 0)),
        axis.title.x.bottom = element_text(size = 108, face = "bold", 
                                         margin = margin(t = 50, r = 0, b = 0, l = 0)),
        axis.text = element_text(size = 108, face = "bold", color = "black"),
        axis.ticks.y.left = element_line(colour = "black", size = 4))

print(p)

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/HMS_Dataset_20348/Figures")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "HMS-GFPT2-BC", "rotated.bar.plot.tiff", sep = " ")
ggsave(filename = file, units = "in", dpi = 300, width = 35, height = 48)

print(p)
