#-----------------------------------------------------------------------------------#

# upload the GFPT2 expression data from CCLE website (all cell lines)

#-----------------------------------------------------------------------------------#
# upload the data downloaded from the CCLE website
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/CCLE")
datGFPT2 <- read.delim("mRNA expression (RNAseq)_ GFPT2.txt")

# transpose the data
datGFPT2 <- as.data.frame(t(datGFPT2), stringsAsFactors = FALSE)

# add one more column with cell line names
datGFPT2$Cell.Line <- rownames(datGFPT2)

# change the order of the column
datGFPT2 <- datGFPT2[, c(2, 1)]

# set up the colnames
colnames(datGFPT2) <- datGFPT2[1, ]

# delete the first row
datGFPT2 <- datGFPT2[-1, ]

# select only breast cell lines
datGFPT2_br <- datGFPT2[grep("BREAST", datGFPT2[, 1], ignore.case = TRUE), ]

# delete the NA rows or cell lines
datGFPT2_br <- datGFPT2_br[-which(is.na(datGFPT2_br[, 2])), ]

# character to numeric
datGFPT2_br[, 2]  <- as.numeric(datGFPT2_br[, 2])

# sort the GFPT2 expression
# datGFPT2_br <- datGFPT2_br[order(datGFPT2_br[, 2], decreasing = TRUE), ]

# delete the "_BREAST" part from the cell line name
datGFPT2_br[, 1] <- sub("_BREAST", "", datGFPT2_br[, 1])

#-----------------------------------------------------------------------------------#

# upload the breast cell line meta-data from the paper: Molecular characterization of breast cancer cell lines through multiple omic approaches

#-----------------------------------------------------------------------------------#
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/HMS_Dataset_20348")
library(readxl)
Cell.Char <- read_excel("13058_2017_855_MOESM2_ESM_SE Smith2017.xlsx", sheet = 1)
Cell.Char <- data.frame(Cell.Char)

# keep the subgroup colulm, row 3 has the information
N.col <- grep("subgroup", Cell.Char[3, ], ignore.case = TRUE)

# keep the first and the subgroup column
Cell.Char <- Cell.Char[, c(1, N.col)]

# delete the NA rows
Cell.Char <- Cell.Char[-which(is.na(Cell.Char[, 2])), ]

# delete the "-" in the cell names
colnames(Cell.Char) <- Cell.Char[1, ]
Cell.Char <- Cell.Char[-1, ]
Cell.Char[, 1] <- sub("-", "", Cell.Char[, 1])
Cell.Char[, 1] <- sub("-", "", Cell.Char[, 1])

#-----------------------------------------------------------------------------------#

# merge the two datasets

#-----------------------------------------------------------------------------------#
library(reshape2)
colnames(datGFPT2_br)[1] <- "Cell.Line"
colnames(Cell.Char)[1] <- "Cell.Line"
dat <- merge(datGFPT2_br, Cell.Char)

# delete the Basal A- HER2 cell lines
dat <- dat[-grep("29", dat[, 3]),]

dat[, 3][grep("Luminal", dat[, 3])] <- "Luminal"

# to replace the "Basal B" with "Mesenchymal" for three cell lines
mesCL <- c("BT549","MDAMB231", "MDAMB157")

for (i in 1:nrow(dat)){
  if (dat$Cell.Line[i] %in% mesCL){
    dat$Subgroup[i] <- "Claudin-low"
  }
}

#-----------------------------------------------------------------------------------#

# plotting

#-----------------------------------------------------------------------------------#
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
p <- ggplot(dat, aes(x = reorder(Cell.Line, GFPT2), y = GFPT2, fill = Subgroup)) + # linetype = 1, 2...
  geom_bar(stat = "identity", position = position_dodge(), color = "black", 
           width = 0.75, size = 4) +
  coord_flip() + 
  scale_fill_manual(values = c("red", "steelblue", "orange1", "#FAE48BFF", "grey100")) + # "grey25", "grey50", "grey75"
  scale_y_continuous(breaks = seq(-7, 7, 1),
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
        axis.text = element_text(size = 96, face = "bold", color = "black"),
        axis.ticks.y.left = element_line(colour = "black", size = 4))

print(p)

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/CCLE/Figures")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "CCLE-GFPT2-BC", "rotated.bar.plot.tiff", sep = " ")
ggsave(filename = file, units = "in", dpi = 300, width = 48, height = 48)
