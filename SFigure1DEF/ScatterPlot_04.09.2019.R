# import proteins found in both LFQ and SILAC
#----------------------------------------------------------------------------------------------#

                                       # load EM data

#----------------------------------------------------------------------------------------------#
library(readxl)
library(ggplot2)
library(ggrepel)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/ScatterPlot")
datEM <- read_excel("ScatterPlot.xlsx", sheet = "D492 VS D492M") # 2059
datEM <- data.frame(datEM)
rownames(datEM) <- datEM$Gene.Name
datEM <- datEM[, -1:-2]
colnames(datEM) <- c("GeneName", "x", "y", "Grouping")
datEM$Grouping <- as.factor(datEM$Grouping)

# to find the gene names appear twice or more
GN <- datEM$GeneName
names(table(GN))[table(GN) > 1] # 0

#----------------------------------------------------------------------------------------------#

                                           # load HE data

#----------------------------------------------------------------------------------------------#
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/ScatterPlot")
datHE <- read_excel("ScatterPlot.xlsx", sheet = "D492HER2 VS D492")
datHE <- data.frame(datHE)
datHE <- datHE[ ,-6:-7]
rownames(datHE) <- datHE$Gene.Name
datHE <- datHE[, -1:-2]
colnames(datHE) <- c("GeneName", "x", "y", "Grouping")
datHE$Grouping <- as.factor(datHE$Grouping)

#----------------------------------------------------------------------------------------------#

                                        # load HM data

#----------------------------------------------------------------------------------------------#
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/ScatterPlot")
datHM <- read_excel("ScatterPlot.xlsx", sheet = "D492HER2 VS D492M")
datHM <- data.frame(datHM)
datHM <- datHM[ ,-6:-7]
rownames(datHM) <- datHM$Gene.Name
datHM <- datHM[, -1:-2]
colnames(datHM) <- c("GeneName", "x", "y", "Grouping")
datHM$Grouping <- as.factor(datHM$Grouping)

#-------------------------------------------------------------------------------------------#

                                        # Plotting

#-------------------------------------------------------------------------------------------#
funScat <- function(dat){
  # for distinguishing the name of the plots
  if (sum(dat$x) == sum(datEM$x)) {
    na <- "EM"
    xtitle <- "(D492/D492M)_LFQ"
    ytitle <-"(D492/D492M)_SILAC"
    cols <- "orange1"
    xequ <- -3
    yequ <- 6
    
  } else if (sum(dat$x) == sum(datHE$x)){
    na <- "HE"
    xtitle <- "(D492/D492HER2)_LFQ"
    ytitle <-"(D492/D492HER2)_SILAC"
    cols <- "red"
    xequ <- -3
    yequ <- 6
    
  } else {
    na <- "HM"
    xtitle <- "(D492M/D492HER2)_LFQ"
    ytitle <-"(D492M/D492HER2)_SILAC"
    cols <- "steelblue"
    xequ <- -4.2
    yequ <- 3
    
  }
  
  # ploting with pearson cor
  cor <- cor.test(dat$x, dat$y, 
                  method = "pearson", conf.level = 0.95)
  cor1 <- as.numeric(format(cor[4], digits = 3))
  cor2 <- paste("Pearson Correlation =", cor1, sep = " ")
  
  xmin <- min(dat[, 2])
  xmax <- max(dat[, 2])
  xtk <- round((xmax - xmin)/5)
  
  
  ymin <- min(dat[, 3])
  ymax <- max(dat[, 3])
  ytk <- round((ymax - ymin)/5)
  
  p <- ggplot(dat, aes(x = x, y = y)) +
    geom_point(color = cols, size = 6) +
    theme_classic() +
    scale_x_continuous(breaks = seq(xmin, xmax, xtk), 
                       labels = scales::number_format(accuracy = 0.01, 
                                                      decimal.mark = ".")) +
    scale_y_continuous(breaks = seq(ymin, ymax, ytk), 
                       labels = scales::number_format(accuracy = 0.01, 
                                                      decimal.mark = ".")) +
    labs(title = "", # LFQ vs SILAC (D492 vs. D492M)
         x = bquote("Log"[2]*.(xtitle)),
         y = bquote("Log"[2]*.(ytitle))) +
    geom_smooth(method = lm, color = "black", size = 3) +
    theme(legend.position = "none",
          axis.line = element_line(colour = 'black', size = 2),
          axis.title.y.left = element_text(size = 24, face = "bold",
                                           margin = margin(t = 0, r = 10, b = 0, l = 0)),
          axis.text = element_text(size = 32, face = "bold", color = "black"),
          axis.title.x.bottom = element_text(size = 32, face = "bold",
                                    margin = margin(t = 20, r = 0, b = 0, l = 0)),
          axis.ticks.y.left = element_line(colour = "black", size = 2)) +
    geom_text(x = xequ, y = yequ, label = cor2, color = "black", fontface = "bold", size = 8)
  # geom_text_repel(aes(x = x, y = y, 
  # label = ifelse(Grouping == "YES", GeneName, "")), size = 3)
  p
  
  # save the plot
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/ScatterPlot")
  file <- paste(format(Sys.time(), "%F %H-%M-%S"), "_", na, "_","box.plot.tiff", sep = " ")
  ggsave(filename = file, units = "in", dpi = 72, width = 10, height = 6)
  
}

funScat(datEM)
funScat(datHE)
funScat(datHM)
