# Upload datasets for D492, D492M and D492HER2
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/14 Knockdowns Experiments/Results/12 Migration/Rfiles")
source("PlottingMigrationRate_D492_2020.01.04.R")
dat[, 1] <- paste0("E_", dat[, 1])
Edat <- dat

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/14 Knockdowns Experiments/Results/12 Migration/Rfiles")
source("PlottingMigrationRate_D492M_2020.01.03.R")
dat[, 1] <- paste0("M_", dat[, 1])
Mdat <- dat

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/14 Knockdowns Experiments/Results/12 Migration/Rfiles")
source("PlottingMigrationRate_HER2_2020.01.04.R")
dat[, 1] <- paste0("H_", dat[, 1])
Hdat <- dat

# combine all datasets for different cell types
dat <- rbind(Edat, Mdat, Hdat)

# add extra column for different cell types
dat$CellType <- c(rep("D492", 8), rep("D492M", 8), rep("D492HER2", 8))
# dat$Gene <- sub(".*si", "", dat[ ,1])

# change column order
dat <- dat[, c(1, 4, 2:3)]

# split datasets for different genes
datGA <- dat[grep("GALE" ,dat[, 1]), ]
datPG <- dat[grep("PGM2L1" ,dat[, 1]), ]
datUG <- dat[grep("UGDH" ,dat[, 1]), ]
datGF <- dat[grep("GFPT2" ,dat[, 1]), ]

# function to factor the column used for coloring later
funFactor <- function(dat){
  
  dat$Treatment <- factor(dat$Treatment, levels = dat$Treatment)
  dat$CellType <- factor(dat$CellType, levels = unique(dat$CellType))
  
  return(dat)
}

datGA <- funFactor(datGA)
datPG <- funFactor(datPG)
datUG <- funFactor(datUG)
datGF <- funFactor(datGF)


# plotting
coltype <- c("steelblue", "skyblue1", "red", "pink", "orange1", "#FAE48BFF")

funPlot <- function(dat, GN){
  # plotting
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
  p <- ggplot(dat, aes(x = CellType, y = dat[, 3], fill = Treatment)) + # linetype = 1, 2...
    geom_bar(stat = "identity", position = position_dodge(), color = "black", 
             width = 0.75, size = 4) +
    geom_errorbar(aes(ymin = dat[, 3] - sd, ymax = dat[, 3] + sd), 
                  width = 0.2,
                  size = 4,
                  position = position_dodge(0.75)) + 
    scale_fill_manual(values = coltype) + 
    scale_y_continuous(breaks = seq(0, 
                                    round(max(dat[, 3])+max(dat[, 4])), 
                                    round(max(dat[, 3])+max(dat[, 4]))/5),
                       labels = scales::number_format(accuracy = 0.01,
                                                      decimal.mark = ".")) +
    theme_classic() +
    labs(title = paste0("Knockdown", " ", "of", " ", GN), 
         x = NULL, 
         y = "Wound Closure Rate") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 120), 
          axis.line = element_line(colour = "black", size = 4),
          legend.position = "none",
          legend.key.size = unit(2.5, "cm"),
          legend.text = element_text(size = 96, face = "bold"),
          legend.title = element_text(size = 96, face = "bold"),
          axis.title.y.left = element_text(size = 108, face = "bold", 
                                           margin = margin(t = 0, r = 50, b = 0, l = 0)),
          axis.text = element_text(size = 108, face = "bold", color = "black"),
          axis.ticks.y.left = element_line(colour = "black", size = 4))
  
  print(p)
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/14 Knockdowns Experiments/Results/12 Migration/Figures")
  file <- paste(format(Sys.time(), "%F %H-%M-%S"), GN, "WoundConfluence", "bar.plot.tiff", sep = " ")
  ggsave(filename = file, units = "in", dpi = 72, width = 26, height = 18)
  
  print(p)
  
}

funPlot(datGA, "GALE")
funPlot(datPG, "PGM2L1")
funPlot(datUG, "UGDH")
funPlot(datGF, "GFPT2")
