# This is a modification of code from "KD_GrowthCurve_EMH_24+48+12.R" by using data "KD_GrowthCurve_EMH_Plotting_2020.04.14.xlsx" - Excel sheet: "D492"&"D492M"
# instead of "KD_GrowthCurve_EMH_Plotting.xlsx" - Excel sheet: "D492"&"D492M"

funGC <- function(sheetName, CellType){
  # upload data after counting cells from ImageJ
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/14 Knockdowns Experiments/Results/9 KD_GrowthCurve")
  library(readxl)
  datNo <- read_xlsx("KD_GrowthCurve_EMH_Plotting_2020.04.14.xlsx", sheet = sheetName)
  datNo <- data.frame(datNo)
  
  colnames(datNo)[1] <- "Sample"
  
  # normalize the cell number to fold changes
  dat <- data.frame(datNo[, 1])
  for (i in 2:ncol(datNo)){
    dat[, i] <- datNo[, i]/datNo[, 2]
  }
  
  colnames(dat) <- colnames(datNo)
  
  # add the column for calculating mean and sd later
  dat$Treatment <- factor(c(rep("WT", 12), 
                            rep("Scr", 12),
                            rep("GALE", 12),
                            rep("PGM2L1", 12),
                            rep("UGDH", 12),
                            rep("GFPT2", 12)), 
                          levels = c("WT", "Scr", "GALE", "PGM2L1", "UGDH", "GFPT2"))
  
  dat <- dat[, c(1, 14, 2:13)]
  
  # Calculate mean and sd
  library(tidyverse)
  dat <- dat %>% 
    gather(key = Time.Point, value = Expression, -c(Sample, Treatment)) %>% 
    unite(col = "Group", Time.Point, Treatment) %>% 
    mutate(Group = factor(Group)) %>% 
    group_by(Group) %>% 
    summarize(avg = mean(Expression), sd = sd(Expression))
  
  dat <- data.frame(dat)
  
  # dat <- dat[c(7:84, 1:6), ]
  
  dat$Time.Point <- factor(rep(paste0("T", seq(24, 90, 6), "h"), times = rep(6, 12)), 
                           levels = paste0("T", seq(24, 90, 6), "h"))
  
  dat$Treatment <- factor(rep(c("siGALE", "siGFPT2", "siPGM2L1", "Scr", "siUGDH", "WT"), 12), 
                          levels = c("WT", "Scr", "siGALE", "siPGM2L1", "siUGDH", "siGFPT2"))
  
  dat <- dat[, c(1, 4:5, 2:3)]
  
  # split the dataset into 4 treatments
  datGA <- dat[which(dat$Treatment %in% c("Scr","siGALE")), ]
  datP <- dat[which(dat$Treatment %in% c("Scr", "siPGM2L1")), ]
  datU <- dat[which(dat$Treatment %in% c("Scr", "siUGDH")), ]
  datGF <- dat[which(dat$Treatment %in% c("Scr", "siGFPT2")), ]
  
  funPlot <- function(dat, GN, s, wt, ht) {
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
    
    # color changes depending on the cell type
    if (CellType == "D492"){
      coltype = c("steelblue", "skyblue1")
    } else if (CellType == "D492M") {
      coltype = c("red", "pink")
    } else {
      coltype <- c("orange1", "#FAE48BFF")
    }
    
    # start plotting
    p <- ggplot(dat, aes(x = Time.Point, 
                         y = dat[, 4], 
                         group = Treatment, 
                         color = Treatment)) + # linetype = 1, 2...
      geom_point(size = 20) +
      geom_line(linetype="solid", size = 6) +
      geom_errorbar(aes(ymin = dat[, 4] - sd, ymax = dat[, 4] + sd), 
                    width = 1,
                    size = 3,
                    position = position_dodge(0)) + 
      scale_color_manual(values = coltype) + 
      scale_y_continuous(breaks = seq(0, 3.5, 0.5)) +
      theme_classic() +
      labs(title = paste0("Knockdown of", " ", GN, " ", "in", " ", CellType), 
           x = "Time Point", 
           y = "Fold change in Cell Growth") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 120), 
            axis.line = element_line(colour = "black", size = 4),
            legend.position = "right",
            legend.key.size = unit(2.5, "cm"),
            legend.text = element_text(size = 96, face = "bold"),
            legend.title = element_text(size = 96, face = "bold"),
            axis.title = element_text(size = 96, face = "bold", 
                                             margin = margin(t = 0, r = 50, b = 0, l = 0)),
            axis.text = element_text(size = 84, face = "bold", color = "black"),
            axis.ticks.y.left = element_line(colour = "black", size = 4))
    
    print(p)
    
    setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/14 Knockdowns Experiments/Results/9 KD_GrowthCurve/Figures")
    file <- paste(format(Sys.time(), "%F %H-%M-%S"), GN, CellType, "KD_GC", "line.plot.tiff", sep = " ")
    ggsave(filename = file, units = "in", dpi = s, width = wt, height = ht)
    
    print(p)
  }
  
  funPlot(datGA, "GALE", 72, 42, 24)
  funPlot(datP, "PGM2L1", 72, 42, 24)
  funPlot(datU, "UGDH", 72, 42, 24)
  funPlot(datGF, "GFPT2", 72, 48, 24)
  
}
