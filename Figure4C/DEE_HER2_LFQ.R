funDeeGN <- function(GN){
  #-----------------------------------------------------------------------------------------#
  
  # LFQ
  
  #-----------------------------------------------------------------------------------------#
  # import the LFQ proteomic dataset from Dundee
  setwd("F:/DundeeAnalysisResults/All datasets_organized/Proteo&Phosphoproteo_QIONG_2017-2018")
  library(readxl)
  lfq_dundee <- read_excel("Proteomics_LFQ_summary_16.04.2018.xlsx", sheet = 1)
  lfq_dundee <- data.frame(lfq_dundee)
  lfq_dundee <- lfq_dundee[, c(5, 8:19)]
  lfq_dundee <- lfq_dundee[, c(1:4, 8:13)] # delete D492M
  lfq_dundee <- lfq_dundee[, c(1:4, 8:10, 5:7)]
  
  # Find the expression of a specific gene
  lfq_dundee <- lfq_dundee[which(lfq_dundee[ ,1] == GN), ]
  
  # Delete one row if there are more than one rows
  lfq_dundee <- lfq_dundee[1, ]
  
  # format the dataframe
  lfq_dundee <- data.frame(t(lfq_dundee))
  lfq_dundee <- data.frame(lfq_dundee[-1, ])
  colnames(lfq_dundee) <- GN
  
  lfq_dundee$"Cell.Line" <- factor(rep(c("D492", "D492DEE", "D492HER2"), times = c(3,3,3)), 
                                   levels = c("D492", "D492DEE", "D492HER2"))
  lfq_dundee <- lfq_dundee[, c(2, 1)]
  lfq_dundee[, 2] <- as.numeric(levels(lfq_dundee[, 2]))[lfq_dundee[, 2]]
  lfq_dundee[, 2] <- lfq_dundee[, 2]/1000000
  
  # function to calculate mean and sem
  mean.sd <- function(dat){
    # mean, sd, sem, number
    colna <- colnames(dat)[2]
    
    # since there are "NA" in the dataset, mean() will return "NA"
    funMean <- function(dat){
      
      m <- mean(dat, na.rm = TRUE)
      
      return(m)
    }
    
    funSD <- function(dat){
      
      s <- sd(dat, na.rm = TRUE)
      
      return(s)
    }
    
    m <- tapply(dat[, 2], dat[, 1], funMean)
    
    s <- tapply(dat[, 2], dat[, 1], funSD)
    # n <- tapply(dat[, 2], dat[, 1], length)
    # sem <- s/sqrt(n)
    
    dat <- data.frame(rbind(m, s))
    dat <- data.frame(t(dat))
    dat$"Cell.Line" <- factor(rownames(dat), levels = rownames(dat))
    
    colnames(dat)[1] <- colna
    colnames(dat)[2] <- "SD"
    
    dat <- dat[, c(3, 1:2)]
    
    return(dat)
    
  }
  
  LFQ.plot <- mean.sd(lfq_dundee)
  
  # Plotting
  ggplot.dat <- function(dat){
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
    
    # plotting
    library(ggplot2)
    p <- ggplot(dat, aes(x = Cell.Line, y = dat[, 2], fill = Cell.Line)) + # linetype = 1, 2...
      geom_bar(stat = "identity", position = position_dodge(), width = 0.8, 
               color = "black", size = 8) +
      geom_errorbar(aes(ymin = dat[, 2] - SD, ymax = dat[, 2] + SD), 
                    width = 0.2,
                    size = 8,
                    position = position_dodge(0.8)) + 
      scale_fill_manual(values = c("steelblue", "steelblue2", "orange1")) + 
      scale_y_continuous(breaks = seq(0, max(dat[, 2]), max(dat[, 2]/5)), 
                         labels = scales::number_format(accuracy = 0.01, 
                                                        decimal.mark = ".")) +
      theme_classic() +
      labs(title = GN, 
           x = NULL, 
           y = paste0("Relative Expression", " ", "(", colnames(dat)[2], ")")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 124), 
            axis.line = element_line(colour = 'black', size = 8),
            legend.position = "none",
            legend.key.size = unit(1.3, "cm"),
            legend.text = element_text(size = 60, face = "bold"),
            legend.title = element_text(size = 60, face = "bold"),
            axis.title.y.left = element_text(size = 96, face = "bold",
                                             margin = margin(t = 0, r = 50, b = 0, l = 0)),
            axis.text = element_text(size = 108, face = "bold", color = "black"),
            axis.ticks.y.left = element_line(colour = "black", size = 6))
    
    print(p)
    
    setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Figures/SupplementaryFig7/Figures")
    file <- paste(format(Sys.time(), "%F %H-%M-%S"), colnames(LFQ.plot)[2], "DEE_HER2_LFQ_bar.plot.tiff", sep = " ")
    ggsave(filename = file, units = "in", dpi = 300, width = 30, height = 20)
    
    print(p)
    
  }
  
  ggplot.dat(dat = LFQ.plot)
  
}

# Plotting all the targets
funDeeGN(GN = "GALE")
funDeeGN(GN = "PGM2L1")
funDeeGN(GN = "UGDH")
funDeeGN(GN = "GFPT2")
funDeeGN(GN = "NQO1")

