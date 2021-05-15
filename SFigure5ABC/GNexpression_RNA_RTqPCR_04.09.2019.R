funRNAexpr <- function(GN){
  # Upload data for RNA expression of targets
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Figures/SupplementaryFig6")
  library(readxl)
  datRNA <- read_excel("ProteomicPaperTargets_RNAexpression.xlsx", sheet = GN)
  datRNA <- data.frame(datRNA)
  datRNA$"Cell.Line" <- factor(rep(c("D492", "D492M", "D492HER2"), times = c(3, 3, 3)), 
                               levels = c("D492", "D492M", "D492HER2"))
  
  # calculate mean and sd
  funCal <- function(dat){
    # function to calculate mean without taking NA into consideration
    funMean <- function(dat){
      temp <- mean(dat, na.rm = TRUE)
      return(temp)
    }
    
    # function to calculate sd without taking NA into consideration
    funSD <- function(dat){
      temp <- sd(dat, na.rm = TRUE)
      return(temp)
    }
    
    m <- tapply(dat[, 2], dat[, 3], funMean)
    sd <- tapply(dat[, 2], dat[, 3], funSD)
    
    dat <- data.frame(rbind(m, sd))
    
    dat <- data.frame(t(dat))
    dat$"Cell.Line" <- factor(rownames(dat), levels = rownames(dat))
    colnames(dat)[1:2] <- c(GN, "SD") 
    dat <- dat[, c(3, 1:2)]
    
    return(dat)
    
  }
  
  datFil <- funCal(datRNA)
  
#------------------------------------------------------------------------------#
  # plotting
  ggplot.dat <- function(dat, wt, ht, s){
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
    p <- ggplot(dat, aes(x = Cell.Line, y = dat[, 2], fill = Cell.Line)) + # linetype = 1, 2...
      geom_bar(stat = "identity", position = position_dodge(), color = "black", 
               width = 0.75, size = 8) +
      geom_errorbar(aes(ymin = dat[, 2] - SD, ymax = dat[, 2] + SD), 
                    width = 0.2,
                    size = 8,
                    position = position_dodge(0.75)) + 
      scale_fill_manual(values = c("steelblue", "red", "orange1")) + 
      scale_y_continuous(breaks = seq(0, max(dat[, 2]), max(dat[, 2]/5)),
                         labels = scales::number_format(accuracy = 0.01,
                                                        decimal.mark = ".")) +
      theme_classic() +
      labs(title = GN, 
           x = NULL, 
           y = paste0("Fold Changes", " ", "(", GN, ")")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 120), 
            axis.line = element_line(colour = "black", size = 8),
            legend.position = "none",
            legend.key.size = unit(2.5, "cm"),
            legend.text = element_text(size = 96, face = "bold"),
            legend.title = element_text(size = 96, face = "bold"),
            axis.title.y.left = element_text(size = 96, face = "bold", 
                                             margin = margin(t = 0, r = 50, b = 0, l = 0)),
            axis.text = element_text(size = 108, face = "bold", color = "black"),
            axis.ticks.y.left = element_line(colour = "black", size = 8)) +
      geom_hline(yintercept = 1, linetype = "dashed", size = 2)
    
    print(p)
    
    setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Figures/SupplementaryFig6")
    file <- paste(format(Sys.time(), "%F %H-%M-%S"), GN, "RT-qPCR", "bar.plot.tiff", sep = " ")
    ggsave(filename = file, units = "in", dpi = s, width = wt, height = ht)
    
    print(p)
    
  }
  
  ggplot.dat(datFil, wt = 26, ht = 18, s = 72)
  
}

funRNAexpr("GALE")
funRNAexpr("PGM2L1")
funRNAexpr("UGDH")
funRNAexpr("GFPT2")
funRNAexpr("SQRDL")
funRNAexpr("NFE2L2")

