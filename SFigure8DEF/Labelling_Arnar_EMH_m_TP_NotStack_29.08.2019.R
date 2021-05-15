# plot different incorporation
iso.m.1 <- function(files = "ARNAR_Targeted_UDP_NaGlu_isotopologues_All_Cell_Lines", 
                 title = "UDP-GlcNAc",
                 TP = "6h",
                 wt = 32,
                 ht = 18,
                 s = 72){
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Metabolomics_Arnar_Labelling_April2019")
  library(readxl)
  dat <- read_xlsx(paste0(files, ".xlsx"), sheet = "No knockdowns")
  dat <- data.frame(dat)
  dat <- dat[, c(1:6, 11)]
  colnames(dat) <- c("Samples", substr(colnames(dat)[2:7], 1, nchar(colnames(dat)[2:7])-4))
  
  # delete the mean enrichment column
  IncoEnrich <- dat[, -7]
  
  # change the name system of the samples in rows
  IncoEnrich[10:27, 1] <- sub("Glc", "Gln", IncoEnrich[10:27, 1], fixed = TRUE) # fix the mistake: should be 1-Gln and 5-Gln NOT Glc
  IncoEnrich[1:27, 1] <- sub("Her2", "D492HER2", IncoEnrich[1:27, 1], fixed = TRUE)
  IncoEnrich[28:57, 1] <- sub("LC", "lc", IncoEnrich[28:57, 1], fixed = TRUE)
  IncoEnrich[28:57, 1] <- sub("LN", "ln", IncoEnrich[28:57, 1], fixed = TRUE)
  IncoEnrich[c(28:30, 43:45), 1] <- sub("-", ",", IncoEnrich[c(28:30, 43:45), 1], fixed = TRUE)
  IncoEnrich[c(28:30, 43:45), 1] <- sub(",2", ",2-", IncoEnrich[c(28:30, 43:45), 1], fixed = TRUE)
  IncoEnrich[, 1] <- c(substr(IncoEnrich[1:27, 1], 1, nchar(IncoEnrich[1:27, 1])-2), IncoEnrich[28:57, 1])
  IncoEnrich[1:27, 1] <- gsub(" ", "_", IncoEnrich[1:27, 1], fixed = TRUE)
  IncoEnrich[ , 1] <- gsub("-", paste0("-", "13C "), IncoEnrich[ , 1], fixed = TRUE)
  IncoEnrich <- IncoEnrich[c(28:42, 43:57, 1:27), ] # change order to D492, D492M and D492HER2

  # change to Percentage format
  percent <- function(x, digits = 2, format = "f", ...) {
    paste0(formatC(100 * x, format = format, digits = digits, ...), "%")
  }
  for (i in 2:6){
    IncoEnrich[, i] <- percent(IncoEnrich[, i]) # one can also download the package "scales" and use the function percent()
  }
  
  # delete the "%" and changed to numeric as percentage is character
  for (i in 2:6){
    IncoEnrich[, i] <- substr(IncoEnrich[, i], 1, nchar(IncoEnrich[, i])-1)
    IncoEnrich[, i] <- as.numeric(IncoEnrich[, i])
  }
  
  # calculate mean and sd for plotting, also added columns for indexing
  lvl1 <- IncoEnrich[, 1][-which(duplicated(IncoEnrich[, 1]))]
  IncoEnrich$"ID" <- factor(IncoEnrich[, 1], levels = lvl1) # add an index column
  
  # move the index column to the second column position
  IncoEnrich <- IncoEnrich[, c(1, 7, 2:6)]
  
  # function for mean and sd of each carbon incorporation
  funCal <- function(dat, eqn, nas){
    ms <- data.frame(matrix(0, 19, 5))
    ms <- cbind(lvl1, ms)
    for (i in 3:7){
      ms[, i-1] <- data.frame(tapply(dat[, i], dat[, 2], eqn))
    }
    colnames(ms) <- c("Samples", paste0(nas, "_", "m", 0:4))
    
    return(ms)
  }
  
  mes <- funCal(dat = IncoEnrich, eqn = mean, nas = "Mean") # Means
  stv <- funCal(dat = IncoEnrich, eqn = sd, nas = "SD") # SD
  
  # function to swap D492HER2_1,2Glc_6h and D492HER2_1Gln_6h
  funChg <- function(dat){
    dat$Samples <- as.character.factor(dat$Samples)
    
    rowna13 <- dat$Samples[13]
    rowna16 <- dat$Samples[16]
    
    dat$Samples[13] <- "a"
    dat$Samples[16] <- "b"
    
    dat$Samples[13] <- rowna16
    dat$Samples[16] <- rowna13
    
    dat <- dat[c(1:12, 16, 14:15, 13, 17:19), ]
    
    return(dat)
  }
  
  mes <- funChg(mes)
  stv <- funChg(stv)
  
  # delete unlabelled samples
  mes <- mes[-c(4:5, 9:10, 12, 15, 18), ]
  stv <- stv[-c(4:5, 9:10, 12, 15, 18), ]
  
  # change numeric to character so the function "melt" will be working
  for (i in 2:6){
    mes[, i] <- as.character(mes[, i])
    stv[, i] <- as.character(stv[, i])
  }
  
  # melt the datasets "mes" and "stv"
  library(reshape2)
  mes.mlt <- melt(mes, id = "Samples")
  mes.mlt[, 3] <- as.numeric(mes.mlt[, 3])
  colnames(mes.mlt)[3] <- "Percentage"
  
  stv.mlt <- melt(stv, id = "Samples")
  stv.mlt[, 3] <- as.numeric(stv.mlt[, 3])
  colnames(stv.mlt)[3] <- "SD"
  
  # Merge (mean melt) and (stv melt) together
  datSum <- cbind(mes.mlt, stv.mlt[, 3])
  colnames(datSum)[2] <- "Incorporation"
  colnames(datSum)[4] <- "SD"
  datSum[, 2] <- substr(datSum[, 2], 6, 7)
  datSum$"Cell.Line" <- rep(c(rep("D492", 3), rep("D492M", 3), rep("D492HER2", 6)), 5)
  datSum <- datSum[, c(1, 5, 2:4)]
  
  # prepare data for plotting
  datSum_6h <- datSum[c(1:6, 8, 10, 12, 
                        13:18, 20, 22, 24, 
                        25:30, 32, 34, 36, 
                        37:42, 44, 46, 48,
                        49:54, 56, 58, 60), ]
  datSum_24h <- datSum[c(1:6, 7, 9, 11, 
                         13:18, 19, 21, 23, 
                         25:30, 31, 33, 35, 
                         37:42, 43, 45, 47,
                         49:54, 55, 57, 59), ]
  
  TP <- tolower(TP)
  if (TP == "6h"){
    mPlot <- datSum_6h
    mPlot[, 1] <- substr(mPlot[, 1], 1, nchar(mPlot[, 1])-3)
    mPlot$"Treatment" <- rep(rep(c("1,2-13C Glc", "1-13C Gln", "5-13C Gln"), 3), 5) 
    mPlot <- mPlot[, c(1:2, 6, 3:5)]
    
  } else if (TP == "24h") {
    mPlot <- datSum_24h
    mPlot[, 1] <- gsub("\\_6.*","", mPlot[, 1])
    mPlot[, 1] <- gsub("\\_24.*","", mPlot[, 1])
    mPlot$"Treatment" <- rep(c("1,2-13C Glc", "1-13C Gln", "5-13C Gln"), 5) 
    mPlot <- mPlot[, c(1:2, 6, 3:5)]
    
  } else {
    cat("please use \"6h\" or \"24h\".")
  }
  
  # separate the three different 13C incorporation
  mPlot_1 <- mPlot[which(mPlot$Treatment == "1,2-13C Glc"), ] # 1,2-13C Glc
  mPlot_2 <- mPlot[which(mPlot$Treatment == "1-13C Gln"), ] # 1-13C Gln
  mPlot_3 <- mPlot[which(mPlot$Treatment == "5-13C Gln"), ] # 5-13C Gln
  
  # function to change the columns into factor for plotting
  funFac <- function(dat){
    dat[, 2] <- factor(dat[, 2], levels = c("D492", "D492M", "D492HER2"))
    dat[, 4] <- factor(dat[, 4], levels = c("m0", "m1", "m2", "m3", "m4"))
    
    return(dat)
  }
  
  mPlot_1 <- funFac(mPlot_1) # 1,2-13C Glc
  mPlot_2 <- funFac(mPlot_2) # 1-13C Gln
  mPlot_3 <- funFac(mPlot_3) # 5-13C Gln
  
  # plotting
  funPlot <- function(datPlot, labeling){
    library(ggplot2)
    library("RColorBrewer")
    
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
    
    # start to plot
    p <- ggplot(datPlot, aes(x = Incorporation, y = datPlot[, 5], fill = Cell.Line)) + # linetype = 1, 2...
      geom_bar(stat = "identity", position = position_dodge(), width = 0.8, 
               color = "black", size = 4) +
      geom_errorbar(aes(ymin = datPlot[, 5] - SD, ymax = datPlot[, 5] + SD), 
                    width = 0.2,
                    size = 4,
                    position = position_dodge(0.8)) + 
      scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
      scale_fill_manual(values = c("steelblue", "red", "orange1")) + 
      theme_classic() +
      labs(title = title, 
           x = labeling, 
           y = "13C incorporation (%)") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 120), 
            axis.line = element_line(colour = 'black', size = 4),
            legend.position = "right",
            legend.key.size = unit(2, "cm"),
            legend.text = element_text(size = 64, face = "bold"),
            legend.title = element_text(size = 64, face = "bold"),
            axis.title.y.left = element_text(size = 108, face = "bold",
                                             margin = margin(t = 0, r = 50, b = 0, l = 0)),
            axis.title.x.bottom = element_text(size = 108, face = "bold",
                                               margin = margin(t = 50, r = 0, b = 0, l = 0)),
            axis.text = element_text(size = 96, face = "bold", color = "black"),
            axis.ticks.y.left = element_line(colour = "black", size = 4))
    
    setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Metabolomics_Arnar_Labelling_April2019/Figures")
    file <- paste(format(Sys.time(), "%F %H-%M-%S"), title, labeling, "bar.plot.tiff", sep = " ")
    ggsave(filename = file, units = "in", dpi = s, width = wt, height = ht)
    
    print(p)
  }
  
funPlot(mPlot_1, labeling = mPlot_1$Treatment[1])
funPlot(mPlot_2, labeling = mPlot_2$Treatment[1])
funPlot(mPlot_3, labeling = mPlot_3$Treatment[1])

}

#------------------------------------------------------------------------#

                 # plot for UDP-GlcNAc and UDP-Glc

#------------------------------------------------------------------------#
iso.m.1()
iso.m.1(files = "ARNAR_Targeted_UDP_Glu_isotopologues_All_Cell_Lines",
        title = "UDP-Glc")
# brewer.pal(n = 5, name = cols)
