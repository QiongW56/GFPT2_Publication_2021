funMean <- function(SheetName, n){
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/14 Knockdowns Experiments/Results/15 KDefficiency_siRNA1&2")
  library(readxl)
  dat <- data.frame(read_excel("RNA_GeneExpression.xlsx", sheet = SheetName))
  
  dat$Treatment <- c(rep("Scramble", n), rep("KD", n))
  dat$Treatment <- factor(dat$Treatment, levels = c("Scramble", "KD"))
  
  # calculate mean and sd
  library(tidyverse)
  dat <- dat %>%
    group_by(Treatment) %>% 
    summarize(avg = mean(GFPT2), sd = sd(GFPT2))
  
  dat <- data.frame(dat)
  
  dat$CellType <- sub("_.*", "", SheetName)
  
  dat$Group <- paste0(dat$CellType, "_", dat$Treatment)
  
  return(dat)
  
}

datE <- funMean("D492_GFPT2", 9)
datM <- funMean("D492M_GFPT2", 9)
datH <- funMean("D492HER2_GFPT2_1", 5)

dat <- rbind(datE, datM, datH)

dat$Treatment <- factor(dat$Treatment, levels = unique(dat$Treatment))
dat$CellType <- factor(dat$CellType, levels = unique(dat$CellType))
dat$Group <- factor(dat$Group, levels = unique(dat$Group))

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

# color changes depending on the cell type
coltype <-  c("steelblue", "skyblue1", "red", "pink", "orange1", "#FAE48BFF")

# start plotting
p <- ggplot(dat, 
            aes(x = CellType, 
                y = dat[, 2], 
                fill = Group)) + # linetype = 1, 2...
  geom_bar(stat = "identity", 
           position = position_dodge(), 
           color = "black", 
           width = 0.75, 
           size = 4) +
  geom_errorbar(aes(ymin = dat[, 2] - sd, ymax = dat[, 2] + sd), 
                width = 0.2,
                size = 4,
                position = position_dodge(0.75)) + 
  scale_fill_manual(values = coltype) + 
  scale_y_continuous(breaks = seq(0, round(max(dat[, 2])+max(dat[, 3])), round(max(dat[, 2])+max(dat[, 3]))/5),
                     labels = scales::number_format(accuracy = 0.01,
                                                    decimal.mark = ".")) +
  theme_classic() +
  labs(title = paste0("Knockdown Efficiency Of GFPT2"), 
       x = NULL, 
       y = paste0("KD/Scr Ratio", " ", "(RT-qPCR)")) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 108), 
        axis.line = element_line(colour = "black", size = 4),
        legend.position = "right",
        legend.key.size = unit(2.5, "cm"),
        legend.text = element_text(size = 60, face = "bold"),
        legend.title = element_text(size = 60, face = "bold"),
        axis.title.y.left = element_text(size = 96, face = "bold", 
                                         margin = margin(t = 0, r = 50, b = 0, l = 0)),
        axis.text = element_text(size = 96, face = "bold", color = "black"),
        axis.ticks.y.left = element_line(colour = "black", size = 4)) +
  geom_hline(yintercept = 1, linetype = "dashed", size = 2)

print(p)

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/14 Knockdowns Experiments/Results/15 KDefficiency_siRNA1&2/Figures")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "GFPT2_1", "KD_RT-qPCR", "bar.plot.tiff", sep = " ")
ggsave(filename = file, units = "in", dpi = 72, width = 36, height = 18)

