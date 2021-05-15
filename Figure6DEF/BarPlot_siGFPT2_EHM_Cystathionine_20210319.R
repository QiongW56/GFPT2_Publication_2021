# Upload the Cystathione data (already normalized to IS and Protein amount)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/ProteomicsManuscript_2021/Cystathionine_pGSK3B")
library(readxl)
dat <- read_excel("ProcessedData_Cystathionine.xlsx")
dat <- data.frame(dat)

# average of the scramble samples
avgE <- mean(dat[grep("E_Scr", dat[, 1]), 2])
avgM <- mean(dat[grep("M_Scr", dat[, 1]), 2])
avgH <- mean(dat[grep("H_Scr", dat[, 1]), 2])

# KD/Scr Ratio
dat[grep("E_", dat[, 1]), 2] <- dat[grep("E_", dat[, 1]), 2]/avgE
dat[grep("M_", dat[, 1]), 2] <- dat[grep("M_", dat[, 1]), 2]/avgM
dat[grep("H_", dat[, 1]), 2] <- dat[grep("H_", dat[, 1]), 2]/avgH

# add treatment column
dat$Treatment <- rep(c(rep("Scramble", 3), rep("KD", 3)), 3)
dat$Treatment <- factor(dat$Treatment, levels = c("Scramble", "KD"))

# add cell type column
dat$CellType <- c(rep("D492", 6), rep("D492M", 6), rep("D492HER2", 6))
dat$CellType <- factor(dat$CellType, levels = unique(dat$CellType))

# combine the (treatment and celltype) two column together
dat$Group <- paste0(dat$Treatment, "_", dat$CellType)
dat$Group <- factor(dat$Group, levels = unique(dat$Group))

# calculate mean and sd
library(tidyverse)
dat <- dat %>%
  group_by(Group) %>% 
  summarize(avg = mean(Cystathionine), sd = sd(Cystathionine))

dat <- data.frame(dat)

dat$CellType <- c(rep("D492", 2), rep("D492M", 2), rep("D492HER2", 2))
dat$CellType <- factor(dat$CellType, levels = unique(dat$CellType))

dat$Treatment <- rep(c("Scramble", "siGFPT2"))
dat$Treatment <- factor(dat$Treatment, levels = unique(dat$Treatment))

dat$Group <- c("D492_Scramble", "D492_KD",
               "D492M_Scramble","D492M_KD",
               "D492HER2_Scramble","D492HER2_KD")
dat$Group <- factor(dat$Group, levels = dat$Group)

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
  labs(title = "Cystathionine", 
       x = NULL, 
       y = "KD/Scr Ratio") +
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

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/ProteomicsManuscript_2021/Cystathionine_pGSK3B/Figures")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "siGFPT2", "Cystathionine", "bar.plot.tiff", sep = " ")
ggsave(filename = file, units = "in", dpi = 72, width = 36, height = 18)

