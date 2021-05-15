#------------------------------------------------------------------------------------------#

                                                # LFQ

#------------------------------------------------------------------------------------------#
# Upload LFQ peptide data
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/2 SILIC_DundeePhosphoproteomicsDataset/Results/Peptide table")
library(readxl)
pepsLFQ <- read_excel("Peptide coverage_LFQ VS SILAC.xlsx", sheet = "LFQ")
pepsLFQ <- data.frame(pepsLFQ)
colnames(pepsLFQ)[2:3] <- c("Start.Position", "End.Position")

# Plotting
library(ggplot2)
p <- ggplot(pepsLFQ) +
  geom_point(aes(pepsLFQ[, 2], pepsLFQ[, 3]), size = 15, color = "steelblue") + # alpha = 0.5
  theme_classic() +
  xlab("Start.Position") +
  ylab("End.Position") +
  scale_x_continuous(limits = c(0, 25000), 
                     labels = scales::number_format(accuracy = 1, big.mark = ",")) + # max(pepsLFQ$Start.Position, na.rm = TRUE)
  scale_y_continuous(limits = c(0, 25000),
                     labels = scales::number_format(accuracy = 1, big.mark = ",")) + # max(pepsLFQ$End.Position, na.rm = TRUE)
  theme(line = element_line(size = 4),
        axis.title.x.bottom = element_text(size = 72, face = "bold",
                                           margin = margin(t = 25 , r = 0, b = 0, l = 0)),
        axis.title.y.left = element_text(size = 72, face = "bold",
                                         margin = margin(t = 0 , r = 25, b = 0, l = 0)),
        axis.text = element_text(size = 64, face = "bold", color = "black"),
        axis.ticks = element_line(size = 8, color = "black")) 

print(p)

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/2 SILIC_DundeePhosphoproteomicsDataset/Results/Peptide table")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "LFQ", "Peptides.tiff", sep = " ")
ggsave(filename = file, units = "in", dpi = 72, width = 32, height = 18)

#------------------------------------------------------------------------------------------#

                                            # SILAC

#------------------------------------------------------------------------------------------#
# Upload SILAC (all three replicates) peptides data
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/2 SILIC_DundeePhosphoproteomicsDataset/Results/Peptide table")
library(readxl)
pepsSILAC <- read_excel("Peptide coverage_LFQ VS SILAC.xlsx", sheet = "SILAC")
pepsSILAC <- data.frame(pepsSILAC)
pepsSILAC <- pepsSILAC[, -c(1, 4, 5, 8, 9)]

T1 <- pepsSILAC[1:nrow(pepsSILAC), 1:2]
T2 <- pepsSILAC[1:nrow(pepsSILAC), 3:4]
T3 <- pepsSILAC[1:nrow(pepsSILAC), 5:6]

colna <- c("Start.Position", "End.Position")
colnames(T1) <- colna
colnames(T2) <- colna
colnames(T3) <- colna

pepsSILAC <- rbind(T1, T2, T3)
pepsSILAC$"Replicates" <- factor(c(rep("Replicate 1", nrow(T1)),
                                   rep("Replicate 2", nrow(T2)),
                                   rep("Replicate 3", nrow(T3))),
                                 levels = c("Replicate 1", "Replicate 2", "Replicate 3"))



# Plotting
library(ggplot2)
p <- ggplot(pepsSILAC) +
  geom_point(aes(Start.Position, End.Position, color = Replicates), size = 15) + # alpha = 0.5
  scale_color_manual(values = c("steelblue", "red", "orange1")) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 40000),
                     labels = scales::number_format(accuracy = 1, 
                                                    big.mark = ",")) + # max(pepsSILAC$Start.Position, na.rm = TRUE)
  scale_y_continuous(limits = c(0, 40000),
                     labels = scales::number_format(accuracy =1, 
                                                    big.mark = ",")) + # max(pepsSILAC$End.Position, na.rm = TRUE)
  theme(line = element_line(size = 4),
        axis.title.x.bottom = element_text(size = 72, face = "bold",
                                           margin = margin(t = 25 , r = 0, b = 0, l = 0)),
        axis.title.y.left = element_text(size = 72, face = "bold",
                                         margin = margin(t = 0 , r = 25, b = 0, l = 0)),
        axis.text = element_text(size = 64, face = "bold", color = "black"),
        axis.ticks = element_line(size = 8, color = "black"),
        legend.text = element_text(size = 64, face = "bold"),
        legend.title = element_text(size = 64, face = "bold")) 

print(p)

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/2 SILIC_DundeePhosphoproteomicsDataset/Results/Peptide table")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "SILAC", "Peptides.tiff", sep = " ")
ggsave(filename = file, units = "in", dpi = 72, width = 32, height = 18)
