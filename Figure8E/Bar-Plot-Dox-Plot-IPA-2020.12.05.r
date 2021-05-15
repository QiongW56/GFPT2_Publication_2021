# import data
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/ProteomicsManuscript_25.11.2020/Figures/Fig.Signaling/Panel3-IPA")
library(readxl)
dat <- read_excel("pCanonicalPathways_HE.xlsx", sheet = 2)
dat <- data.frame(dat)

# delete the "molecules" column
dat <- dat[, -ncol(dat)]

# rename column names
colnames(dat)[2:3] <- c("P.value", "Z.score")

# reorder the pvalue
dat <- dat[order(dat[, 2], decreasing = FALSE), ]

# factor the pathway column, so later the plotting will be in the order wanted
pathway.level <- dat$Ingenuity.Canonical.Pathways
dat$Ingenuity.Canonical.Pathways <- factor(dat$Ingenuity.Canonical.Pathways, 
                                              levels = pathway.level)

# plotting
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/ProteomicsManuscript_25.11.2020/Figures/Fig.Signaling/Panel3-IPA")
library(ggplot2)
p <- ggplot(data = dat, 
            aes(x = Ingenuity.Canonical.Pathways, 
                y = P.value, 
                group = 1)) +
  geom_bar(stat="identity", 
           width = 0.5, 
           fill = ifelse(dat$Z.score > 0, "orange1", "steelblue"),
           size = 2, 
           alpha = 1) + 
  geom_point(data = dat, 
             aes(x = Ingenuity.Canonical.Pathways,
                 y = abs(Z.score)), 
             size = 3) +
  geom_line(data = dat, 
            aes(x = Ingenuity.Canonical.Pathways,
                y = abs(Z.score)), 
            size = 1) +
  coord_flip() + 
  ggtitle('') +
  theme_classic() +
  ylab("-Log(p.value)") +
  scale_y_continuous(breaks = seq(0, 5, 1),
                     labels = scales::number_format(accuracy = 0.1)) +
  theme(panel.grid = element_blank(),
        legend.position = 'bottom',
        legend.title = element_blank(),
        legend.text = element_text(size = 7),
        legend.key.size = unit(0.5, "lines"),
        plot.title = element_text(hjust = 1, size = 14, face = "bold"),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title.x = element_text(size = 12, face = "bold"),
        axis.title.y = element_text(size = 14, face = "bold"),
        axis.line = element_line(size = 1))

file <- paste(format(Sys.time(), "%F %H-%M-%S"), "BarPlot_IPA_D492HER2_D492.tiff", sep = " ")
ggsave(filename = file, units = "in", dpi = 300, width = 5, height = 5)

