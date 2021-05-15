funFrac <- function(No){
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/2 SILIC_DundeePhosphoproteomicsDataset/Results/Peptide table/PerseusOutputTable_SILAC")
  dat <- read.delim(paste0("R", No, ".txt"))
  dat <- data.frame(dat)
  dat <- dat[-1, ]
  
  dat <- dat[, c(24, 1:10)]
  
  for (i in 2:ncol(dat)){
    dat[, i] <- as.numeric(levels(dat[, i]))[dat[, i]]
  }
  
  table(duplicated(dat$Sequence))
  
  rownames(dat) <- dat[, 1]
  dat <- dat[, -1]
  
  library(ComplexHeatmap)
  library(grDevices)
  
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/2 SILIC_DundeePhosphoproteomicsDataset/Results/Peptide table/Heatmap")
  file <- paste(format(Sys.time(), "%F %H-%M-%S"), "Fraction_HeatmapPlot.tiff")
  tiff(filename = file, units="in", width = 7, height = 10, res = 72)
  
  max(dat) # 2.84605
  min(dat) # -2.84605
  col_fun = circlize::colorRamp2(c(-2.85, 0, 2.85), c("steelblue", "white", "red"))
 
  ht <- Heatmap(dat,
                col_fun,
                heatmap_legend_param = list(title = "Expression",
                                            title_gp = gpar(fontsize = 25),
                                            labels_gp = gpar(fontsize = 25),
                                            at = c(--2.85, 0, 2.85),
                                            labels = c("low", "zero", "high"),
                                            legend_height = unit(4, "cm")),
                column_names_side = "top",
                cluster_columns = FALSE,
                show_row_names = FALSE,
                show_row_dend = FALSE,
                column_names_gp = gpar(cex = 2),
                row_names_gp = gpar(cex = 2),
                column_dend_height = unit(3, "cm"),
                km = 10,
                combined_name_fun = NULL, 
                row_title = paste0("Peptites", 1:10, " "),
                row_title_gp = gpar(fontsize = 25),
                column_title = "Fractions",
                column_title_gp = gpar(fontsize = 25))
 
   print(ht)
  
  # export the heatmap plot into .tiff
  dev.off()
  
}