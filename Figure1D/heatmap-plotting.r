# Heatmap
library(ComplexHeatmap)
library(grDevices)

setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/EMTMetaboliteMarker_heatmap/Figures")
file <- paste(format(Sys.time(), "%F %H-%M-%S"), "HeatmapPlot.tiff")
tiff(filename = file, units="in", width = 12, height = 12, res = 300)

col_fun_1 = circlize::colorRamp2(c(-1.88, 0, 1.77), c("steelblue", "white", "red"))
ht_1 <- Heatmap(datLFQ, 
                col = col_fun_1,
                heatmap_legend_param = list(title = "",
                                            title_gp = gpar(fontsize = 30),
                                            labels_gp = gpar(fontsize = 30),
                                            at = c(-1.88, 0, 1.77),
                                            labels = c("low", "zero", "high"),
                                            legend_height = unit(5, "cm")),
                column_names_side = "top", 
                row_names_side = "left",
                row_dend_side = "left",
                show_row_dend = FALSE,
                column_names_gp = gpar(cex = 2),
                row_names_gp = gpar(cex = 2),
                row_dend_width = unit(6,"cm"),
                column_dend_height = unit(3, "cm"),
                km = 4,
                # combined_name_fun = NULL,
                row_title = NULL,
                row_title_gp = gpar(fontsize = 30),
                column_title = "LFQ",
                column_title_gp = gpar(fontsize = 30))

col_fun_2 = circlize::colorRamp2(c(-1.15, 0, 1.16), c("steelblue", "white", "red"))
ht_2 <- Heatmap(datSILAC, 
                col = col_fun_2,
                heatmap_legend_param = list(title = "",
                                            title_gp = gpar(fontsize = 30),
                                            labels_gp = gpar(fontsize = 30),
                                            at = c(-1.15, 0, 1.16),
                                            labels = c("low", "zero", "high"),
                                            legend_height = unit(4, "cm")),
                column_names_side = "top", 
                row_names_side = "left",
                row_dend_side = "left",
                show_row_dend = FALSE,
                column_names_gp = gpar(cex = 2),
                row_names_gp = gpar(cex = 2),
                row_dend_width = unit(6,"cm"),
                column_dend_height = unit(3, "cm"),
                km = 4,
               # combined_name_fun = NULL,
                row_title_gp = gpar(fontsize = 30),
                column_title = "SILAC",
                column_title_gp = gpar(fontsize = 30),
                show_heatmap_legend = F,
                show_row_names = F)

ht_list <- ht_1 + ht_2
draw(ht_list, ht_gap = unit(10, "mm"))

dev.off() # export the heatmap plot into .tiff