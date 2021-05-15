enrichment_chart <- function(result_df,
                             top_terms = 15,
                             plot_by_cluster = FALSE,
                             num_bubbles = 4,
                             even_breaks = TRUE) {
  necessary <- c("Term_Description", "Fold_Enrichment", "lowest_p",
                 "Up_regulated", "Down_regulated")
  
  if (!all(necessary %in% colnames(result_df))) {
    stop("The input data frame must have the columns:\n",
         paste(necessary, collapse = ", "))
  }
  
  if (!is.logical(plot_by_cluster)) {
    stop("`plot_by_cluster` must be either TRUE or FALSE")
  }
  
  if (!is.numeric(top_terms) & !is.null(top_terms)) {
    stop("`top_terms` must be either numeric or NULL")
  }
  
  if (!is.null(top_terms)){
    if (top_terms < 1) {
      stop("`top_terms` must be > 1")
    }
  }
  
  # sort by lowest adj.p
  result_df <- result_df[order(result_df$lowest_p), ]
  
  ## Filter for top_terms
  if (!is.null(top_terms)) {
    if (plot_by_cluster & "Cluster" %in% colnames(result_df)) {
      keep_ids <- tapply(result_df$ID, result_df$Cluster, function(x)
        x[seq_len(min(top_terms, length(x)))])
      keep_ids <- unlist(keep_ids)
      result_df <- result_df[result_df$ID %in% keep_ids, ]
    } else if (top_terms < nrow(result_df)){
      result_df <- result_df[seq_len(top_terms), ]
    }
  }
  
  num_genes <- vapply(result_df$Up_regulated,
                      function(x) length(unlist(strsplit(x, ", "))), 1)
  num_genes <- num_genes + vapply(result_df$Down_regulated,
                                  function(x) length(unlist(strsplit(x, ", "))), 1)
  
  result_df$Term_Description <- factor(result_df$Term_Description,
                                       levels = rev(unique(result_df$Term_Description)))
  
  g <- ggplot2::ggplot(result_df, ggplot2::aes_(~Fold_Enrichment, ~Term_Description))
  g <- g + ggplot2::geom_point(ggplot2::aes(color = -log10(result_df$lowest_p),
                                            size = num_genes), na.rm = TRUE)
  g <- g + ggplot2::theme_bw()
  g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, color = "black", face = "bold"),
                          axis.text.y = ggplot2::element_text(size = 12, color = "black", face = "bold"),
                          plot.title = ggplot2::element_blank())
  g <- g + ggplot2::xlab("Fold Enrichment")
  g <- g + ggplot2::theme(axis.title.y = ggplot2::element_blank())
  g <- g + ggplot2::theme(axis.title.x = ggplot2::element_text(color = "black", size = 14, face = "bold"))
  g <- g + ggplot2::labs(size = "# Genes", color = expression(-log[10](p)))
  g <- g + ggplot2::theme(legend.title = ggplot2::element_text(color = "black", size = 14),
                          legend.text = ggplot2::element_text(color = "black", size = 10, face = "bold"))
  
  ## breaks for # genes
  if (max(num_genes) < num_bubbles) {
    g <- g + ggplot2::scale_size_continuous(breaks = seq(0, max(num_genes)))
  } else {
    
    if (even_breaks) {
      brks <- base::seq(0, max(num_genes),
                        round(max(num_genes) / (num_bubbles + 1)))
    } else {
      brks <- base::round(base::seq(0, max(num_genes),
                                    length.out = num_bubbles + 1))
    }
    g <- g + ggplot2::scale_size_continuous(breaks = brks)
  }
  
  g <- g + ggplot2::scale_color_continuous(low = "#f5efef", high = "red")
  
  if (plot_by_cluster & "Cluster" %in% colnames(result_df)) {
    g <- g + ggplot2::facet_grid(result_df$Cluster ~ .,
                                 scales = "free_y", space = "free", drop = TRUE
    )
  } else if (plot_by_cluster) {
    message("For plotting by cluster, there must a column named `Cluster` in the input data frame!")
  }
  
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Annotations/pathfindR/BubbleChart_2020.04.07/Figures")
  file <- paste(format(Sys.time(), "%F %H-%M-%S"), 
                "BubbleChart.tiff", 
                sep = " ")
  
  ggplot2::ggsave(filename = file, units = "in", dpi = 360, width = 7.6, height = 4.8)
  
  return(g)
  
}
