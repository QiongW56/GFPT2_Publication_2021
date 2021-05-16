# MetPlot function, need the classical data frame structure. i.e. need 
# to have timepoints and class columns, with correct names and datatype!

metPlot <- function(dataset, groups, xaxis ,metname, colors, pointsize, isBox){
  
  if (missing(colors))
  {
    colors = c("") # If missing, colors will be an empty set, forcing the algorithm to choose ggplot default
  }
  if (missing(pointsize))
  {
    pointsize = 2.5
  }
  if (missing(isBox))
  {
    isBox = F
  }
  
  if (missing(groups))
  {
    groups = "Class"
  }
 
  
  # find the column corresponding to the selected metabolite, group and xaxis
  colID <- which(colnames(dataset) == metname)
  groupID <- which(colnames(dataset) == groups)
  xaxisIdx <- which(colnames(dataset) == xaxis)
  
  # If group selected is not a factor in the dataset, then change it to a factor:
  if (!is.factor(dataset[,groupID]))
  {
    dataset[,groupID] <- factor(dataset[,groupID])
  }
  
  # If the metabolite is not found:
  if (any(colID) == FALSE ) {
    print("No such metabolite in the data set! Try again.")
    return()
  }
  
  # Check if colors selected match numbers of groups.
  if (length(unique(dataset[,groupID])) < length(colors))
  {
  bColor = TRUE
  }
  else
  {
    bColor = FALSE
    cat("More groups than colors, ggplot default colors will be used.\n")
  }

  
  # Create the plots  
  p <- ggplot(data = dataset, aes(x = dataset[,xaxisIdx], dataset[,colID], colour = dataset[,groupID]))+
    geom_point(size = pointsize) 
  
  # Turn it into boxplot if selected.
  if (isTRUE(isBox)) {
    dataset[,xaxisIdx] <- as.factor(dataset[,xaxisIdx])
    p <- ggplot(data = dataset, aes(x = dataset[,xaxisIdx], dataset[,colID], fill = dataset[,groupID]))+
      geom_boxplot()
  }
  
  # Add colors if selected. 
  if (isTRUE(bColor)) {
    if (isBox)
    {
      p <- p +scale_fill_manual(name = groups, values = alpha(colors, .7) ) +
          theme_bw()
    }
    else
    {
      p <- p  +  scale_color_manual(name = groups, values = colors) +
                  theme_bw()
    }
     
      
      
  }
  else
  {
    p <- p + theme_bw()
  }
  
  
  
  cat("Metabolite plotted:", metname)
  return(p)   
}









