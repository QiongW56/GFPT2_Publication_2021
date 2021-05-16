# MetPlot function, need the classical data frame structure. i.e. need 
# to have timepoints and class columns, with correct names and datatype!

metPlot <- function(dataset, groups, xaxis ,metname, colors,type, pointsize,linesize){
  
  
  # ADD: point type!
  # ADD: a combined style, that combines average with scatter.
  # ADD: a combined style box and scatter, with transparent boxes.
  # ADD: Error of? Error as standard error? 
  
  if (missing(colors))
  {
    colors = c("") # If missing, colors will be an empty set, forcing the algorithm to choose ggplot default
  }
  if (missing(pointsize))
  {
    pointsize = 2.5
  }
  if (missing(linesize))
  {
    linesize = 1.5
  }
  
  
  if (missing(type))
  {
    type = "scatter"
    # Add message
  }
  # Add a warning here:
  if (missing(groups))
  {
    print("No group selected! Try again.")
    return()
  }
  if (missing(xaxis))
  {
    print("No x-axis variable selected! Try again.")
    return()
  }
  
 
  
  # find the column corresponding to the selected metabolite, group and xaxis
  colID <- which(colnames(dataset) == metname)
  groupID <- which(colnames(dataset) == groups)
  xaxisID <- which(colnames(dataset) == xaxis)
  
  # If the metabolite is not found:
  if (any(colID) == FALSE ) {
    print("No such metabolite in the data set! Try again.")
    return()
  }
  
  # If group selected is not a factor in the dataset, then change it to a factor:
  if (!is.factor(dataset[,groupID]))
  {
    dataset[,groupID] <- factor(dataset[,groupID])
    message(paste(groups,"not a factor, temporarily changed for plotting purposes.\n"))
  }

  
  # Check if colors selected match numbers of groups.
  if (length(unique(dataset[,groupID])) < length(colors))
  {
  bColor = TRUE
  }
  else
  {
    bColor = FALSE
    cat("No colors selected or more groups than colors, ggplot default colors will be used.\n")
  }

  
  # Create the plots 
  # If type is not specified, plot scatter:
  if (type == "scatter")
  {
  p <- ggplot(data = dataset, aes(x = dataset[,xaxisID], dataset[,colID], colour = dataset[,groupID]))+
    geom_point(size = pointsize) 
  }
  
  # Type box:
  else if (type == "box") {
    dataset[,xaxisID] <- as.factor(dataset[,xaxisID])
    p <- ggplot(data = dataset, aes(x = dataset[,xaxisID], dataset[,colID], fill = dataset[,groupID]))+
      geom_boxplot()
  }
  
  # Type average:
  else if (type == "average")
  {
       # Here we simply add the average into each sample, and replot 
       # it n times, n beeing the number in each average category,
       # This is easier than creating a new, smaller dataset.
       # Possible to make faster, as the average calculation needs only
       # to be made once, but I am lazy.
      dataset$average <- rep(0,length(dataset[,1]))
      dataset$stdev <- rep(0,length(dataset[,1]))
      for (i in 1:length(dataset[,1]))
      {
        temp <- subset(dataset[,colID], dataset[,groupID] == dataset[i,groupID] &
                         dataset[,xaxis] == dataset[i,xaxis])
        dataset$average[i] <- mean(temp, na.rm = T)
        dataset$stdev[i] <- sd(temp, na.rm = T)
      }
      
      p <- ggplot(dataset, aes(x=dataset[,xaxisID], y = average, color = dataset[,groupID],
                               shape = dataset[,groupID])) +
        # geom_point(aes(x = Time, y = Glutamine), size = 3, alpha = .2)+
        geom_point( size = pointsize)+
        geom_line(size = linesize)+
        geom_errorbar(aes(ymin=average-stdev, ymax=average+stdev),
                      width=.1, size = linesize*0.8)
        }
    
    
  
  else 
  {
    message(paste(type,"is not a valid option. Scatter plot created.\n"))
    p <- ggplot(data = dataset, aes(x = dataset[,xaxisID], dataset[,colID], colour = dataset[,groupID]))+
      geom_point(size = pointsize) 
  }
  
  
  
  
  
  
  
  # Add colors if selected. 
  if (isTRUE(bColor)) {
    if (type == "box")
    {
      p <- p +scale_fill_manual(name = groups, values = alpha(colors, .7) ) +
          theme_bw()
    }
    else
    {
      p <- p  +  scale_color_manual( values = colors) +
                theme_bw()+
                  theme(legend.title = element_blank())
                  
    }
     
      
      
  }
  else
  {
    p <- p + theme_bw()
  }
  
  
  
  cat("Metabolite plotted:", metname, "\n")
  return(p)   
}









