

metBarO <- function(dataset,group,fillgroup,met,colors, upperError = FALSE, normby = NULL, normsample = NULL, groupby = NULL)
{
  bCol <- TRUE
  if (missing(colors)){
    bCol <- FALSE
  }
  
  
  fillidx <- which(colnames(dataset) == fillgroup)
  groupidx <- which(colnames(dataset) == group)
  metidx   <- which(colnames(dataset) == met)
  nameofmet <- colnames(dataset)[metidx]
  cat("Metabolite: ", nameofmet,"\n")
  groupbyidx   <- which(colnames(dataset) == groupby)
  
  if (!is.null(normby)){
    if (normby != fillgroup & normby != group){
      cat("You must normalized by either input group \n")
      return()
    }
    else{
      normidx <- which(colnames(dataset) == normby)
      if (!(normsample %in% dataset[,normidx])){
        cat("The normsample must be included in the normby group\n")
        return
      }
      normF <- lapply(levels(dataset[,groupbyidx]), function(x){mean(dataset[ dataset[,groupbyidx] == x & dataset[,normidx] == normsample, metidx])})
      print(normF)
      
      for (i in 1:length(normF)){
        if (!is.nan(normF[[i]])){
        dataset[dataset[,groupbyidx] == levels(dataset[,groupbyidx])[i], metidx]  <- dataset[dataset[,groupidx] == levels(dataset[,groupbyidx])[i] ,metidx]/normF[[i]]
        }
      }  
    }
    
  }
  
  # Create the sd and mean columns for the errorbars
  dataset$sd <- rep(0,length(dataset[,1]))
  dataset$mean <- rep(0,length(dataset[,1]))
  for (i in 1:length(dataset[,1]))
  {
    whichgroup <- dataset[i,groupidx]
    whichfill  <- dataset[i,fillidx]
    groupvalues <- subset(dataset,dataset[,groupidx] == whichgroup & dataset[,fillidx] == whichfill)
    dataset$sd[i] <- sd(groupvalues[,metidx]) 
    dataset$mean[i] <- mean(groupvalues[,metidx]) 
  }
  
  if (upperError)
  {
    limits <- aes(ymax = dataset$mean + dataset$sd,
                  ymin = dataset$mean)
  }
  else{
    
    limits <- aes(ymax = dataset$mean + dataset$sd,
                  ymin = dataset$mean - dataset$sd)
  }
  
  # To ensure all groups are factors:
  if (!is.factor(dataset[,groupidx])){
    dataset[,groupidx] <- as.factor(dataset[,groupidx])
  }
  if (!is.factor(dataset[,fillidx])){
    dataset[,fillidx] <- as.factor(dataset[,fillidx])
  }
  
  
  p <- ggplot(data = dataset, aes(x = dataset[,groupidx], y = mean, fill =dataset[,fillidx]))+
    geom_bar(stat="identity", position = position_dodge(0.9), color = "black", size = 1)+
    geom_errorbar(limits, position = position_dodge(0.9),width = 0.25, size = 1) 
  if (bCol){
    p <- p + scale_fill_manual(values = colors)
  }
  p <- p + theme_bw()
  p <- p + ylab("") + xlab("")
  p <- p + ggtitle(nameofmet)
  p <- p + theme(plot.title   = element_text(size=20, hjust= .5, vjust = 2, face = "bold"),
                 legend.title = element_blank(),
                 legend.spacing.y = unit(0, "mm"), 
                 panel.border = element_rect(colour = "white", fill=NA),
                 aspect.ratio = 0.9, axis.text = element_text(colour = 1, size = 12),
                 legend.background = element_blank(),
                 legend.box.background = element_rect(colour = "black"),
                 axis.title.y = element_text(size=16, face = "bold", colour = "black",
                                             margin = margin(t = 0, r = 15, b = 0, l = 0)),
                 axis.text.x  = element_text(size=14, face = "bold", colour = "black"),
                 axis.text.y  = element_text(size=14, face = "bold", colour = "black"),
                 axis.line = element_line(colour ="black", size = 1),
                 axis.ticks = element_line(colour ="black", size = 1),
                 legend.text  = element_text(size =14))
  
  return(p)
}
