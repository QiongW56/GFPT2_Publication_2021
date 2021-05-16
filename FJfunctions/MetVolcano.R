# A function that plots a volcano plot 
# Input: dataframe, column indecies, optional: threshold values(p and fold),
# different type of statistical tests? Label of significant values?
# log base?
# ATM: input only numberical data with colum names, samples in rows, features in columns.
# Use the numdata <- alldata[,-c(1,2..,n)] where 1,2..n are non-numerical and/or 
# not to be included in the volcano plot.
# pt = p value threshold, defaul at 0.05
# b = base for log transformation, default 2
# ft = fold change threshold, the default is 2, i.e. twofold change.
# bLable,  a boolean, default FALSE, if selected i.e. bLabel = TRUE, then it displays
# labels on significant points (labels = column names), colors if you want to change the 
# color scheme, order: not significant, larger, smaller.
# For p value adjustment, to account for mutliple testing use the argument padj = "".
# Available methods are = "holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
#   "fdr", "none". See: https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html
# for more information.
# The default is none.

metVolcano <- function(dataset,group1idx, group2idx, pt, ft,b,bLabel,colors, padj, boutp, bori)
{
  
  if(missing(colors))
  {
    colors <- c("grey69", "#BC3C29FF", "#0072B5FF")
  }
  if (missing(pt))
  {
    pt <- 0.05
  }
  if (missing(b))
  {
    b <- 2
  }
  if (missing(bLabel))
  {
    bLabel <- FALSE
  }
  if (missing(padj))
  {
    padj <- c("none")
  }
  if (missing(ft))
  {
    ft <- b
  }
  if (missing(boutp))
  {
    boutp = FALSE
  }
  if (missing(bori))
  {
    bori = FALSE
  }
  
  #if (minvalue)
  #{
#    pcadata[pcadata <= 0] <- NA
#    minValue <- min(pcadata[pcadata >= 0], na.rm = T)
#    pcadata[is.na(pcadata)] <- minValue/2
  #}
  
  
  rowlen <- length(group1idx)+length(group2idx)
  pvalues <- matrix(0,length(dataset[1,]),1)
  foldchange <- matrix(0,length(dataset[1,]),1)
  for (i in 1:length(dataset[1,]))
  {
    loggroup1 <- log(mean(dataset[group1idx,i], rm.na = T),b)
    loggroup2 <- log(mean(dataset[group2idx,i], rm.na = T),b)
    foldchange[i] <- loggroup1-loggroup2
  }


  for (i in 1:length(dataset[1,]))
  {
   data1 <- c(dataset[group1idx,i])
   data2 <- c(dataset[group2idx,i])
   t <- t.test(data1,data2)
   pvalues[i] <- t$p.value 
  }
  

  if (padj != "none")
  {
     pvalues =  p.adjust(pvalues,method = padj, n = length(dataset[1,]))
  }
 
  
  # Create a dataframe used for plotting:
  group <- rep(1,length(dataset[1,]))
  Labels <- matrix("",length(dataset[1,]),1)

  plotdata <- data.frame(Foldchange = foldchange, pvalues = pvalues, group = group, Labels == Labels)
  plotdata$group[plotdata$Foldchange >= log(ft,b) & -log(pvalues,10) >= -log(pt,10) ] <- 2
  plotdata$group[plotdata$Foldchange <= -log(ft,b) & -log(pvalues,10) >= -log(pt,10)] <- 3
  # 
  if (!any(plotdata$group == 2))
  {
    colors <- colors[-2]
  }
 
  alllabels <- colnames(dataset)
  
  for (j in 1:length(dataset[1,]))
  {
    if (plotdata$group[j] != 1)
    {
      plotdata$Labels[j] <- alllabels[j]
    }
    else
    {
      plotdata$Labels[j] <- c("")
    }
  }

  plotdata$group <- factor(plotdata$group)

  maxF <- max(abs(plotdata$Foldchange))
  maxF <- round(maxF)+1
  
 p <- ggplot(data= plotdata, aes(x = Foldchange, y=-log(pvalues,10)))+
    geom_point(aes(fill = group), shape = 21,size = 2.8,
               ,show.legend = FALSE)+
    scale_colour_manual(values=colors)+
    geom_hline(yintercept = -log(pt,10), linetype = 5)+
    geom_vline(xintercept = c(-log(ft,b), log(ft,b)), linetype = 5)+
    scale_x_continuous(breaks = seq(-maxF,maxF, by = 1))+
    scale_y_continuous( breaks = seq(0,10, by = 1)) +
    labs(y= expression(-log[10]~"(p-value)"), x = expression(log[2]~"(fold change)"))+
    theme_bw()+
    theme(axis.text=element_text(size=16),
          axis.title=element_text(size=16,face="bold"))+
   scale_fill_manual(values=colors,name = "")

 if (bLabel)
 {
   p <- p + geom_text_repel(aes(label=plotdata$Labels),
                            segment.color = 'black',
                            segment.size = 0.1)
 }
 
  if (boutp)
  {
    sign <- matrix(0,length(dataset[1,]),1)
    direction <-  matrix("",length(dataset[1,]),1)
    for (i in 1:length(dataset[1,]))
    {
      if (plotdata$pvalues[i] <= pt)
      {
        sign[i] <- 1
      }
      else
      {
        sign[i] <- 0
      }
      if (plotdata$Foldchange[i] >= log(ft,b))
      {
        direction[i] <- "UP"
      }
      else if (plotdata$Foldchange[i] <= -log(ft,b))
      {
        direction[i] <- "DOWN"
      }
      else
      {
        direction[i] <- "No Change"
      }
      
    }
    
    outpdata <- data.frame(Metabolite = alllabels,p_value = plotdata$pvalues, Foldchange = plotdata$Foldchange,
                           Significance = sign, Direction = direction)
    if (bori)
    {
      tempdata <- subset(outpdata, outpdata$Significance == 1 & outpdata$Direction != "No Change")
      outpdata <- dataset[, which(colnames(dataset) %in% tempdata$Metabolite)]
      p <- outpdata
    }
    else
    {
      p <- outpdata
    }
  }
  
  return(p)
}





