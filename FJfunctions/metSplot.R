# MetPlot function, need the classical data frame structure. i.e. need 
# to have timepoints and class columns, with correct names and datatype!

metSplot <- function(dataset, xaxis ,metname, myColors, pointsize, isBox){

  nargin <- length(as.list(match.call())) -1;
  bColor <- TRUE;
  
  if (nargin < 5)
  {
    isBox = F;
  }
  
  if (nargin < 4)
  {
    pointsize = 2.5;
    print("Default point size: 2.5")
  }
  else
  {
      cat("Point size: ",pointsize, "\n")
  }
  
  if (nargin < 3) {
    print("Using ggplot default color scheme. To change, input color vector.")
    bColor <- FALSE;
    myColors <- matrix(0,1,nlevels(dataset$Class));
  }
  
  # find the column corresponding to the selected metabolite
  colID <- which(colnames(dataset) == metname);
  
  # If the metabolite is not found:
  if (any(colID) == FALSE ) {
    print("No such metabolite in the data set! Try again.")
    return()
  }
  
if (missing(xaxis)){
  xaxis <- c("Time")
}
 # if (nlevels(dataset$Class) != length(myColors)) 
#  {
#    print("Number of Colors and Classes don't match, using default ggplot color scheme.")
#    bColor <- FALSE;
#  }
  
  # For the x-axis, can change the tick-interval, max and min values.
  # default, first time point to the last with an interval of 1.
  xaxisIdx <- which(colnames(dataset) == xaxis)

  
  minTime <- min(dataset[,xaxisIdx])
  maxTime <- max(dataset[, xaxisIdx])
  interv <- 1;
  Timeaxis <- c("Day")
  
  timeintcheck <- seq(minTime, maxTime, by = interv)
  if (length(timeintcheck) > 20 )
  {
    interv <- interv*(length(timeintcheck)/14)
    Timeaxis <- c("Time")
  }

  
 p <- ggplot(data = dataset, aes(x = dataset[,xaxisIdx], dataset[,colID], colour = Class))+
      geom_point(size = pointsize)+
   #   geom_smooth(method = "lm")+
      scale_x_continuous(name = Timeaxis, breaks = seq(minTime, maxTime, by = interv))+
      scale_y_continuous(name = c("Abundance"))
 
 if (isTRUE(isBox)) {
   dataset$Time <- as.factor(dataset$Time)
   p <- ggplot(data = dataset, aes(x = Time, dataset[,colID], colour = Class))+
     geom_boxplot()
 }
 
 
 if (isTRUE(bColor)) {
  p <- p  +  scale_color_manual(name = "Class", values = myColors) +
             theme_bw()
 }
 else
 {
   p <- p + theme_bw()
 }



cat("Metabolite plotted:", metname)
return(p)   
}









