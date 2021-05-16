# read in the target metabolites...


MetTarget <- function(polarity)
{
  
  if (missing(polarity))
  {
    cat("You must select a polarity, positive selected as a default. \n")
    polarity = c("POS")
  }
  
  packages <- c("readxl")
  #if any are not, the missing package(s) will be installed and loaded
  package.check <- lapply(packages, FUN = function(x) 
  {
    if (!require(x, character.only = TRUE)) 
    {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })
    
   
    
    outpdata <- read_excel("TargetMetabolites.xlsx", sheet = polarity) 
    outpdata <- as.data.frame(outpdata)
    
    return(outpdata)  
  
}