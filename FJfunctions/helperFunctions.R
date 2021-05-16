# A list of simple usefull functions,

# A helper function for others
excludeAmets <- function(dataset)
{
  allnames  <- colnames(dataset)
  excludeID <- grep("A_", allnames)
  excludeB  <- 1:(length(dataset[1,])) %in% excludeID
  newdataset <- dataset[,!excludeB]
  
}

# A function that gets rid of the average colums 
# in the colnames, as well as time bag ect..
# remember: the time column should be the last one before
# the metabolites.
# this function was created before the helper exlcudeAmets.
# Old Version:
#listMets <- function(dataset)
#{
#  len <- length(dataset[1,])
#  metstart <- which(colnames(dataset) == "Time") + 1;
#  allnames <- colnames(dataset[metstart:len])
  # now find and exlude the average names:
#  excludeID <- grep("A_", allnames) # find the index of average columns
#  excludeB <- 1:(length(allnames)) %in% excludeID # a boolean returns true if value is in excl.
#  outpname <- allnames[!excludeB]  # leave out the exluded
#  return(outpname)
#}

listMets <- function(dataset, start){
  len <- ncol(dataset)
  if (missing(start)){
     metstart <- grep("Time",colnames(extra), ignore.case = T)
     if (any(metstart)){
       metstart <- metstart + 1
     }
     else{
       cat("You need either a column called 'Time' or instert number for first met column \n")
       return()
     }
  }
  else {
    metstart <- start
  }   
  
  outp <- colnames(dataset[metstart:len])
  return(outp)
}






infoMets <- function(dataset)
{
  newdataset <- excludeAmets(dataset)
  rowlen <- length(newdataset[1,])
  collen <- length(newdataset[,1])
  # find start
  metstart <- which(colnames(newdataset) == "Time");
  metnumber <- rowlen-metstart;
  timenumber <- length(unique(newdataset$Time))
  classnumber <- length(unique(newdataset$Class))
  expnumber <- length(unique(newdataset$Experiment))
  

  print(paste("                    Number of metabolites:",metnumber))
  print(paste("           Number of distinct time points:", timenumber))
  print(paste(" Number of different classes (conditions):", classnumber))
  print(paste(" Number of different experiments included:", expnumber))
  
}

# function to load data....
# the number n, is for the first n non-numerical column
loadData <- function(filename, sheet, n)
{
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
  
  outpdata <- read_excel(filename, sheet = sheet) 
  outpdata <- as.data.frame(outpdata)
  
  for (i in 1:n)
  {
    outpdata[,i] <- as.factor(outpdata[,i])
  } 
  
  
  return(outpdata)  
}






