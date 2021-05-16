# Function to rearange data for isocor
# Input: excel sheet, plus number of carbons of the metabolite
# Excel sheet: first column, name of sample, each sample in one row, with 
# increasing M+n, starting with M+0. Each sample contains the same number
# of isotopologues, put zeros where necessary.
IsocorRe <- function(filename, sheetname, metname, numCarb, outpname)
{  
  nargin <- length(as.list(match.call())) -1;
  
  if (nargin < 5)
  {
    print("No output name, set as isocorR.txt")
    outpname <- c("isocorR.txt")
  }
  
  isodata <- read_excel(filename, sheet = sheetname, col_names = FALSE) 
  isodata <- as.data.frame(isodata)
  
  # Find length of isotopologues (len), and samples (nsamp)
  len   <- length(isodata[1,])-1
  nsamp <- length(isodata[,1])
  if (len > numCarb+1)
  {
    print("Either the number of carbons selected is wrong or the input file contains to 
           many rows per sample. \n");
    return()
  }
  
  # To add missing isotopologues if needed: 
  if (len < (numCarb+1))
  {
    zero <- rep(0,nsamp)
    dif <- numCarb + 1 - len
    for (i in 1:dif)
    {
      isodata[len+i+1] <- zero  
    }
  }
  
  # Create new data frame:
  sampcol <- NULL
  metcol  <- NULL
  isocol  <- NULL
  for (i in 1:nsamp) 
  {
    sampcol <- c(sampcol,c(isodata[i,1], rep( c(""), (length(isodata)-2) )))
    metcol  <- c(metcol, c(metname, rep( c(""), (length(isodata)-2) )))
    isocol  <- c(isocol, t(isodata[i, 2:(length(isodata))]))
  }
  nullcol <- rep(c(""), length(isocol))
  
  isocordata <- data.frame(sampcol,metcol, nullcol,isocol)
  write.table(isocordata, file = outpname, sep = "\t", quote = F,row.names = F, col.names = F )
  cat("Data successfully written to:", outpname, "\n")
}
  