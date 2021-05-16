

# A code that parses the XML output from targetlynx into a nice
# dataframe/excelfile...


ParseTargetXML <- function(filename, outpname, bCSV) {

  
if (missing(bCSV))
{
  bCSV = F
}
  
# This package is needed, add one of those things that checks if it 
# in place.  
packages = c("XML")
  if (!bCSV)
  {
    packages <- c(packages, "xlsx")
  }

  #use this function to check if each package is on the local machine
  #if a package is installed, it will be loaded
  #if any are not, the missing package(s) will be installed and loaded
  package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })


# Read the XML file:
qdata <- xmlParse(filename)  
# Record the root:
top <- xmlRoot(qdata)
# The structure of top of the file is as follows:
# Root - QUANDATASET
# Children: XMLFILE, DATASET, GROUPDATA
# We go into groupdata, there is only one node, GROUP, 
# --- the xmlfile and dataset are not of interest
# In GROUP there are 3 nodes: METHODDATA, SAMPLELISTDATA, CALIBRATIONDATA
# We are only interested in SAMPLELISTDATA, although the calibration data
# might be of interest.
# SAMPLELISTDATA contains the sample nodes, as many as the number of samples, all
# call sample
# Record all the sampledata
sampledata <- top[["GROUPDATA"]][["GROUP"]][["SAMPLELISTDATA"]]

# Each sample node, has children nodes, COMPOUND and USERDATA, the number of 
# compound nodes are as many as the number of compounds.
compounds <- names(sampledata[[1]]) # record the number of compounds, should be the same for all.
outpmatrix <- matrix(0,length(names(sampledata)), length(names(compounds)) - 1)
# For the normalized to IS data, if wanted.
noutpmatrix <- matrix(0,length(names(sampledata)), length(names(compounds)) - 1)
rawfile <- xmlSApply(sampledata, xmlGetAttr, "name") # This gets the raw file name
rawremove <- NULL
isremove <- NULL
samplename <- NULL
Allsamplename  <- xmlSApply(sampledata, xmlGetAttr, "text") # This gets the sample names

# -------------------------------------------------------------------------
# Lets look at each sample at a time, compound by compound.
# -------------------------------------------------------------------------
for (s in 1:length(names(sampledata)))
{
  
  # First get sample name, and raw file name. 
  samplename <- c(samplename, Allsamplename[[s]])
  rawname <- rawfile[s]
 
  
#  compoundname <- sample$METHOD[12]
  # the reason for geting the raw file name is to make sure this is actually a
  # measurement (Targetlynx allows emty lines in the samplelist, this gets recorded in
  # the exported output). 
   if (rawname == "") {
     rawremove <- c(rawremove, s)
     next
   }
   else
   {
     for (j in (1:(length(compounds) - 1))) # the -1 because of the USERDATA note in SAMPLE
          {
            # get the info on sample s, compound c
            temp <- xmlSApply(sampledata[[s]][[j]], xmlAttrs)
            area <- temp$PEAK[[6]]
            outpmatrix[s,j] <- as.numeric(area)
            # Get IS information:
            tempIS <- xmlSApply(sampledata[[s]][[j]][[1]], xmlAttrs)
            isarea <- tempIS[[1]]
            # Now a conditional so we won't try to normalize the IS.
            if (isarea == "")
            {
              if ( !(j %in% isremove))
              {
                isremove <- c(isremove,j)
              }
              
            }
            else
            {
              isarea <- as.numeric(isarea)
              area  <- as.numeric(area)
              noutpmatrix[s,j] <- area/isarea
            }
            
            
     }
   }
}




# Get the compound names:
compnames <- xmlSApply(sampledata[[1]], xmlGetAttr, "name")
# create the column names:
namescol <- c("Sample")
for (i in 1:(length(compnames)-1))
{
  namescol <- c(namescol, compnames[[i]])
}

# Now remove unvanted bullshit
if (!is.null(rawremove))
{
  samplename <- samplename[-rawremove]
  outpmatrix <- outpmatrix[-rawremove,]
}


# create the dataframe:
outpdata <- data.frame(samplename, outpmatrix)
colnames(outpdata) <- namescol



if ( (length(isremove)) != (length(compounds) - 1))
{
  if (!is.null(rawremove)){
  noutpmatrix <- noutpmatrix[-rawremove,]
  }
  noutpdata <- data.frame(samplename, noutpmatrix)
  colnames(noutpdata) <- namescol
  if (!is.null(isremove))
  {
   isremove <- isremove+1
   noutpdata <- noutpdata[,-c(isremove)]
  }
}


  




if (!bCSV)
{
  write.xlsx(outpdata, outpname, sheetName = "Raw data", 
             col.names = TRUE, row.names = FALSE, append = TRUE)
  # If the length of the compounds to remove is equal to the number of compounds
  # then no normalized data sheet is created, other wise it is.
  if ( (length(isremove)) != (length(compounds) - 1))
  {
    write.xlsx(noutpdata, outpname, sheetName = "IS Normalized", 
               col.names = TRUE, row.names = FALSE, append = TRUE)
  }
}
else
{
  write.csv(outpdata, outpname, row.names = FALSE)
  # If the length of the compounds to remove is equal to the number of compounds
  # then no normalized data sheet is created, other wise it is.
  if ( (length(isremove)) != (length(compounds) - 1))
  {
    noutpname <- paste("IS_normalized",outpname, sep = "_")
    write.csv(noutpdata, noutpname, row.names = FALSE)
  }
}

}


