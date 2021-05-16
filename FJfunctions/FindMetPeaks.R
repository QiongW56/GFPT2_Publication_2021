# Automatic peak annotation by using targetlynx method

# Use the mzMin and mzMax from the CAMERA output, add some mz difference to allow for 
# inacuracies in measurements. 
# Works as follows: go through every line of the CAMERA output, until the value in mzMax is 
# larger than mzMax + mz difference.
# If a match is found, check the retention time, if it matches we have a winner.
# after looking at the partent ion, look for fragments if available, do they need to be in the 
# same CAMERA peak group as the parent ion?
# Possible bugs: more than one found, 
# Do we search for fragments if the parent ion is not found?

# INPUT: the CAMERA peaklist dataframe. 
# Columns relevant for search: 1: mz, 2: mzmin, 3: mzmax, 4: rt, last columnn: pcgroup

# OUTPUT values for integrated peaks for each sample, if fragment is there add it to the
# result, do however also output separate values for the parent ion and fragments.
# Columns relevant for output: 7: first sample, number of columns - 3: last sample.



FindMetPeaks <- function(peaklist, metname, nClass ,polarity)
{
  
  if (missing(polarity))
  {
    polarity = "POS"
  }
  
  # First find mz and RT values
  isFrag <- FALSE
  massD  <- 0.005   # This will be customizable
  mzRT   <- AutofindMZ(metname, polarity)
  RT     <- mzRT[length(mzRT) - 1]
  RTdiff <- mzRT[length(mzRT)]
  RTMin <- RT - RTdiff
  RTMax <- RT + RTdiff
  outp  <- NULL
  nFrag <- 0
  
  if (length(mzRT) > 3)
  {
    isFrag <- TRUE
    nFrag <- length(mzRT) - 3
  }

  outp <- matrix(0,nClass,(nFrag+2))
    
for (n in 1:(1+nFrag))  {
  # Do the search...
  parent <- mzRT[n]
  i = 1 # our counter
  while(TRUE) 
  {
    mzMin <- peaklist$mzmin[i] - massD
    mzMax <- peaklist$mzmax[i] + massD
    
    # if it matches we check the RT
    if (parent >= mzMin & parent <= mzMax) 
    {
      # if the RT matches we are done...
      cat("Found something \n")
      cat(peaklist$rt[i]/60, "\n")
      cat(peaklist$mz[i],"\n")
      print(i)
      if (peaklist$rt[i]/60 >= RTMin & peaklist$rt[i]/60 <= RTMax)
        { 
          cat("we found the sucker! \n")
          for (k in 1:nClass){
          outp[k,n] <- peaklist[i,((length(peaklist)-nClass)+k-3)]
          }
          break
        }
    }
    
    # We want to stop if the mz in the ordered peaklist is larger than mz + massD
    if (peaklist$mz[i] > parent+massD)
    {
      break
    }
    
    # if we are at the end of peaklist we want to stop
    if(i == length(peaklist$mz))
    {
      break
    }
    
    # Otherwise we check the next line
    i = i + 1 
  }

}
  
# Get the sum
for (f in 1:(nFrag+1)){
outp[,nFrag+2] <- outp[,nFrag+2] + outp[,f]
}

# Name the columns
outp <- data.frame(outp) 
for (j in 1:(nFrag+1))
{
  colnames(outp)[j] <- mzRT[j]
}
colnames(outp)[length(outp)] <- "Total"  

# Row names:
for (k in 1:nClass){
  row.names(outp)[k] <- colnames(peaklist) [ ( length(peaklist)-nClass )+k-3]
}
  
return(outp)  
}













