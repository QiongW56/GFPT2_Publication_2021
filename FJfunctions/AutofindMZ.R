# Find mz RT based on name and polarity:

AutofindMZ <- function(name, polarity)
{
  
  mets <- MetTarget(polarity)
  metnumber <- matrix(0,length(mets$mz))
  frag <- matrix(0,4)
  bfrag <- FALSE
  outp <- NULL
  
 for (i in 1:length(mets$mz))
  {
    if (mets$Metabolite[i] == name)
    {
      metnumber[i] <- 1
    }
  }
  
 parent <- which(metnumber != 0)
  
 if (length(parent) != 0)
   {
   
    for (f in 5:8)
    {
      if (mets[parent,f] != 0)
      {
        frag[f] <- 1
        bfrag <- TRUE
      }
    }
    
  # If there are no fragments:  
  if (!bfrag)
    {
      # the output is [mz, RT, RTdiff]
       outp <- c(mets$mz[parent],mets$RT[parent], as.double(mets$RTdiff[parent]))
    }
 else 
    {
      fragnumber <- which(frag != 0)
      fragout <- matrix(0,length(fragnumber))
      count <- 1
      for (j in fragnumber)
      {
        fragout[count] <- mets[parent,j]
        count = count + 1
      }
      
      # the output [parent, frag 1,..., RT]
      outp <- c(mets$mz[parent], fragout,mets$RT[parent], as.double(mets$RTdiff[parent]))
      
      
    }
    
}
  return(outp)

}




