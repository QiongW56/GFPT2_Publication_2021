# Find mz RT based on name and polarity:

findMZ <- function(name, polarity)
{
  if (missing(polarity))
  {
    polarity = "POS"
  }
  
  mets <- MetTarget(polarity)
  metnumber <- matrix(0,length(mets$mz))
  frag <- matrix(0,4)
  bfrag <- FALSE
  
  for (i in 1:length(mets$mz))
  {
    if (mets$Metabolite[i] == name)
    {
      metnumber[i] <- 1
    }
  }
  
  outp <- which(metnumber != 0)
  
  
  if (length(outp) != 0)
  {
    for (f in 5:8)
    {
      if (mets[outp,f] != 0)
      {
        frag[f] <- 1
        bfrag <- TRUE
      }
    }
    if (!bfrag)
      {
       cat("m/z:", mets$mz[outp], "    RT:",mets$RT[outp],"min.\n")
      }
    else 
    {
      fragnumber <- which(frag != 0)
      cat("m/z:", mets$mz[outp])
      for (j in fragnumber)
      {
        cat(",",mets[outp,j])
      }
      cat("    RT:",mets$RT[outp],"min.\n")
      
    }
    
  }
  
  
  else
  {
    cat("Metabolite not found. \n")
  }
  
}