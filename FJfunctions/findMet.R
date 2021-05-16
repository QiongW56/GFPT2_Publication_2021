

# Find met based on RT and mz...

findMet <- function(mz,pol, massD)
{

  
  if (missing(pol))
  {
    pol = "POS"
  #  cat("Polarity missing, positive selected as default. \n")
  }
  
  if (missing(massD))
  {
    massD = 0.005
  }
  
  mets <- MetTarget(pol)
  mzMin <- mz - massD
  mzMax <- mz + massD
  
  metnumber <- matrix(0,length(mets$mz))
  fragline <- matrix(0,length(mets$mz))
  fragnumber <- matrix(0,4)

  for (i in 1:length(mets$mz))
  {
    if (mets$mz[i] >= mzMin &  mets$mz[i] <= mzMax)
    {
      metnumber[i] <- 1
    }
    
  }
  outp <- which(metnumber != 0)
  
  if (length(outp) != 0)
  {
    for (j in 1:length(outp))
   {
      cat(mets$Metabolite[outp[j]], "RT:", mets$RT[outp[j]], "min.\n")
    }

  }
  
  for (i in 1:length(mets$mz))
    {
      for (j in 5:8)
      {
        if (mets[i,j] >= mzMin &  mets[i,j] <= mzMax)
        {
          fragline[i] <- 1
        }
      }
    }
    
    i <- which(fragline != 0)
    
    if (length(fragline) != 0 )
    {
      for (j in i)
      {
      cat("Fragment of",mets$Metabolite[j],"RT:", mets$RT[j], "min.\n")
      }
    }
    
    if (length(i) == 0 & length(outp) == 0)
    {
      cat("Metabolite not found. \n")
    }
   
}
  







