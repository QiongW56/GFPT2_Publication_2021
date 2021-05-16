


MetPeakPick <- function(peaklist,nSamples,Polarity) 
{
  outp <- data.frame(row.names = colnames(peaklist)[ (length(peaklist)-3-nSamples+1):(length(peaklist)-3)    ])
  target <- MetTarget(Polarity)
  for (i in 1:length(target$Metabolite))
  {
    met  <- target$Metabolite[i]
    temp <- FindMetPeaks(peaklist,met, nSamples, Polarity)
    outp <- cbind(outp,temp$Total)
    colnames(outp)[i] = met
  }
  
  return(outp)
  
}