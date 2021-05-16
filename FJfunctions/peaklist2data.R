# CAMERA peaklist helper functions

# Turn a camera peaklist into a standard dataframe
peaklist2data <- function(peaklist,xcmsset)
{
  
  samplenames  <- row.names(xcmsset@phenoData)
  samplenr <- length(xcmsset@phenoData[,1])
  SampleClass  <- xcmsset@phenoData[,1]
  classnr <- length(unique(SampleClass))
  
  # str of peaklist:
  # 1: mz,  2: mzmin, 3: mzmax, 4: rt, 5: rtmin, 6: rtmax, 7: npeaks, 8:(8+classnr-1): sampleclass,
  # (8+classnr):(8+classnr+samplenr-1): samples
  # Plan: retreive mz, sampleclass, samples, transpose...
  
  newdata <- data.frame(Sample = samplenames, Class = SampleClass, Time = c(1:samplenr) )
  for (i in 1:length(peaklist[,1]))
  {
    cname <- as.character( round(peaklist[i,1],5))
    newdata <- cbind(newdata, t(peaklist[i, (8+classnr):(8+classnr+samplenr-1)]))
    colnames(newdata)[3+i] = cname
  }
  
  return(newdata)
  
}

peakRT <- function(peaklist, mz)
{
  e <- 0.0001
  id <- which(peaklist$mz >= mz-e & peaklist$mz <= mz+e)
  RT = peaklist$rt[id]
  group = peaklist$pcgroup[id]
  cat("RT: ", RT, "\n")
  cat("PCgroup: ", group, "\n")
  return(id)
}







