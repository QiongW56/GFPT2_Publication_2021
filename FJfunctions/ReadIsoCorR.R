# Something read isocor results

ReadIsoCorR <- function(filename) {

readdata <- read.delim(filename, header = T, sep = "\t")
isodata  <- data.frame(readdata$Isotopic.cluster, readdata$Peak.index, 
                       readdata$Isotopologue.distribution)

nIsot   <- max(isodata$readdata.Peak.index)+1         # number of isotopologues (+1 cause M+0)
sampleN <- readdata$Sample[readdata$Sample != c("")]  # Sample names
nsamples <- length(unique(sampleN))                    # number of samples
nreps    <- length(isodata[,1])/nIsot/nsamples          # number of reps.
namecol  <- NULL

for (i in 1:nsamples )
{
  for (j in 1:nreps)  
  {
    namecol <- c(namecol,rep(as.character(sampleN[i]), nIsot)) 
  }
}



isodata[,length(isodata)+1]     <- namecol
isodata[,length(isodata)]       <- as.factor(isodata[,length(isodata)])
names(isodata)[length(isodata)] <- "Sample" 
isodata$readdata.Peak.index <- as.factor(isodata$readdata.Peak.index)
names(isodata)[1] <- "Isotopic.Cluster"
names(isodata)[2] <- "Isotopologue"
names(isodata)[3] <- "Isotopologue.distribution"

newIsotop <- NULL
for (i in 1:nIsot)
{
  newIsotop[i] <- paste("M+",i-1, sep="")
}

isodata$Isotopologue <- rep(newIsotop, nsamples*nreps)

return(isodata)
}









