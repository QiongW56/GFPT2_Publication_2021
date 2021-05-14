# upload the datset from IncuCyte for D492M (Wound Confluency)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/14 Knockdowns Experiments/Results/12 Migration")
library(readxl)
dat <- data.frame(read_excel("GlycanKD_Migration_IncuCyto_ValidImages_2020.01.01.xlsx", 
                             sheet = 2))

# delete unnecessory columns
dat <- dat[, -(which(is.na(dat[1, ])) : ncol(dat))]

# delete unnecessory rows
dat <- dat[-(which(is.na(dat[, 1])) : nrow(dat)), ]

# change the colnames
colnames(dat) <- dat[1, ] 

# delete the first row after changing the colnames
dat <- dat[-1, ]

# delete the first column with date in it
dat <- dat[, -1]

# change from character to numeric
for (i in 2:ncol(dat)){
  dat[, i] <- as.numeric(dat[, i])
}

# delete the unwanted columns based on the tag 
dat <- dat[, c(1, which(as.numeric(dat[1, 2:ncol(dat)]) == 1)+1)]

# delete the tag row (1 or 2)
dat <- dat[-1, ]

# delete the empty rows
# dat <- dat[-which(rowSums(is.na(dat)) == ncol(dat)), ]

WT <- grep("B", names(dat))
Scr <- grep("C", names(dat))
GALE <- grep("D", names(dat))
PGM2L1 <- grep("E[0-9]", names(dat))
UGDH <- grep("F", names(dat))
GFPT2 <- grep("G", names(dat))

# rearrange the columns
dat <- dat[, c(1, 
               WT, # WT B
               Scr, # Scr C
               GALE, # GALE D
               PGM2L1, # PGM2L1 E
               UGDH, # UGDH F
               GFPT2)] # GFPT2 G

# change the column with time points with a "T" in them
dat$Elapsed <- paste0("T", dat$Elapsed, "h")

# change the colnames based on treatments
colnames(dat)[2:ncol(dat)] <- c(paste("WT", 1:length(WT), sep = "_"),
                                paste("Scr", 1:length(Scr), sep = "_"),
                                paste("siGALE", 1:length(GALE), sep = "_"),
                                paste("siPGM2L1", 1:length(PGM2L1), sep = "_"),
                                paste("siUGDH", 1:length(UGDH), sep = "_"),
                                paste("siGFPT2", 1:length(GFPT2), sep = "_"))

# split the datasets
datWT <- dat[, c(1, grep("WT", names(dat)))]
datScr <- dat[, c(1, grep("Scr", names(dat)))]
datGA <- dat[, c(1, grep("GALE", names(dat)))]
datPG <- dat[, c(1, grep("PGM2L1", names(dat)))]
datUG <- dat[, c(1, grep("UGDH", names(dat)))]
datGF <- dat[, c(1, grep("GFPT2", names(dat)))]

# a function to calculate the mean for each dataset
funMean <- function(dat){
  
  ncol <- ncol(dat)
  
  for (i in 1:nrow(dat)){
    dat$Mean[i] <- mean(as.numeric(dat[i, 2:ncol]), na.rm = TRUE)
  }
  
  colnames(dat)[ncol(dat)] <- paste("Mean", 
                                    substr(colnames(dat)[2], 1, (nchar(colnames(dat)[2])-2)),
                                    sep = "_")
    
  return(dat)
  
}

# add a mean column for each dataset
datWT <- funMean(datWT)
datScr <- funMean(datScr)
datGA <- funMean(datGA)
datPG <- funMean(datPG)
datUG <- funMean(datUG)
datGF <- funMean(datGF)

# Merge Scr with KDs
datScrGA <- merge(datScr, datGA)
datScrPG <- merge(datScr, datPG)
datScrUG <- merge(datScr, datUG)
datScrGF <- merge(datScr, datGF)

# function to find the time point with the most difference
funMax <- function(dat){
  
  Avg <- dat[, grep("Mean", colnames(dat))[1]] - dat[, grep("Mean", colnames(dat))[2]]
  
  maxT <-  dat[ ,1][which(Avg == max(Avg))]
  
  maxT <- substr(maxT, 2, nchar(maxT)-1)
  
  maxT <- as.numeric(maxT)
  
  return(maxT)
  
}

maxTGA <- funMax(datScrGA) # T46
maxTPG <- funMax(datScrPG) # T32
maxTUG <- funMax(datScrUG) # T46
maxTGF <- funMax(datScrGF) # T48

funDel <- function(dat, maxT){
  
  dat$No <- substr(dat[, 1], 2, nchar(dat[, 1])-1)
  
  dat$No <- as.numeric(dat$No)
  
  dat <- dat[, c(ncol(dat), 1:(ncol(dat)-1))]
  
  dat <- dat[which(dat[ ,1] <= maxT), ]
  
  dat <- dat[, -grep("Mean", colnames(dat))]
  
  return(dat)
  
}

datScrGA <- funDel(datScrGA, maxTGA)
datScrPG <- funDel(datScrPG, maxTPG)
datScrUG <- funDel(datScrUG, maxTUG)
datScrGF <- funDel(datScrGF, maxTGF)

funSlope <- function(dat){
 
   slope <- NULL
  
   for (i in 3:ncol(dat)){
    temp <- lm(dat[, i] ~ dat[, 1])$coefficients[2]
    slope <- c(slope, temp)
  }
  
    slope <- data.frame(slope)
    
    slope <- cbind(colnames(dat)[3:ncol(dat)], slope)
    
    colnames(slope) <- c("Sample", "Slope")
   
   return(slope)
   
}

datScrGA.S <- funSlope(datScrGA)
datScrPG.S <- funSlope(datScrPG)
datScrUG.S <- funSlope(datScrUG)
datScrGF.S <- funSlope(datScrGF) 
# note: different row numbers in Scramble dataset will change the slope even though same dataset

# function to calculate mean and sd for Scr and KD
funMeanSD <- function(dat){
  
  #dat[, 1] <- as.character.factor(dat[, 1])
  dat$Treatment <- substr(dat[, 1], 1, (nchar(dat[, 1])-2))
  dat <- dat[, c(1, 3, 2)]
  dat$Treatment <- factor(dat$Treatment, levels = unique(dat$Treatment))
  
  library(tidyverse)
  datMsd <- dat %>% 
    group_by(Treatment) %>% 
    summarize(mean = mean(Slope, na.rm = TRUE),
              sd = sd(Slope, na.rm = TRUE))
  
  datMsd <- data.frame(datMsd)
  
  datMsd[ ,1] <- as.character.factor(datMsd[ ,1])
  
  datMsd[1, 1] <- paste0(datMsd[1,1],  "_", datMsd[2,1])
  
  return(datMsd)
  
}

datGAmeanSD <- funMeanSD(datScrGA.S)
datPGmeanSD <- funMeanSD(datScrPG.S)
datUGmeanSD <- funMeanSD(datScrUG.S)
datGFmeanSD <- funMeanSD(datScrGF.S)

# Combined all four knockdowns
dat <- rbind(datGAmeanSD, datPGmeanSD, datUGmeanSD, datGFmeanSD)

#----------------------------------------------------------------------------------------------#

# Delete the specific Gene

#----------------------------------------------------------------------------------------------#
# Delete UGDH
# dat <- dat[-grep("UGDH", dat[ ,1]), ]
