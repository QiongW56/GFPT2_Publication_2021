# upload the datset from IncuCyte for D492M (Wound Confluency)
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/4 Experiments documents_QIONG/14 Knockdowns Experiments/Results/12 Migration")
library(readxl)
dat <- data.frame(read_excel("GlycanKD_Migration_IncuCyto_ValidImages_2020.01.01.xlsx", 
                             sheet = 9))
# change the colnames
colnames(dat) <- dat[1, ] 

# delete the first row after changing the colnames
dat <- dat[-1, ]

# delete the first column with date in it
dat <- dat[, -1]

# delete the empty rows
dat <- dat[-which(rowSums(is.na(dat)) == ncol(dat)), ]

# change from character to numeric
for (i in 2:ncol(dat)){
  dat[, i] <- as.numeric(dat[, i])
}

# delete the unwanted columns based on the tag 
dat <- dat[, c(1, which(as.numeric(dat[1, 2:ncol(dat)]) == 1)+1)]

# delete the tag row (1 or 2)
dat <- dat[-1, ]

# rearrange the columns
dat <- dat[, c(1, 
               2:3, 14:15, 26:27, # WT B
               4:5, 16:17, 28:29, # Scr C
               6:7, 18:19, 30:31, # GALE D
               8:9, 20:21, 32:33, 38:39, # PGM2L1 E
               10:11, 22:23, 34:35, 40:41, # UGDH F
               12:13, 24:25, 36:37)] # GFPT2 G

# change the column with time points with a "T" in them
dat$Elapsed <- paste0("T", dat$Elapsed)

# change the colnames based on treatments
colnames(dat)[2:ncol(dat)] <- c(paste("WT", 1:6, sep = "_"),
                                paste("Scr", 1:6, sep = "_"),
                                paste("siGALE", 1:6, sep = "_"),
                                paste("siPGM2L1", 1:8, sep = "_"),
                                paste("siUGDH", 1:8, sep = "_"),
                                paste("siGFPT2", 1:6, sep = "_"))

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
  
  maxT <- substr(maxT, 2, nchar(maxT))
  
  maxT <- as.numeric(maxT)
  
  return(maxT)
  
}

maxTGA <- funMax(datScrGA) # T44
maxTPG <- funMax(datScrPG) # T36
maxTUG <- funMax(datScrUG) # T36
maxTGF <- funMax(datScrGF) # T44

funDel <- function(dat, maxT){
  
  dat$No <- substr(dat[, 1], 2, nchar(dat[, 1]))
  
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
