#----------------------------------------------------------------------------------------------#

# function 1 for organizing the temp file from TargetLynx

#----------------------------------------------------------------------------------------------#
funSplit <- function(sheetNo, MetNo){ # x is the sheet number and y is the no. of metabolite
  # upload the temp file from TargetLynx with 3 metabolites in each excel sheet
  setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentBasic_03.12.2019/Metabolomics_KD_EMH_Basic_03.12.2019")
  library(readxl)
  rawDat <- data.frame(read_excel("RawData.xlsx", col_names = FALSE, sheet = sheetNo))
  
  # delete the rows with only "NA"
  rawDat <- rawDat[-which(rowSums(is.na(rawDat)) == ncol(rawDat)), ]
  
  # extract all the compounds
  Compd <- gsub("^.*:  ",  
                "", 
                rawDat[which(grepl(":", rawDat[, 1]) & grepl("Compound", rawDat[, 1])), 1])
  # to("IS", , from = 1, same.size = FALSE) https://rdrr.io/cran/lessR/man/to.html
  
  # delete all the necessary rows
  rawDat <- rawDat[-which(rowSums(is.na(rawDat)) > 1), ] # delete all the "comment rows"
  colnames(rawDat) <- rawDat[1, ] # change the dataset colnames
  rawDat <- rawDat[-which(is.na(rawDat[, 1])), ] # delete all the "title rows"
  
  # split rawDat into small datasets based on the metabolites
  rawDat1 <- rawDat[1:25, ][-22:-25, ]
  colnames(rawDat1)[which(colnames(rawDat1) == "Area")] <- Compd[1]
  rawDat2 <- rawDat[26:50, ][-22:-25, ]
  colnames(rawDat2)[which(colnames(rawDat2) == "Area")] <- Compd[2]
  rawDat3 <- rawDat[51:75, ][-22:-25, ]
  colnames(rawDat3)[which(colnames(rawDat3) == "Area")] <- Compd[3]
  
  if (MetNo == 1){
    rawDat <- rawDat1
  } else if (MetNo == 2){
    rawDat <- rawDat2
  } else {
    rawDat <- rawDat3
  }
  
  return(rawDat)
  
  }

#---------------------------------------------------------------------------------#

#---------------------------------------------------------------------------------#
# D492+IS
EdatIS1 <- funSplit(1, 1)
EdatIS2 <- funSplit(1, 2)
EdatIS3 <- funSplit(1, 3)

# D492+Metabolites
EdatMet1 <- funSplit(2, 1)
EdatMet2 <- funSplit(2, 2)
EdatMet3 <- funSplit(2, 3)

# D492M+IS
MdatIS1 <- funSplit(3, 1)
MdatIS2 <- funSplit(3, 2)
MdatIS3 <- funSplit(3, 3)

# D492M+Metabolites
MdatMet1 <- funSplit(4, 1)
MdatMet2 <- funSplit(4, 2)
MdatMet3 <- funSplit(4, 3)

# D492HER2+IS
HdatIS1 <- funSplit(5, 1)
HdatIS2 <- funSplit(5, 2)
HdatIS3 <- funSplit(5, 3)

# D492HER2+Metabolites
HdatMet1 <- funSplit(6, 1)
HdatMet2 <- funSplit(6, 2)
HdatMet3 <- funSplit(6, 3)

#--------------------------------------------------------------------------------------------#

# function 2 to organize the data into the template as in the AcidicNeg run

#--------------------------------------------------------------------------------------------#
# upload the Acidic_Neg data from before with protein concentration information
setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/Glycan precursor analysis/KDexperimentRepeat_16.07.2019/QIONG")
proData <- data.frame(read_excel("AcidicNeg_Knockdown_Metabolomics_3Scr_18.07.2019.xlsx"))
dat <- proData[, 1:2]
dat[, 1] <- gsub("acid", "basic", dat[, 1])

# split into cell types
datE <- dat[1:24, ]
datM <- dat[25:48, ]
datH <- dat[49:72, ]

#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
# Function to add the columns for new metabolites
funPaste <- function(datCell, datMet1, datMet2, datMet3, datMet4, datMet5, datMet6){
  
  funPaste <- function(dat, datMet){
    for (i in 1:nrow(dat)){
      
      dat$met[i] <- datMet[which(grepl(substr(dat$Sample.ID[i], 1, 7), 
                                       datMet[, 4], ignore.case = TRUE)), 5] 
    }
    
    colnames(dat)[ncol(dat)] <- colnames(datMet)[5]
    
    return(dat)
    
  }
  
  for (i in 1:6){
    
    if (i == 1){
      datMet <- datMet1
    } else if (i == 2) {
      datMet <- datMet2
    } else if (i == 3) {
      datMet <- datMet3
    } else if (i == 4){
      datMet <- datMet4
    } else if (i == 5){
      datMet <- datMet5
    } else (
      datMet <- datMet6
    )
    
    datCell <- funPaste(datCell, datMet)
    
    for (i in 2:ncol(datCell)){
      datCell[, i] <- as.numeric(datCell[, i])
    }
    
  }
  
  return(datCell)
  
}

datEfinal <- funPaste(datE, EdatIS1, EdatIS2, EdatIS3, EdatMet1, EdatMet2, EdatMet3)
datMfinal <- funPaste(datM, MdatIS1, MdatIS2, MdatIS3, MdatMet1, MdatMet2, MdatMet3)
datHfinal <- funPaste(datH, HdatIS1, HdatIS2, HdatIS3, HdatMet1, HdatMet2, HdatMet3)

glycanKD <- rbind(datEfinal, datMfinal, datHfinal)
