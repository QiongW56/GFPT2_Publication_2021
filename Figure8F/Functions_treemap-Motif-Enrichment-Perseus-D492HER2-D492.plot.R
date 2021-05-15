setwd("C:/Users/lenovo/OneDrive - Háskóli Íslands/PC-HI/5 ProteomicPaper/Figures&Tables in the paper/ProteomicsManuscript_25.11.2020/Figures/Fig.Signaling/Panel4-Perseus-MotifEnrichment")
library(readxl)
dat <- read_excel("HE_motif enrichment_all.xlsx")
dat <- data.frame(dat)

# function for treemap
treemap.plot <- function(dat, 
                         col.index = 1, 
                         col.size = 2, 
                         col.color = 3, 
                         color = "Reds",
                         title = "",
                         lower.limit = 0,
                         upper.limit = 0.01,
                         mid.value = 0.005,
                         label.size = 12,
                         lgd = 16,
                         height = 5,
                         width = 6){
  library(treemap)
  library(grDevices)
  # upload data
  dat <- data.frame(dat)
  dat[, col.size] <- as.numeric(dat[, col.size])
  dat[, col.color] <- as.numeric(dat[, col.color])
  
  file <- paste(format(Sys.time(), "%F %H-%M-%S"), "treemap.plot.tiff", " ")
  tiff(filename = file, res = 300, height = height, width = width, units = "in")
  treemap(dat, #Your data frame object
          index = c(colnames(dat)[col.index]),  #A list of your categorical variables
          vSize = colnames(dat)[col.size],  #This is your quantitative variable
          vColor = colnames(dat)[col.color],
          type="value", #Type sets the organization and color scheme of your treemap
          palette = color,  #Select your color palette from the RColorBrewer presets or make your own.
          title = title, #Customize your title
          fontsize.title = 16, #Change the font size of the title
          position.legend = 'bottom',
          mapping = c(upper.limit, mid.value, lower.limit),
          range = c(lower.limit, upper.limit),
          fontsize.labels = label.size,
          fontface.labels = 2,
          border.lwds = 2,
          fontsize.legend = lgd)
  dev.off()
  
  library(treemap)
  # upload data
  dat <- data.frame(dat)
  dat[, col.size] <- as.numeric(dat[, col.size])
  dat[, col.color] <- as.numeric(dat[, col.color])
  
  treemap(dat, #Your data frame object
          index = c(colnames(dat)[col.index]),  #A list of your categorical variables
          vSize = colnames(dat)[col.size],  #This is your quantitative variable
          vColor = colnames(dat)[col.color],
          type="value", #Type sets the organization and color scheme of your treemap
          palette = color,  #Select your color palette from the RColorBrewer presets or make your own.
          title = title, #Customize your title
          fontsize.title = 16, #Change the font size of the title
          position.legend = 'bottom',
          mapping = c(upper.limit, mid.value, lower.limit),
          range = c(lower.limit, upper.limit),
          fontsize.labels = label.size,
          fontface.labels = 2,
          border.lwds = 2,
          fontsize.legend = lgd)
}

treemap.plot(dat = dat)
