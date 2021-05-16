# loadAll, load all relevant functions from the FJfunctions folder

loadAll <- function() {
  
  
packages = c("ggplot2","scales", "ggrepel","readxl") 
#use this function to check if each package is on the local machine
#if a package is installed, it will be loaded
#if any are not, the missing package(s) will be installed and loaded
  
package.check <- lapply(packages, FUN = function(x) {
if (!require(x, character.only = TRUE)) {
    install.packages(x, dependencies = TRUE)
    library(x, character.only = TRUE)
    }
})
  
files <- list.files("./FJfunctions", pattern = "[.][rR]") 
for (i in files)
{
  cat("Loading:", i, "\n")
  source(file.path("./FJfunctions", i))
}
}