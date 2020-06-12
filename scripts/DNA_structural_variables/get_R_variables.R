# install packages
if (!requireNamespace("DNAshapeR", quietly = TRUE)){
  print("Installing R modules - reminder to set Rscript path in main Matlab script.")
  if (!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
  }
  BiocManager::install("DNAshapeR", version = "3.8")
}

# run
library(DNAshapeR)
filename <- 'R_tmp.fa'
data <- getShape(filename, shapeType = 'Default', parse = TRUE)

# process data and save
data[5] <- NULL
x <- data.frame(data)
x_omit <- t(na.omit(as.data.frame(t(x))))
write.csv(x_omit, file = "R_variables.csv", row.names=FALSE)