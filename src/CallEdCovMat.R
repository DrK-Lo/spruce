setwd("/Users/katie/Desktop/AdaptreeData/spruce")
source("MakeEdCovMat.R")
files <- read.table("SpFileInfo.txt",header=TRUE)

filens<- paste(files$fileBaseName, ".bayenv.covmat",sep="")
filens[1]<- "COMBINED.table.contig_flt10.bayenv"  #will need to fix this later

for (i in 1:length(filens)){
  print(i)
  EditCovMat(filens[i])
}  