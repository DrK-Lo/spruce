### bayenv2: Each environmental variable should be standardized, i.e. subtract 
### the mean and then divided through by the standard deviation of the variable 
### across populations. Each line of the file is an environmental
### variable across populations, values should be tab separated. The variables should 
### be listed in the same population order as they appear in the allele count files.

### lfmm The variable file is a vector composed of n lines and D columns. 
### Each line is the values of the D variables for the corresponding individual. 
### Below, an example of variable file for n = 3 individuals and D = 1 covariable.
### Warning: If you set several covariables, the program will be launched for 
### each covariable sequentially and independently.
### 0.252477
### 0.216618
### -0.47509

### load environment file from Tongli
setwd("~/Desktop/AdaptreeData")
setwd("spruce")

envi <- read.csv("SpruceAllEnvis.csv")
tail(envi)

goodpops=NULL
for (i in 1:nrow(envi)){
  s1 <- sum(which(envi[i,]==-9999), is.na(envi[i,]))
  if (s1==0){goodpops <- c(goodpops,i)}
}

popIDs <- paste("Sp_",envi$WorkingSL,sep="")
goodpopIDs <- popIDs[goodpops]
envicols<-5:ncol(envi)
goodpopIDs

std.envi <- envi
for (a in 1:length(envicols)){
  ec <- envi[goodpops,envicols[a]]
  ave <- mean(ec)
  stdev <- sd (ec)
  std.envi[goodpops,envicols[a]]<- (ec-ave)/stdev
}

bayenv.poporder.seqcap <- read.table("seqcap.bayenv.pop_order")
bayenv.poporder.gbs <- read.table("gbs.bayenv.pop_order")
lfmm.poporder.seqcap <- read.table("seqcap.lfmm.pop_order.txt")
lfmm.poporder.gbs <- read.table("gbs.lfmm.pop_order.txt")
tail(lfmm.poporder.gbs)

GetEnviMat <- function(poplist){
  # poplist is a vector from the dataset
  # this function returns the environment matrix corresponding to poplist
  # Note that some variables are defined above
  
  goodInds <- which(poplist %in% goodpopIDs)
  popNames <- popIDs[goodInds]
  linesRemoved <- 1:length(poplist)
  linesRemoved<-linesRemoved[-goodInds]
  ncols <- length(envicols)
  envimat <- matrix(NA, nrow=length(goodInds), ncol=ncols)
  
  for (i in 1:length(goodInds)){
    index<-which(goodpopIDs==poplist[i])
    envimat[i,]<-as.numeric(std.envi[index, envicols])
  }

  return(list(envimat=envimat, popNames=popNames, linesRemoved=linesRemoved))
}


### write the column names for the environments
  write.table(colnames(envi)[envicols], "envifiles/enviNamesAllAnalyses.txt", col.names=FALSE,
            row.names=FALSE, sep="\t", eol="\n")
  

### For bayes factor calculation, can only use populations with environments
### Will have to edit the covariance matrix to remove last 4 populations
  bayenv.seqcap.envi<-GetEnviMat(bayenv.poporder.seqcap$V2)
  sum(is.na(bayenv.seqcap.envi$envimat))
  dim(bayenv.seqcap.envi$envimat)
  length(bayenv.poporder.seqcap$V2)
  bayenv.seqcap.envi$linesRemoved
  bayenv.poporder.seqcap$V2[bayenv.seqcap.envi$linesRemoved]
  # my envimat has popluations in rows and envis in columns
  # Bayenv2 wants envis in rows and populations in columns
  write.table(t(bayenv.seqcap.envi$envimat), "envifiles/seqcap.bayenv.BF.envi", col.names=FALSE,
              row.names=FALSE, sep="\t", eol="\t\n")
  write.table(bayenv.seqcap.envi$popNames, "envifiles/seqcap.bayenv.BF.envi.popnames",col.names=FALSE,
              row.names=FALSE)
  write.table(bayenv.seqcap.envi$linesRemoved, "envifiles/seqcap.bayenv.BF.envi.linesRemoved",col.names=FALSE,
            row.names=FALSE)
  ### For XTX calculation, can use all pops (don't need envi)

  bayenv.gbs.envi <- GetEnviMat(bayenv.poporder.gbs$V2)
  sum(is.na(bayenv.gbs.envi$envimat))
  dim(bayenv.gbs.envi$envimat)
  length(bayenv.poporder.gbs$V2)
  bayenv.gbs.envi$linesRemoved
  bayenv.poporder.gbs$V2[bayenv.gbs.envi$linesRemoved]
  # my envimat has popluations in rows and envis in columns
  # Bayenv2 wants envis in rows and populations in columns
  write.table(t(bayenv.gbs.envi$envimat), "envifiles/gbs.bayenv.BF.envi", col.names=FALSE,
            row.names=FALSE, sep="\t", eol="\t\n")
  write.table(bayenv.gbs.envi$popNames, "envifiles/gbs.bayenv.BF.popnames", col.names=FALSE,
            row.names=FALSE)
  write.table(bayenv.gbs.envi$linesRemoved, "envifiles/gbs.bayenv.BF.envi.linesRemoved",col.names=FALSE,
            row.names=FALSE)


### For lfmm calculation, can only use populations with environments
### Will have to edit the lfmm file to remove the relevant lines
  lfmm.seqcap.envi <- GetEnviMat(lfmm.poporder.seqcap$V2)
  sum(is.na(lfmm.seqcap.envi$envimat))
  dim(lfmm.seqcap.envi$envimat)
  length(lfmm.poporder.seqcap$V2)
    lfmm.seqcap.envi$linesRemoved
  lfmm.poporder.seqcap$V2[lfmm.seqcap.envi$linesRemoved]
  # lfmm wants individuals in rows and environments in columns: matches my data
  write.table(lfmm.seqcap.envi$envimat, "envifiles/seqcap.lfmm.BF.envi", col.names=FALSE,
            row.names=FALSE, sep="\t", eol="\t\n")
  write.table(lfmm.seqcap.envi$popNames, "envifiles/seqcap.lfmm.BF.popnames", col.names=FALSE,
            row.names=FALSE)
  write.table(lfmm.seqcap.envi$linesRemoved, "envifiles/seqcap.lfmm.BF.envi.linesRemoved",col.names=FALSE,
            row.names=FALSE)

  lfmm.gbs.envi <- GetEnviMat(lfmm.poporder.gbs$V2)
  sum(is.na(lfmm.gbs.envi$envimat))
  dim(lfmm.gbs.envi$envimat)
  length(lfmm.poporder.gbs$V2)
    lfmm.gbs.envi$linesRemoved
  lfmm.poporder.gbs$V2[lfmm.gbs.envi$linesRemoved]
  # lfmm wants individuals in rows and environments in columns: matches my data
  write.table(lfmm.gbs.envi$envimat, "envifiles/gbs.lfmm.BF.envi", col.names=FALSE,
            row.names=FALSE, sep="\t", eol="\t\n")
  write.table(lfmm.gbs.envi$popNames, "envifiles/gbs.lfmm.BF.popnames", col.names=FALSE,
            row.names=FALSE)
  write.table(lfmm.gbs.envi$linesRemoved, "envifiles/gbs.lfmm.BF.envi.linesRemoved",col.names=FALSE,
            row.names=FALSE)
