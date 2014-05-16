#klott@rogue:/data/seqcap/spruce/bwa_pseudo/final_tables/bayenv2
library(fields)
source("http://bioconductor.org/biocLite. R") 
biocLite("made4")
require(made4)

setwd("~/Desktop/AdaptreeData/spruce/covmats")
files <- read.table("../SpFileInfo.txt",header=TRUE)
filens<- paste(files$fileBaseName, ".bayenv.covmatED",sep="")
filens[1]<- "COMBINED.table.contig_flt10.bayenvED"  #will need to fix this later

pops <- read.table("../seqcap.bayenv.pop_order")

## Make heatmaps 

for (i in 1:length(filens)){
  npops <- length(scan(filens[i], nlines=1, skip=1))
  covmat <- matrix(scan(filens[i]), ncol=npops, byrow=TRUE)
  dim(covmat)
  rownames(covmat)<-pops$V2
  colnames(covmat)<-pops$V2

  pdf(paste("pdfs/covmats/",filens[i],".pdf",sep=""),width=10, height=9, bg="white")
    par(mar=c(4,4,4,7))
    image(covmat, axes=FALSE, col=heat.colors(100)) 
    every10 <- seq(1,npops,10)
    axis(1,at=seq(0,1-1/npops, 1/(npops/10)),lab=as.character(pops$V2[every10]), 
         tick=FALSE, hadj=1, padj=0, line=0, las=2, cex=0.5)
    axis(2,at=seq(0,1-1/npops, 1/(npops/10)),lab=as.character(pops$V2[every10]), 
         tick=FALSE, hadj=1, padj=0, line=0, las=2, cex=0.5)
    mtext(as.character(filens[i]),side=3,outer=FALSE, line=2)
    box()
    image.plot(covmat, legend.only=TRUE,col=heat.colors(100))
  dev.off()
}

cor.table <- matrix(NA, nrow=length(filens), ncol=length(filens))
colnames(cor.table)=filens
rownames(cor.table)=filens

slope.table <- cor.table

par(mar=c(4,4,1,1))

nx<-rep(NA, length(filens)^2) 

v1<- nx
v2<- nx
cor.out<- nx
intercept.out<- nx
slope.out<- nx
i.info <- data.frame(MinimumData=nx, Imputation=nx, SNPs=nx, Ascertainment=nx) 
j.info <- data.frame(MinimumData=nx, Imputation=nx, SNPs=nx, Ascertainment=nx) 

### Need to edit covariance matrix from GBS ###
extraPop <- read.table("../bayenv_gbs_extraPop.pop_order", header=FALSE)
xtraPop <- as.numeric(extraPop[1])

k=0
for (i in 1:9){
  if (i<9){
    covmat_i <- matrix(scan(filens[i]), ncol=253, byrow=TRUE)
  }
  if (i==9){
    covmat_gbs <- matrix(scan(filens[i]), ncol=254, byrow=TRUE)
    # remove the extra population from GBS
    covmat_i <- covmat_gbs[-xtraPop,-xtraPop]
  }
  
  for (j in 1:9){
      k=k+1
      i.info[k,]<-unlist(files[i,2:5])
      j.info[k,]<-unlist(files[j,2:5])
      v1[k]<-filens[i]
      v2[k]<-filens[j]
        if (j<9){
          covmat_j <- matrix(scan(filens[j]), ncol=253, byrow=TRUE)
        }
      if (j==9){
          covmat_gbs <- matrix(scan(filens[j]), ncol=254, byrow=TRUE)
          # remove the extra population from GBS
          covmat_j <- covmat_gbs[-xtraPop,-xtraPop]
        }
    c_out <- cor.test(covmat_i, covmat_j, method="pearson")
    cor.table[i,j]<-c_out$est
    cor.out[k]<-c_out$est
    if(i!=j){
      m<-lm(as.numeric(covmat_j)~as.numeric(covmat_i))
      intercept.out[k]<-m$coef[1]
      slope.out[k]<-m$coef[2]
      slope.table[i,j]<-m$coef[2]
      jpeg(paste("pdfs/regressions/",filens[i],"_VS_",filens[j],".jpeg", sep=""),width=1000,height=800)
        plot(covmat_i,covmat_j, xlim=c(0.4,0.6), ylim=c(0.4,0.6), xlab=filens[i], 
              ylab=filens[j], bty="l")
        abline(m, col="red", lwd=2)
        abline(0,1, lwd=5)
     dev.off()
    }
  }
}

colnames(i.info)<- paste("i_",colnames(i.info))
colnames(j.info)<- paste("j_",colnames(j.info))

out<- data.frame(i.info, j.info, v1x=v1,v2y=v2, cor.out, intercept.out, slope.out)
write.csv(out, file="CompareCovMatsResults.csv")

cor.table2<- data.frame(files$MinimumData,files$Imputation,files$SNPs, 
                        files$Ascertainment, files$fileBaseName, cor.table)
slope.table2<- data.frame(files$MinimumData,files$Imputation,files$SNPs, 
                        files$Ascertainment, files$fileBaseName, slope.table)

write.csv(cor.table2, file="CorTable.csv")
write.csv(slope.table2, file="SlopeTable.csv")
