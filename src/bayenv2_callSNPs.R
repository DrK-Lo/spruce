#scp -r klott@rogue.zoology.ubc.ca:/data/seqcap/spruce/bwa_pseudo/final_tables/spruce_GBS_all.table.contig_flt10.bayenv* .
#scp -r klott@rogue.zoology.ubc.ca:/data/seqcap/spruce/bwa_pseudo/final_tables/var_out_spruce_all_beagle10_p99.tab.bayenv* .
#also transfer edited covariance matrix from laptop
 #var_out_588ind_all_v2_COMBINED.table.contig_flt80.noncoding.bayenv.covmatED.goodPops
#first source bayenv2_RemovePopsMissingData.R

args <- commandArgs(trailingOnly = TRUE) 
startSNP <- as.numeric(args[1])
endSNP <- as.numeric(args[2])
groupnum <- as.numeric(args[3])
rm(args)

#take argument for SNPpath

SNPpath<-"~/SpruceBayenv2/var_out_spruce_all_beagle10_p99.tab.bayenv"
#ROGUE SNPpath<-"/data/seqcap/spruce/bwa_pseudo/final_tables/var_out_spruce_all_beagle10_p99.tab.bayenv"

  # this is not edited to match environment.  editing happens below.
envipath1<-"~/SpruceBayenv2/envifiles/seqcap.bayenv.BF.envi1"
envipath2<-"~/SpruceBayenv2/envifiles/seqcap.bayenv.BF.envi2"
#ROGUE envipath1<-"/data/seqcap/spruce/bwa_pseudo/final_tables/envifiles/seqcap.bayenv.BF.envi1"
# envipath2<-"/data/seqcap/spruce/bwa_pseudo/final_tables/envifiles/seqcap.bayenv.BF.envi2"

covpath<- "~/SpruceBayenv2/var_out_588ind_all_v2_COMBINED.table.contig_flt80.noncoding.bayenv.covmatED.goodPops"
#ROGUE covpath <- "/data/seqcap/spruce/bwa_pseudo/final_tables/bayenv2/var_out_588ind_all_v2_COMBINED.table.contig_flt80.noncoding.bayenv.covmatED.goodPops"
  #MCMC 85500
  #this matrix is edited to match environment
SNPnames <- "~/SpruceBayenv2/var_out_spruce_all_beagle10_p99.tab.bayenv.loci"
#ROGUE SNPnames<- "/data/seqcap/spruce/bwa_pseudo/final_tables/var_out_spruce_all_beagle10_p99.tab.bayenv.loci"

PopsMissingData <- "~/SpruceBayenv2/envifiles/seqcap.bayenv.BF.envi.linesRemoved"
#ROGUE PopsMissingData <- "/data/seqcap/spruce/bwa_pseudo/final_tables/envifiles/seqcap.bayenv.BF.envi.linesRemoved"

enviDims1 <- dim(read.table(envipath1))
npops <- enviDims1[2]
nenvi1 <- enviDims1[1]

enviDims2 <- dim(read.table(envipath2))
nenvi2 <- enviDims2[1]

missingPops<- scan(PopsMissingData)

npopsUNED<- length(scan(SNPpath, nlines=1))

OUTFILE1<- paste("var_out_spruce_all_beagle10_p99.tab.bayenv.envi1.OUT",groupnum, sep="")
write(NULL, paste(OUTFILE1,".bf",sep=""), append=FALSE)
write(NULL, paste(OUTFILE1,".xtx",sep=""), append=FALSE)

OUTFILE2<- paste("var_out_spruce_all_beagle10_p99.tab.bayenv.envi2.OUT",groupnum, sep="")
write(NULL, paste(OUTFILE2,".bf",sep=""), append=FALSE)
write(NULL, paste(OUTFILE2,".xtx",sep=""), append=FALSE)

for (i in startSNP:endSNP){
  skp <- i*2-2
  
  thisSNPname <- scan(SNPnames, what="character", skip=(i-1), 
                      nlines=1, quiet=TRUE)[2]
  
  # get SNP data
  SNPdata1<- matrix(scan(SNPpath, skip=skp, nlines=2,
                        quiet=TRUE), byrow=TRUE, ncol=npopsUNED)
    #note that not all populations have environmental data
  SNPdata<-SNPdata1[,-missingPops]
  
  # write SNP file
  write.table(SNPdata, file=thisSNPname, sep="\t", eol="\t\n",
              col.names=FALSE, row.names=FALSE)
  
  RANDOM <- sample(1:99999,1)
  
  # call bayenv2
  #  script<-paste("./bayenv2 -i", thisSNPname, "-m", covpath, "-e", envipath, "-p", npops,
  #  "-k 10000", "-n", nenvi, "-t -c -X -r", RANDOM, "-o", OUTFILE)
  script.envi1<-paste("./bayenv2 -i", thisSNPname, "-m", covpath, "-e", envipath1, "-p", npops,
    "-k 10000", "-n", nenvi1, "-t -c -X -r", RANDOM, "-o", OUTFILE1)
    print(script.envi1)
  system(script.envi1)
  
    script.envi2<-paste("./bayenv2 -i", thisSNPname, "-m", covpath, "-e", envipath2, "-p", npops,
    "-k 10000", "-n", nenvi2, "-t -c -X -r", RANDOM, "-o", OUTFILE2)
    print(script.envi2)
    system(script.envi2)
  
  # remove SNP file
  system(paste("rm ", thisSNPname, "*", sep=""))
  
}