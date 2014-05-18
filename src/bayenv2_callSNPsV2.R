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

SNPnames1 <- matrix(scan("~/SpruceBayenv2/var_out_spruce_all_beagle10_p99.tab.bayenv.loci", skip=(startSNP-1), nlines=(endSNP-startSNP+1), what="character"), ncol=2, byrow=TRUE)
	#this takes about 10 minutes to read in
SNPnames<- SNPnames1[,2]
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
#write(NULL, paste(OUTFILE1,".bf",sep=""), append=FALSE)
#write(NULL, paste(OUTFILE1,".xtx",sep=""), append=FALSE)

OUTFILE2<- paste("var_out_spruce_all_beagle10_p99.tab.bayenv.envi2.OUT",groupnum, sep="")
#write(NULL, paste(OUTFILE2,".bf",sep=""), append=FALSE)
#write(NULL, paste(OUTFILE2,".xtx",sep=""), append=FALSE)

SNPdata1<- matrix(scan(SNPpath, skip=(startSNP-1), nlines=2*(endSNP-startSNP+1), quiet=TRUE), byrow=TRUE, ncol=npopsUNED)

if (file.exists(paste(OUTFILE2,".xtx",sep=""))){
	o <- read.table(paste(OUTFILE2,".xtx",sep=""))
	vect<- which((SNPnames %in% o$V1)==FALSE)
}else{
		vect<- 1:(endSNP-startSNP+1)
}


for (i in vect){
	
	print("____________________________________________________")
	print(i+startSNP-1))
  
  thisSNPname <- SNPnames[i]
  
  # get SNP data
	SNPdata2<- SNPdata1[(i*2-1):(i*2),]
    #note that not all populations have environmental data
  SNPdata<-SNPdata2[,-missingPops]
  
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