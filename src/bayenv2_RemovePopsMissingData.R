

covpath<- "~/SpruceBayenv2/var_out_588ind_all_v2_COMBINED.table.contig_flt80.noncoding.bayenv.covmatED"
#covpath<-"~/Desktop/AdaptreeData/spruce/covmats/var_out_588ind_all_v2_COMBINED.table.contig_flt80.noncoding.bayenv.covmatED"

envipath<-"~/SpruceBayenv2/envifiles/seqcap.bayenv.BF.envi"
missingPopsPath <- "~/SpruceBayenv2/envifiles/seqcap.bayenv.BF.envi.linesRemoved"
#missingPopsPath <- "~/Desktop/AdaptreeData/spruce/envifiles/seqcap.bayenv.BF.envi.linesRemoved"

covmat <- read.table(covpath)
missingPops<- scan(missingPopsPath)
covmat2<- covmat[-missingPops,-missingPops]

dim(covmat2)

  write.table(covmat2, file=paste(covpath,".goodPops", sep=""), sep="\t", eol="\t\n", col.names=FALSE, row.names=FALSE)
 