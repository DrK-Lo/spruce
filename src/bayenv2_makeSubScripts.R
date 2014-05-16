### Write submission script for Westgrid
# ssh klott@bugaboo.westgrid.ca
# cd SpruceBayenv2

# for 5 loci, takes less than 1 minute
# for 20,000 loci, should take 2.8 days or 90 hours

# 15122420 lines in var_out_spruce_all_beagle10_p99.tab.bayenv file
# 7561210 snps
# 379 jobs with 20,000 loci each
# 

for (j in 1:379){

  startSNP<-(j-1)*20000 +1
  endSNP<- startSNP + 19999
  
  if (j==379){endSNP<-7561210 }

### practice snp  
  #  for (j in 1:1){  
  # startSNP=1
  #  endSNP=5
  SubmitScriptName <- paste("SubmitBayenvSNPS_Set",j,"_", startSNP,"_to_",endSNP, sep="") 

#practice  write(	c("#!/bin/bash", "#PBS -r n", "#PBS -m a", "#PBS -M klott@zoology.ubc.ca", "#PBS -l walltime=00:10:00", "#PBS -l procs=1", "#PBS -l mem=1gb", "" ,"cd $PBS_O_WORKDIR", paste("R --vanilla < bayenv2_callSNPs.R --args", startSNP, endSNP, j, sep=" ")), ncolumns=1, append=FALSE, file=SubmitScriptName)
	
		write(	c("#!/bin/bash", "#PBS -r n", "#PBS -m a", "#PBS -M klott@zoology.ubc.ca", "#PBS -l walltime=90:00:00", "#PBS -l procs=1", "#PBS -l mem=1gb", "" ,"cd $PBS_O_WORKDIR", paste("R --vanilla < bayenv2_callSNPs.R --args", startSNP, endSNP, j, sep=" ")), ncolumns=1, append=FALSE, file=SubmitScriptName)
	
		#Submit Submission Script to system	
		system(paste("qsub ", SubmitScriptName))

	}