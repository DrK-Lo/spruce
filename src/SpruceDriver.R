### First: run bayenv2 for covariance matrices
bayenv2_spruce_makeCovarAllPops.txt

## First make the environment files for Spruce
MakeEnviFilesSpruce.R
  # note that for bayenv2, the environment file has to be split in 1/2 (130/2 = 65 variables each)

### Next edit the covariance matrices from Bayenv2
ssh klott@rogue.zoology.ubc.ca
cd /data/seqcap/spruce/bwa_pseudo/final_tables/bayenv2
CallEdCovMat.R

### Next check correlations among Bayenv2 covariance matrices
CompareCovmatBayenv2.R

### Next edit covariance matrix for missing environmental data
bayenv2_RemovePopsMissingData.R

### Next want to call SNPs on bugaboo
ssh klott@bugaboo.westgrid.ca
cd SpruceBayenv2
  # need script bayenv2_callSNPs.R in this folder
  # need covariance matrix var_out_588ind_all_v2_COMBINED.table.contig_flt80.noncoding.bayenv.covmatED.goodPops
  # need envi files seqcap.bayenv.BF.envi1 and seqcap.bayenv.BF.envi2
    # note that the envi variables are split in half, because bayenv2 can't handle all of them
  # need the list of SNPs var_out_spruce_all_beagle10_p99.tab.bayenv
R
source("bayenv2_makeSubScripts.R")

