# Instructions for running [RepeatExplorer2](http://repeatexplorer.org/) on [MetaCentrum](https://metavo.metacentrum.cz/en/index.html) (especially for off-target Hyb-Seq reads)
(see [REonMetaCentrum.txt](REonMetaCentrum.txt) for detailed descriptions)
---
__[outputUNMAPPEDfomBAMasFASTQ.sh](outputUNMAPPEDfomBAMasFASTQ.sh)__
  * creates FASTQ.gz files with unmapped reads only
  * works on BAM files from [HybPhyloMaker](https://github.com/tomas-fer/HybPhyloMaker)

__[makeUnmappedReadsSummaryMETA.sh](makeUnmappedReadsSummaryMETA.sh)__
  * summary of the previous step
  * creates a file with the number of read pairs for every sample

__[runREonMETA_parallel.sh](runREonMETA_parallel.sh)__
  * creates job files for parallel RE run of multiple samples
  * needs [list.txt](list.txt)

__[makeREsummaryMETA.sh](makeREsummaryMETA.sh)__
  * summary of the previous step
  * needs [listGS.txt](listGS.txt)
  * creates tables and plots based on automated annotation (using [REplot.R](REplot.R))

__[runREComparativeOnMETA.sh](runREComparativeOnMETA.sh)__
  * comparative RE run
  * needs [listComparative.txt](listComparative.txt)
  * creates a comparative plot and a network based on the number of edges among samples (using [REphylo.R](REphylo.R))
