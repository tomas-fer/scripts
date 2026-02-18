# Instructions for running RepeatExplorer2 on MetaCentrum
(see REonMetaCentrum.txt for detailed descriptions)
---
__outputUNMAPPEDfomBAMasFASTQ.sh__
  * creates FASTQ.gz files with unmapped reads only
  * works on BAM files from HybPhyloMaker

__makeUnmappedReadsSummaryMETA.sh__
  * summary of the previous step
  * creates a file with the number of read pairs for every sample

__runREonMETA_parallel.sh__
  * creates job files for parallel RE run of multiple samples
  * needs list.txt

__makeREsummaryMETA.sh__
  * summary of the previous step
  * needs listGS.txt
  * creates tables and plots based on automated annotation (using REplot.R)

__runREComparativeOnMETA.sh__
  * comparative RE run
  * needs listComparative.txt
  * creates a comparative plot and a network based on the number of edges among samples (using REphylo.R)
