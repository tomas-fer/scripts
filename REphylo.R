#!/usr/bin/env Rscript
#--------------------------------------------------------------------------------------------
# HybPhyloMaker: create consensus network from RepeatExplorer2 comparative results
# https://github.com/tomas-fer/HybPhyloMaker
# v.1.8.0
# Tomas Fer, 2025
# tomas.fer@natur.cuni.cz
#--------------------------------------------------------------------------------------------
#go to the folder seqclust/clustering/clusters and run following commands to create a folder with all renamed matrices:
#mkdir OEnumberEdges
#for i in $(ls -d dir*); do cp $i/observed_expected_number_of_edges.csv OEnumberEdges/${i}.csv; done
#cd OEnumberEdges
#run as ./createTrees.R nrsamples

library(ape)
#read arguments from command line
args <- commandArgs()
nrsamples <- as.numeric(args[6])
mat_files <- dir(pattern="*csv") #create a list of all files

#all trees with at least nrsamples
TreeFromDistMat <- function(file) {
  clno <- strsplit(file, "_")
  clno <- strsplit(clno[[1]][2], "\\.")[[1]][1]
  #print(clno)
  mat = read.delim(file)
  rownames(mat) <- mat[,1] #get rownames from the first column
  mat <- mat[,!names(mat) %in% c("species")] #remove column species
  mat[mat == 0] <- NA #replace 0s by NAs
  mat <- mat[rowSums(is.na(mat)) != ncol(mat), ] #remove lines entirely with NAs
  mat <- mat[,colSums(is.na(mat)) != nrow(mat), ] #remove columns entirely with NAs
  mat <- 1/mat #inverse all numbers
  if(length(mat)<nrsamples){
    return(NA)
  }
  else {
    tryCatch(
        {
        tree <- njs(as.matrix(mat))
        return(tree)
        },
        error=function(e) {
            message('An Error Occurred')
            print(e)
            return(NA)
        },
        warning=function(w) {
            message('A Warning Occurred')
            print(w)
            return(NA)
        }
    )
  }
}

#run function over all RE comparative matrices
print("Creating NJ trees...")
treesRE <- lapply(mat_files, TreeFromDistMat) #run the function over the list of trees
treesRE <- treesRE[!is.na(treesRE)] #remove NAs produced by errors/warnings
write.tree(treesRE, file=paste0("trees",nrsamples,".tre"))

#calculate/plot consensus network (for thresholds from 0.05 to 0.15, step 0.01)
library(phangorn)
tr <- read.tree(file=paste0("trees",nrsamples,".tre"))
for (x in seq(0.05, 0.15, 0.01)) {
  print(c("Consensus network with threshold: ",x))
  cnet <- consensusNet(tr, x)
  write.nexus.networx(cnet,file=paste0("cnet",x,".nex"))
  a=15
  par(mar=c(a,a,a,a),xpd=T)
  pdf(file=paste0("cnet",x,".pdf"))
  plot(cnet,cex=0.7,tip.color='firebrick4',edge.width=0.3,edge.lty=1,direction='axial')
  dev.off()
}
