#Processing trees from quartet sampling (Pease et al., 2018)
#output trees are modified with sed&grep and plotted in R to resemble trees in the publication
#WORKING VERSION, MAY BE NOT PROPERLY WORKING YET
#T. Fer, tomas.fer@natur.cuni.cz, 2019

#1. Modify the tree from command line
#tree with 'qc' values (used later for coloring nodes)
sed -e 's/\[[^][]*\]//g' RESULT.labeled.tre.qc > RESULT.labeled.tre.qc.modif #removes everything within '[ ]', i.e. df values after species names
sed -i 's/qc=//g' RESULT.labeled.tre.qc.modif #removes 'qc='
#tree with all values (used for plotting the tree and three scores
sed 's/\[\&[^][]*\,//g' RESULT.labeled.tre.figtree > RESULT.labeled.tre.figtree.modif #removes everything within [] except 'score=...' and '['
sed -i 's/\[[^][]*\]//g' RESULT.labeled.tre.figtree.modif #removes everything within '[ ]', i.e. df values after species names
sed -i 's/score=//g' RESULT.labeled.tre.figtree.modif #removes 'score='
sed -i 's/\]//g' RESULT.labeled.tre.figtree.modif #removes ']'
grep tree1 RESULT.labeled.tre.figtree.modif | sed 's/^.*=//' > RESULT.labeled.tre.figtree.modif.nwk #take only tree, i.e., make newick file

#2. Continue in R (requires 'phyloch' package and 'devtools')
#Label trees by node color points according to node label values
library(devtools) #for installing 'phyloch' package
library(phyloch) #for tree labelling
library(phytools) #for improved tree reading
#install_github("fmichonneau/phyloch")
treeqc<-read.newick("RESULT.labeled.tre.qc.modif")
# root tree (edgelabel=T for correct labels according to Czech et al.)
outgroup <- c("Phoenix-dactylifera_cds", "Hanguana-malayana_S373")
treeqc<-root(treeqc, edgelabel=T, outgroup)
tree<-read.newick("RESULT.labeled.tre.figtree.modif.nwk")
tree<-root(tree, edgelabel=T, outgroup)
plot.phylo(tree, align.tip.label=T, cex=0.8) #plot tree
nodelabels(tree$node.label, frame ="none", adj=-0.5, cex=0.5) #add labels (the three
#color nodes according to 'qc' values [-1,1]
dotsize=1.2
node.support(treeqc$node.label, cutoff = -1, cex=dotsize, mode = "dots", col="red") #all labelled by red
node.support(treeqc$node.label, cutoff = -0.05, cex=dotsize, mode = "dots", col="darkorange") #> -0.05 labelled by dark orange
node.support(treeqc$node.label, cutoff = 0, cex=dotsize, mode = "dots", col="lightgreen") #> 0 labelled by light green
node.support(treeqc$node.label, cutoff = 0.2, cex=dotsize, mode = "dots", col="darkgreen") #> 2 labelled by dark green

