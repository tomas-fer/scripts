#Processing scored tree from ASTRAL (with the option -t 4)
#This creates a PNG image of the rooted ASTRAL species tree with pie charts on branches
#See ASTRAL manual for explanation
#T. Fer, tomas.fer@natur.cuni.cz, 2020

require(treeio) #latest version from GitHub necessary? (implements 'read.astral')
require(ape)

#Import the tree generated with ASTRAL -t 4
x <- read.astral("Astral_50_75_t4.tre")
#Save posterior probabilities 
pp1 <- as.numeric(x[['pp1']])
pp2 <- as.numeric(x[['pp2']])
pp3 <- as.numeric(x[['pp3']])
pc <- cbind( pp1, pp2, pp3)

tree <- root(x@phylo, "Phoenix-dactylifera_cds") #root tree

#Plot/save the tree (modify 'adj' and 'cex' values according to your needs)
png("Astral50_t4.png", width = 2480, height = 2480)
#par(lty = "blank") #no lines in pie charts
plot(tree, cex = 3, label.offset = 1)
nodelabels(pie = pc, adj = c(-0.5, 1.1), piecol = c("cornflowerblue","cadetblue3","lightblue"), cex = 0.5)
#nodelabels(pie = pc, adj = c(-0.5, 1.1), piecol = c("blue","red","green"), cex = 0.5)
dev.off()
