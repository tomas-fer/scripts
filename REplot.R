#plot barplots for RE results (multiple species, general groups) in basic R
#i.e., combine 4th columns from every name_sumpropsimple.txt, keep repeat names as 1st column
#first line = species names without spaces (or every name has to be quoted)

a <- read.table("RE_sumpropgenMb.txt", row.names = 1, header=T)
m <- as.matrix(a,rownames = TRUE)
m[is.na(m)] = 0 #change NA to 0
#remove lines with small clusters, contamination and singlets
remove <- c("contamination","smallClusters","singlets")
m2 <- m[!rownames(m) %in% remove, ]
m2perc <- apply(m2, 2, function(x){x*100/sum(x,na.rm=T)}) #calculate percentage
cols <- c("darkorange2","cornflowerblue","burlywood4","yellow3","forestgreen","pink3","red3","grey","grey89","white")
icols <- c("#e69570","#7aa3a7","#A77B7E","#d3cc33","#6a943d","#dd8da1","#e75a68","#99a099","#d5dbde","#f8fcfb")

#plot all
pdf(file="REgen.pdf")
par(mar=c(5,10,5,3)) #set margins (bottom, left, top, right) - might need to adjust
par(xpd=TRUE) #allow plotting legend outside the plot area
barplot(m, horiz = TRUE, col=icols, las = 1) #las=1 - labels horizontally
legend("topright", legend=rownames(m), fill=icols, ncol=3, inset=c(0,-0.18), bty="n") #inset need to be adjusted
title(xlab = "Genome size [Mb/1C]")
dev.off()

#plot percentage
pdf(file="REgenPerc.pdf")
par(mar=c(5,10,5,3)) #set margins (bottom, left, top, right) - might need to adjust
par(xpd=TRUE) #allow plotting legend outside the plot area
barplot(m2perc, horiz = TRUE, col=icols, las = 1) #las=1 - labels horizontally
legend("topright", legend=rownames(m2), fill=icols, ncol=3, inset=c(0,-0.18), bty="n") #inset need to be adjusted
title(xlab = "Repeat proportion [%]")
dev.off()

a <- read.table("RE_sumpropsimpleMb.txt",row.names = 1, header=T)
m <- as.matrix(a,rownames = TRUE)
m[is.na(m)] = 0 #change NA to 0
#remove lines with small clusters, contamination and singlets
remove <- c("contamination","smallClusters","singlets")
m2 <- m[!rownames(m) %in% remove, ]
m2perc <- apply(m2, 2, function(x){x*100/sum(x,na.rm=T)}) #calculate percentage
cols <- c("darkorange2","cornflowerblue","burlywood4","yellow3","forestgreen","pink3","red3","grey","grey89","white")
icols <- c("#e69570","#7aa3a7","#A77B7E","#d3cc33","#6a943d","#dd8da1","#e75a68","#99a099","#d5dbde","#f8fcfb")

#plot all
pdf(file="REsimple.pdf")
par(mar=c(5,10,5,3)) #set margins (bottom, left, top, right) - might need to adjust
par(xpd=TRUE) #allow plotting legend outside the plot area
barplot(m, horiz = TRUE, col=icols, las = 1) #las=1 - labels horizontally
legend("topright", legend=rownames(m), fill=icols, ncol=3, inset=c(0,-0.18), bty="n") #inset need to be adjusted
title(xlab = "Genome size [Mb/1C]")
dev.off()

#plot percentage
pdf(file="REsimplePerc.pdf")
par(mar=c(5,10,5,3)) #set margins (bottom, left, top, right) - might need to adjust
par(xpd=TRUE) #allow plotting legend outside the plot area
barplot(m2perc, horiz = TRUE, col=icols, las = 1) #las=1 - labels horizontally
legend("topright", legend=rownames(m2), fill=icols, ncol=3, inset=c(0,-0.18), bty="n") #inset need to be adjusted
title(xlab = "Repeat proportion [%]")
dev.off()

