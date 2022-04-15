#install.packages("gplots")

library(gplots)


data1 <- read.csv("SARS_CoV_2_Origin_project_FullTable_1.csv", header=T)

data1f <- t(data1[20050:21100, 2:15])
data2 <- data.matrix(data1[21050:21100,2:15])
data2t <- t(data2)

colors <- c("white", "lightgreen", "pink", "lightblue", "yellow")
#lmat <- rbind( c(5,3,4), c(2,1,4) )
lhei <- c(1.5, 4)

#heatmap(data2t, Rowv = NA, Colv = NA, col=colors)

pdf("seq_1.pdf")
heatmap.2(data2t, notecex=0.7, density.info="none", cellnote = data1f, notecol="black", cexRow = 0.5,lhei = c(1.5, 4.0), keysize = 1, trace="none", margins = c(12,9), col = colors, hline = NA, vline = NA, dendrogram = "none",Rowv = FALSE, Colv = FALSE)

dev.off()