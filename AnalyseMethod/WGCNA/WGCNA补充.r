### WGCNA补充1
ModuleTree = hclust(dist(t(net$MEs)), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/ModuleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(ModuleTree, main = "Module clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 1.2, col = "red");
# Determine cluster under the line
# 保存图片 Module clustering to detect outliers
Moduleclust = cutreeStatic(ModuleTree, cutHeight = 1.2, minSize = 10)
table(Moduleclust)





### WGCNA补充2
library(ape)
#calculate eigengenes
MEs = moduleEigengenes(WGCNA_matrix, colors = module.colours, excludeGrey = FALSE)$eigengenes

#calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

#cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = 'average');

#plot the result with phytools package
par(mar=c(2,2,2,2))
plot.phylo(as.phylo(METree),type = 'fan',show.tip.label = FALSE, main='')
tiplabels(frame = 'circle',col='black', text=rep('',length(unique(modules))), bg = levels(as.factor(module.colours)))