# source("https://bioconductor.org/biocLite.R")
# biocLite(c("AnnotationDbi", "impute","GO.db", "preprocessCore"))
# site="https://mirrors.tuna.tsinghua.edu.cn/CRAN"
# install.packages(c("WGCNA", "stringr", "reshape2"), repos=site)

library(WGCNA)
library(dynamicTreeCut)
library(fastcluster)


## 输入表达矩阵
exprMat <- "tumor_sample_DEG_Matrix.txt"
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)



# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)

# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)

##导入数据##
dataExpr <- read.table(exprMat, sep='\t', row.names=1, header=T, 
                     quote="", comment="", check.names=F)

dim(dataExpr)
head(dataExpr)[,1:8]



## 数据筛选
## 筛选中位绝对偏差前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
# apply函数只能用于处理矩阵类型的数据，也就是说所有的数据必须是同一类型。因此要使用apply函数的话，需要将数据类型转换成矩阵类型。
# apply函数一般有三个参数，第一个参数代表矩阵对象，第二个参数代表要操作矩阵的维度，1表示对行进行处理，2表示对列进行处理。第三个参数就是处理数据的函数。apply会分别一行或一列处理该矩阵的数据。
# MAD（Median absolute deviation, 中位数绝对偏差）是单变量数据集中样本差异性的稳健度量。mad是一个健壮的统计量，对于数据集中异常值的处理比标准差更具有弹性，可以大大减少异常值对于数据集的影响。
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad > 
                 max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]

## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))

for (i in 1:ncol(dataExpr)){
	dataExpr[,i] = as.numeric(dataExpr[,i])
}

## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(dataExpr)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ", ")));
 # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}




## 查看是否有离群样品
sampleTree = hclust(dist(dataExpr), method = "average");
sizeGrWindow(12,9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# 保存图片 Sample clustering to detect outliers




abline(h = 8000, col = "red");
clust = cutreeStatic(sampleTree, cutHeight = 8000, minSize = 10)
table(clust)
# clust
# 0			
# 1077

keepSamples = (clust==0) ## 剔除离群样本
dataExpr = dataExpr[keepSamples, ]
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# 保存图片Scale independence & Mean connectivity


## 一般不使用推荐power
# power = sft$powerEstimate
# power
# [1] 6



## 经验power (无满足条件的power时选用)
# 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使
# 无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
if (is.na(power)){
  power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
          ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
          ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
          ifelse(type == "unsigned", 6, 12))       
          )
          )
}




### 网络构建
## 一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值
# maxBlockSize: 计算机能处理的最大模块的基因数量 (默认5000)；
#  4G内存电脑可处理8000-10000个，16G内存电脑可以处理2万个，32G内存电脑可
#  以处理3万个
#  计算资源允许的情况下最好放在一个block里面。
# corType: pearson or bicor
# numericLabels: 返回数字而不是颜色作为模块的名字，后面可以再转换为颜色
# saveTOMs：最耗费时间的计算，存储起来，供后续使用
# mergeCutHeight: 合并模块的阈值，越大模块越少
# 必须做数值化
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType, 
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = "expr_d5_TOM",
                       verbose = 3)

# 根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
# **0 (grey)**表示**未**分入任何模块的基因。 
table(net$colors)
  # 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15 
# 767 258 168 142 127  70  64  63  57  53  52  48  46  46  37  35 

length(net$colors)
length(table(net$colors))
mergedColors = labels2colors(net$colors)
table(mergedColors)


## 层级聚类树展示各个模块
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
# 保存图片 Cluster Dendrogram





# WGCNA补充1
# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs

### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
library(stringr)
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

ModuleTree = hclust(dist(t(MEs_col)), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/ModuleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(ModuleTree, main = "Module clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 1.2, col = "red");
Moduleclust = cutreeStatic(ModuleTree, cutHeight = 1.2, minSize = 10)
table(Moduleclust)
# Determine cluster under the line
# 保存图片 Module clustering to detect outliers




#source("http://bioconductor.org/biocLite.R")
#biocLite("marray")
# biocLite("limma")
library(marray)
library(limma)
write.list(net,"expredata_WGCNA_net.txt")
write.table(dataExpr,"expredata_WGCNA_result.txt",quote=FALSE,sep="\t")
### 模块表格建立方法：net.txt第一行(colors)含有3180个数字，代表模块号，分列后一一匹配对应到result.txt的第一行(Gene)即可得到模块列表
# 得到文件Module.txt




moduleColors = labels2colors(net$colors)
TOM = 1-TOMsimilarityFromExpr(dataExpr, power = 5);
plotTOM = TOM^7;
diag(plotTOM) = NA;
sizeGrWindow(9,9)
geneTree = net$dendrograms[[1]];
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
# 保存图片 Network heatmap plot, all genes




# Read in the annotation file
# annot = read.csv(file = "GeneAnnotation.csv");
# Select modules需要修改，选择需要导出的模块颜色
modules = unique(mergedColors);
# Select module probes选择模块探测
probes = colnames(dataExpr)
inModule = is.finite(match(mergedColors, modules));
modProbes = probes[inModule];
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("expredata_edges",".txt", sep=""),
                               nodeFile = paste("expredata_node",".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.1,
                               nodeNames = modProbes,
                               #altNodeNames = modGenes,
                               nodeAttr = mergedColors[inModule]);
# 使用expredata_node.txt建立模块基因文本Module_Gene.txt





## WGCNA补充2
#biocLite("ape")
library(ape)
par(mar=c(2,2,2,2))
 # par(mar=c(0.1,0.1,0.1,0.1))
#calculate eigengenes
MEs = moduleEigengenes(dataExpr, colors = moduleColors, excludeGrey = FALSE)$eigengenes

#calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);

#cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = 'average');

#plot the result with phytools package
# par(mar=c(2,2,2,2))
plot.phylo(as.phylo(METree),type = 'fan',show.tip.label = FALSE, main='')
tiplabels(frame = 'circle',col='black', text=rep('',length(unique(modules))), bg = levels(as.factor(moduleColors)))
# 保存图片 





### 建立模块基因文本Module_Gene.txt
## 剔除灰色（grey）未分类到模块的基因
expredata_node <- read.table("expredata_node.txt",sep = "\t", stringsAsFactors = F,header = T, fill=TRUE)
expredata_node = expredata_node[,-2]
colnames(expredata_node) = c("Symbol","type")

Module = data.frame(data.frame(table(mergedColors))[order(data.frame(table(mergedColors))[,2],decreasing=T),])## 默认升序，decreasing=T时降序
Module <- data.frame(Module=unique(Module[,1]))

library(dplyr)
Module = filter(Module, Module !='grey')

Module_colour_Gene = data.frame(Symbol="0",type="0",Module="0")
for (i in 1:nrow(Module)){
	data_1 = expredata_node[expredata_node$type == Module[i,1],]
	data_2 = data.frame(data_1,Module=paste("m",i,sep=""))
	Module_colour_Gene = rbind(Module_colour_Gene,data_2)
}
Module_colour_Gene = Module_colour_Gene[-1,]
write.table(Module_colour_Gene,"Module_colour_Gene.txt",row.names = F,quote = F,sep="\t", col.names = T)

Module_Gene = Module_colour_Gene[,c(3,1)]
write.table(Module_Gene,"Module_Gene.txt",row.names = F,quote = F,sep="\t", col.names = T)




## 寻找每个模块的关键基因
Module_colour = unique(Module_colour_Gene[,2:3])
HubGenes <- data.frame(colour=rownames(data.frame(HubGenes=chooseTopHubInEachModule(dataExpr,moduleColors))),HubGenes=chooseTopHubInEachModule(dataExpr,moduleColors))
HubGenes = merge(HubGenes, Module_colour, by.x = "colour",by.y = "type",all=FALSE)
write.table(HubGenes,"HubGenes.txt",row.names = F,quote = F,sep="\t", col.names = T)



### 建立Crosstalk文本Module_crosstalk.txt
Module <- data.frame(unique(Module_Gene[,1]))#对模块列去重
dir.create("Crosstalk_text")
for(i in 1:nrow(Module)){
	a = Module_Gene[Module_Gene$Module == Module[i,1],]
	b = data.frame((a[,2]))
	d = t(b)
	write.table(d,paste("Crosstalk_text/",Module[i,1], '.txt', sep=""),row.names = F,quote = F,sep="\t")#生成txt文件	
}

 
 

 
 
 
 



#################################
########    关联表型     ########
#################################

### 加入关联表型数据
trait <- "traitData.txt"
# 读入表型数据，不是必须的
if(trait != "") {
  traitData <- read.table(file=trait, sep='\t', header=T, row.names=1,
                          check.names=FALSE, comment='',quote="")
  sampleName = rownames(dataExpr)
  traitData = traitData[match(sampleName, rownames(traitData)), ]
}

### 模块与表型数据关联
if (corType=="pearsoon") {
  modTraitCor = cor(MEs_col, traitData, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, traitData, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}
## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.




### 绘制模块之间相关性图
# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
library(stringr)
MEs = net$MEs

### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)

# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", 
                      marDendro = c(3,3,2,4),
                      marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                      xLabelsAngle = 90)			
# 保存图片 Eigengene adjacency heatmap



## 如果有表型数据(traitData)，也可以跟ME数据放一起，一起出图
MEs_colpheno = orderMEs(cbind(MEs_col, traitData))
plotEigengeneNetworks(MEs_colpheno, "Eigengene adjacency heatmap", 
                     marDendro = c(3,3,2,4),
                     marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                     xLabelsAngle = 90)					
# 保存图片 Eigengene adjacency heatmap




# signif表示保留几位小数
textMatrix = paste(signif(modTraitCor, 2), "\n(", signif(modTraitP, 1), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
labeledHeatmap(Matrix = modTraitCor, xLabels = colnames(traitData), 
               yLabels = colnames(MEs_col), 
               cex.lab = 0.5, 
               ySymbols = colnames(MEs_col), colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, setStdMargins = FALSE, 
               cex.text = 0.5, zlim = c(-1,1),
               main = paste("Module-trait relationships"))	
# 保存图片 Module-trait relationships




## 从上图可以看到MEmagenta与Insulin_ug_l相关
## 模块内基因与表型数据关联
# 性状跟模块虽然求出了相关性，可以挑选最相关的那些模块来分析，
# 但是模块本身仍然包含非常多的基因，还需进一步的寻找最重要的基因。
# 所有的模块都可以跟基因算出相关系数，所有的连续型性状也可以跟基因的表达
# 值算出相关系数。
# 如果跟性状显著相关基因也跟某个模块显著相关，那么这些基因可能就非常重要。

### 计算模块与基因的相关性矩阵
if (corType=="pearsoon") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
             as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}


# 计算性状与基因的相关性矩阵
## 只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵。
if (corType=="pearsoon") {
  geneTraitCor = as.data.frame(cor(dataExpr, traitData, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
             as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, traitData, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}
## Warning in bicor(x, y, use = use, ...): bicor: zero MAD in variable 'y'.
## Pearson correlation was used for individual columns with zero (or missing)
## MAD.


# 最后把两个相关性矩阵联合起来,指定感兴趣模块进行分析
table(mergedColors)
module = "yellow"
pheno = "WT"
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(traitData))
# 获取模块内的基因
moduleGenes = moduleColors == module

sizeGrWindow(7, 7)
par(mfrow = c(1,1))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)		
# 保存图片				   
