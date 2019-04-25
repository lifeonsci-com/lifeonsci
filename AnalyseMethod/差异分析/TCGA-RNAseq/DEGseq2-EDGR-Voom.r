## 


CHOL_SampleMatrix_sort <- read.table("CHOL_SampleMatrix_sort.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE, check.names = F)
rownames(CHOL_SampleMatrix_sort) = CHOL_SampleMatrix_sort[,1]



###  rlogTransformation DESeq2的标准化
library(DESeq2)
expr =CHOL_SampleMatrix_sort[,2:ncol(CHOL_SampleMatrix_sort)]
group_list=ifelse(substr(colnames(expr),14,15)=='01','tumor','normal')
exprSet=na.omit(expr)

colData <- data.frame(row.names=colnames(exprSet), 
 group_list=group_list)
						 
dds <- DESeqDataSetFromMatrix(countData = exprSet,
colData = colData,
design = ~ group_list)

dds <- DESeq(dds)
res <- results(dds, 
 contrast=c("group_list","tumor","normal"))
resOrdered <- res[order(res$padj),]
head(resOrdered)
DEG =as.data.frame(resOrdered)
DEG =data.frame(DEG=rownames(DEG),DEG)
DEG_sort_p <- na.omit(DEG[DEG$pvalue < 0.01 ,])

write.table(DEG_sort_p, file = "DEG_sort_p.txt", append = FALSE, quote = F, sep = "\t",row.names = F,col.names = TRUE);
write.table(DEG, file = "DegData_DESeq2.txt", append = FALSE, quote = F, sep = "\t",row.names = F,col.names = TRUE);

# 可以看到程序非常好用！
# 它只对RNA-seq的基因的reads的counts数进行分析，请不要用RPKM等经过了normlization的表达矩阵来分析。
rld <- rlogTransformation(dds)  ## 得到经过DESeq2软件normlization的表达矩阵！
# rlog() may take a long time with 50 or more samples,
# vst() is a much faster transformation
# 数据集小于30使用rlog(),大数据集vst()
exprSet_new=assay(rld)
par(cex = 0.7)
n.sample=ncol(exprSet)
if(n.sample>40) par(cex = 0.5)
cols <- rainbow(n.sample*1.2)
par(mfrow=c(2,2))
boxplot(exprSet, col = cols,main="expression value",las=2)
boxplot(exprSet_new, col = cols,main="expression value",las=2)
hist(exprSet)
hist(exprSet_new)













## calcNormFactors 标准化
library(edgeR)
expr =CHOL_SampleMatrix_sort[,2:ncol(CHOL_SampleMatrix_sort)]
group_list=ifelse(substr(colnames(expr),14,15)=='01','tumor','normal')
exprSet=na.omit(expr)

d <- DGEList(counts=exprSet,group=factor(group_list))
keep <- rowSums(cpm(d)>1) >= 2
table(keep)
d <- d[keep, , keep.lib.sizes=FALSE]
d$samples$lib.size <- colSums(d$counts)
d <- calcNormFactors(d)
d$samples
dge=d
design <- model.matrix(~0+factor(group_list))
rownames(design)<-colnames(dge)
colnames(design)<-levels(factor(group_list))
dge=d
dge <- estimateGLMCommonDisp(dge,design)
dge <- estimateGLMTrendedDisp(dge, design)
dge <- estimateGLMTagwiseDisp(dge, design)

fit <- glmFit(dge, design)
# https://www.biostars.org/p/110861/
lrt <- glmLRT(fit,contrast=c(-1,1)) 
nrDEG=topTags(lrt, n=nrow(dge))
nrDEG=as.data.frame(nrDEG)
head(nrDEG)
edgeR_DEG =nrDEG 
edgeR_DEG =data.frame(DEG=rownames(edgeR_DEG),edgeR_DEG)
DEG_sort_p <- edgeR_DEG[edgeR_DEG$PValue < 0.01 ,]

write.table(DEG_sort_p, file = "DEG_sort_p.txt", append = FALSE, quote = F, sep = "\t",row.names = F,col.names = TRUE);
write.table(edgeR_DEG, file = "DegData_edgeR.txt", append = FALSE, quote = F, sep = "\t",row.names = F,col.names = TRUE);












##
library(limma)
library(edgeR)
expr =CHOL_SampleMatrix_sort[,2:ncol(CHOL_SampleMatrix_sort)]
group_list=ifelse(substr(colnames(expr),14,15)=='01','tumor','normal')
exprSet=na.omit(expr)

design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)
design

dge <- DGEList(counts=exprSet)
# dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

v <- voom(dge,design,plot=TRUE, normalize="quantile")#标准化
# exprSet_new=v$E #获取标准化的表达谱
fit <- lmFit(v, design)

group_list
cont.matrix=makeContrasts(contrasts=c('tumor-normal'),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

tempOutput = topTable(fit2, coef='tumor-normal', n=Inf)
DEG_limma_voom = na.omit(tempOutput)
head(DEG_limma_voom)
DEG_limma_voom = data.frame(DEG=rownames(DEG_limma_voom),DEG_limma_voom)
DEG_sort_p <- DEG_limma_voom[DEG_limma_voom$P.Value < 0.01 ,]

write.table(DEG_sort_p, file = "DEG_sort_p.txt", append = FALSE, quote = F, sep = "\t",row.names = F,col.names = TRUE);
write.table(DEG_limma_voom, file = "DegData_limma.txt", append = FALSE, quote = F, sep = "\t",row.names = F,col.names = TRUE);





















