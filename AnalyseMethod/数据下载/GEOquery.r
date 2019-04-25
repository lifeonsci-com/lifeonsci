### GEO数据下载
library(GEOquery) 
eSet <- getGEO("GSE42872",destdir = ".",AnnotGPL = F,getGPL = F) 
save(eSet,file = "GSE42872_eSet.Rdata")

### ID转换
ob <- eSet[[1]] 
exprSet <- exprs(ob) 
samples <- sampleNames(ob) 
pdata <- pData(ob) 
group_list <- as.character(pdata[,2]) 



### 其次，下载GPL注释信息：
## ①与上面类似，使用getGEO下载：
    library(GEOquery)
    gpl <- getGEO("GPL16570",destdir = ".") 
    probe2gene <- Table(gpl)[,c(1,11)]
    save(probe2gene,file = "probe2gene.Rdata")

## ②如果可以找到该平台所对应的R包(http://www.bio-info-trainee.com/1399.html)，则可以这样下载：
    library(BiocInstaller)
    BiocInstaller::biocLite("hugene10sttranscriptcluster.db") 
library(hugene10sttranscriptcluster.db) 



### hcluster图
colnames(new_exprSet) <- paste(group_list,1:ncol(new_exprSet)) 
nodepar <- list(lab.cex=0.6,pch=c(NA,19),
                cex=0.7,col="blue")
hc <- hclust(dist(t(new_exprSet)))
par(mar=c(5,5,5,10))
plot(as.dendrogram(hc),nodePar=nodepar,horiz=T)

### PCA图
library(ggfortify)
library(ggplot2)
df=as.data.frame(t(new_exprSet))
df$group=group_list
autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = 'group')


### 设置分组信息
group_list <- factor(c(rep("control",3), rep("experiment",3)))  


### 使用limma包进行差异基因的分析
suppressMessages(library(limma))
design <- model.matrix(~0+factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(new_exprSet)
contrast.matrix <- makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
##step1
fit <- lmFit(new_exprSet,design)
##step2
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) ## default no trend !!!
##step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput)
write.csv(nrDEG2,"limma_notrend.results.csv",quote = F)


### 画热图
## apply函数一般有三个参数，第一个参数代表矩阵对象，第二个参数代表要操作矩阵的维度，1表示对行进行处理，2表示对列进行处理。第三个参数就是处理数据的函数。apply会分别一行或一列处理该矩阵的数据。
library(pheatmap)
choose_gene=names(tail(sort(apply(new_exprSet,1,mad)),50))
choose_matrix=new_exprSet[choose_gene,]
choose_matrix=t(scale(t(choose_matrix)))
pheatmap(choose_matrix)








