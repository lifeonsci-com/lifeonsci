## 构建表达谱
## CEL文件处理
source("http://bioconductor.org/biocLite.R")
cels = list.files("GSE64614_RAW/",pattern='[gz]')
## 安装R.utils并加载
# biocLite("R.utils")
# biocLite("R.oo")
# biocLite("R.methodsS3")
library(R.utils)
library(R.oo)
library(R.methodsS3)
library(hgu133plus2cdf)
sapply(paste('GSE64614_RAW',cels,sep='/'),gunzip)
celpath=paste(getwd(),'GSE64614_RAW',sep='/')
oldWD=setwd(celpath)
## 数据预处理
# biocLite("affy")
library(affy)
library(BiocGenerics)
library(parallel)
library(Biobase)
raw_data=ReadAffy()
setwd(oldWD)
# 删除GSE16515_1文件及其内容
unlink('GSE64614_RAW',recursive=TRUE)
# 用rma方法处理原始数据，其结果是经过对数变换的；也可以用mas5处理，其结果是原始信号强度
# mas_data = mas5(raw_data)
rma_data=rma(raw_data)
rma_exp=exprs(rma_data)
write.table(rma_exp, 'probeExpData/GSE64614.txt', sep='\t',row.names = TRUE, col.names = TRUE, quote=F)


## 注释
# biocLite('annotate')
# biocLite('hgu133plus2.db')
library(annotate)
library(AnnotationDbi)
library(stats4)
library(IRanges)
library(S4Vectors)
library(XML)
affydb=annPkgName(rma_data@annotation,type='db')
library(affydb,character.only=TRUE)
library(org.Hs.eg.db)
## 去除一对空的探针和一对多的探针
# raw_symbols = as.matrix(getSYMBOL(rownames(rma_exp),affydb))
raw_geneid = as.matrix(getEG(rownames(rma_exp),affydb))
colnames(raw_geneid) = c('geneid')
new_geneid = as.numeric(as.matrix(raw_geneid))
rma_exp2 = cbind(new_geneid, rma_exp)
# 去除NA值
geneid=na.omit(new_geneid)
# loc=match(rownames(rma_exp),rownames(geneid))
dim(rma_exp2)
rma_exp3 = na.omit(rma_exp2)[,-1]
dim(rma_exp3)
## 对多个探针对应一个基因的探针集取均值
geneid_factor=factor(geneid)
gene_exp_matrix=apply(rma_exp3, 2, function(x) tapply(x,geneid_factor,mean))
dim(gene_exp_matrix)
rownames(gene_exp_matrix) = levels(geneid_factor)
gene_exp_matrix_1 = data.frame(Entrez_ID = rownames(gene_exp_matrix),gene_exp_matrix)
write.table(gene_exp_matrix_1, 'normalize_data/GSE69438_1.txt',sep = '\t', quote=F, row.names=F)