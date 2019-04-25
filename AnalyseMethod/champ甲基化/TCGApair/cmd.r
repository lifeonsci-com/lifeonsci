
library("impute")
library("ChAMP")
library(pheatmap)

RawNARemove<-function(data,missratio=0.3){
  threshold<-(missratio)*dim(data)[2]
  NaRaw<-which(apply(data,1,function(x) sum(is.na(x))>threshold))
  zero<-which(apply(data,1,function(x) all(x==0))==T)
  NaRAW<-c(NaRaw,zero)
  if(length(NaRAW)>0){
    data1<-data[-NaRAW,]
  }else{
    data1<-data;
  }
  data1
}

data <- read.table('LUAD_JHU_USC__HumanMethylation450_beta.txt',stringsAsFactors=F,sep="\t",check.names=F)
methdata=RawNARemove(data)
newmethdata<-impute.knn(methdata)$data
rownames(newmethdata)=rownames(methdata)
colnames(newmethdata)=colnames(methdata)

class <- c('01','11')
data.t <- newmethdata[,substr(colnames(newmethdata),14,15) %in% class]

sample_id <- substr(colnames(data.t) ,1 ,12)
data.p<- data.t[,which(sample_id %in% sample_id[duplicated(sample_id)])]
tumors <- as.data.frame(data.p[,which(substr(colnames(data.p) ,14 ,15) %in% class[1])])
normals <- as.data.frame(data.p[,paste(substr(colnames(tumors),1,12),class[2],sep="-")])
count.table <- cbind(tumors , normals)

rownames(count.table) <- gsub(":.*",'',rownames(count.table))
groups=as.factor(c(rep("Tumor", length(colnames(tumors))),rep("Normal", length(colnames(normals)))))
myDMP <- champ.DMP(beta = count.table ,pheno=groups)  ## same result with limma
myDMR <- champ.DMR(beta = as.matrix(count.table),pheno=groups,method="DMRcate")
if(!dir.exists('methy')){
	dir.create('methy')
}
#colnames(count.table) <- gsub(".",'-',colnames(count.table),fixed=T)
write.table(file="methy/paired_methylation_betavalue.txt",data.frame(Probe=rownames(count.table),count.table, check.names=F,stringsAsFactors=F),row.names=F,sep="\t",quote=F)
write.table(file="methy/paired_methylation_DMP.txt",data.frame(Probe=rownames(myDMP$Tumor_to_Normal),myDMP$Tumor_to_Normal, check.names=F,stringsAsFactors=F),row.names=F,sep="\t",quote=F)
write.table(file="methy/paired_methylation_DMR.txt",data.frame(DMR=rownames(myDMR$DMRcateDMR), myDMR$DMRcateDMR, check.names=F,stringsAsFactors=F),row.names=F,sep="\t",quote=F)

data <- myDMP$Tumor_to_Normal
data.f <- data[which(data$cgi == 'island' & !data$feature %in% c('Body','IGR',"3'UTR")),]
beta.f <- count.table[rownames(data.f),] 
write.table(file="methy/paired_methylation_DMP_filter.txt",data.frame(Probe=rownames(data.f), data.f, check.names=F,stringsAsFactors=F),row.names=F,sep="\t",quote=F)
write.table(file="methy/paired_methylation_betavalue_filter.txt",data.frame(Probe=rownames(beta.f), beta.f, check.names=F,stringsAsFactors=F),row.names=F,sep="\t",quote=F)

annotation_col = data.frame(CellType = factor(rep(c("Tumor", "Normal"), each=29)))
rownames(annotation_col) = colnames(beta.f)
tiff('methy/paired_methlation_beta_heatmap.tiff',width = 600, height = 600)
pheatmap(as.matrix(beta.f), cluster_cols =F, show_rownames=F, show_colnames=F, annotation_col = annotation_col)
dev.off()

