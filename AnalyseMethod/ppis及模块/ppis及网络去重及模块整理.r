### 5.PPIs
setwd(paste(oldwd,"6. ppis\\3.sig",sep="/"))
#读入蛋白网路及筛选得分
#自己写

#plotgene为所需要构造ppis的基因文本
sig_ppis=String_human_PPIs_Symbol[(String_human_PPIs_Symbol$protein1 %in% plotgene$Gene) & (String_human_PPIs_Symbol$protein2 %in% plotgene$Gene),]
#通过table查看网络是否有误
table(sig_ppis$protein1 %in% plotgene$Gene)
table(sig_ppis$protein2 %in% plotgene$Gene)

#网路去重
data2=sig_ppis
colnames(data2)[1:2]=c("Gene1","Gene2")
data2$Gene1 <-as.character(data2[,1])
data2$Gene2 <-as.character(data2[,2])

if (nrow(data2) > 0 ){
	table(data2$Gene1>data2$Gene2)
	z=data2[data2$Gene1>data2$Gene2,]
	zz=data2[data2$Gene1<data2$Gene2,]
	zzz=zz[,c(2,1,3)]
	colnames(zzz)=colnames(z)
	zzzz=unique(rbind(z,zzz))
	table(zzzz$Gene1>zzzz$Gene2)
	network=zzzz
}
write.table(network,"sig_ppis.txt",row.names = F,quote = F,sep="\t")



## pre_Module_Gene是通过算法得到的基因集合,注意读进来前做筛选pval这一步，即算法得到结果的基因集合需要至少需要p<0.05 
pre_Module_Gene <- read.table("keep_ppis.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = F, fill=TRUE)
Module_Gene = data.frame(Module = 0,Symbol = 0)

for(i in 1:nrow(pre_Module_Gene)){
	gene <- pre_Module_Gene[i,]
	gene_1 = data.frame(t(gene))
	colnames(gene_1)="Symbol"
	gene_2 = data.frame(gene_1[gene_1$"Symbol" != "",])
	colnames(gene_2)=colnames(gene_1)
	gene_3 = data.frame(Module=paste("m",i,sep=""),gene_2)
	Module_Gene = rbind(Module_Gene,gene_3)
}
Module_Gene = Module_Gene[-1,]
write.table(Module_Gene,"Module_Gene.txt",row.names = F,quote = F,sep="\t")
