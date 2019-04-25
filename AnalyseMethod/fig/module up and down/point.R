rm(list = ls())
gc()

setwd("E:\\CompanyMession\\others\\Figure")
list.files(getwd())

modules=read.table("Module_Gene.txt",header=T,sep="\t")
head(modules)

FC=read.table("DegData_limma.txt",header=T,sep="\t",row.names = NULL)
head(FC)

index=match(modules$Symbol,FC$row.names)
head(FC[index,2])

#runif(nrow(modules),min=-0.05,max = 0.05)
res=data.frame(module=modules$Module+runif(nrow(modules),min=-0.1,max = 0.1),logFC=FC[index,2])
head(res)
res[which(res$logFC<(-2.5)),2]=(-1.78)
sum(res$logFC<(-2.5))
library(ggplot2)
direction=as.factor(ifelse(abs(res$logFC)>0.5,ifelse(res$logFC>0.5,"UP","DW"),"NO"))
res$direct=direction;
ggplot(res,aes(x=logFC,y=module,color=direct))+geom_point(size=2,shape=20)+scale_y_continuous(breaks = seq(1,9,by=1))+scale_x_continuous(breaks = c(-1.5,-0.5,0,0.5,1,2.5))+geom_vline(xintercept = -0.5,linetype="dashed",color="#9B9898")+geom_vline(xintercept = 0,color="#9B9898")+geom_vline(xintercept = 0.5,linetype="dashed",color="#9B9898")+scale_color_manual(values=c("#4388F7","#514F4F","#F72F8F"))
