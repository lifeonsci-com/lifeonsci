rm(list=ls())
gc()

setwd("D:/Corporation/SixArticle/Pancreatic_Cancer/Tumor_Exp/GSE62452")
som_svd_res=read.csv("Som_Svd_output.csv",header=T)
colnames(som_svd_res)


plotData=som_svd_res[,5:9]
colnames(plotData)
head(rownames(plotData))

plotData=apply(plotData, 2, as.numeric)
sum(is.na(plotData))
head(plotData)
library(fpc)

set.seed(252964)

norm_plotData=scale(t(plotData))
norm_plotData=t(norm_plotData)
sum(is.na(norm_plotData))
table(som_svd_res[,3])
cluster=som_svd_res[,3]
#kmeans <- kmeans(norm_plotData, 6) 

#plotcluster(norm_plotData, kmeans$cluster) 
#kmeans$cluster
#cluster=kmeans$cluster
#cluster[which(cluster==3)]=1.5
#cluster
write.table(cbind(som_svd_res[,1:4],kmeans$cluster),"clusterInfor.txt",row.names = F,col.names = F,sep="\t",quote=F)
kmdata=norm_plotData[order(cluster),]
library(pheatmap)
pheatmap(kmdata, cellwidth = 56,color = colorRampPalette(c("green3","orange", "red"))(100),cluster_row=F,cluster_col=F, fontsize=9, fontsize_row=6)

library(gplots)
pheatmap(kmdata, cellwidth = 56,color = greenred(100),cluster_row=F,cluster_col=F, fontsize=9, fontsize_row=6)

#save.image("pltData.Rdata")
#rm(list=ls())
#gc()
c1=norm_plotData[which(cluster==1),]
c2=norm_plotData[which(cluster==2),]
c3=norm_plotData[which(cluster==3),]
c4=norm_plotData[which(cluster==4),]
c5=norm_plotData[which(cluster==5),]
c6=norm_plotData[which(cluster==6),]
x=1:5
c1median=apply(c1, 2, median)
stage=colnames(c1)
ef1=lm(c1median~x)
#b=-0.9314  k=0.2915
predict(ef1)
plot1=data.frame(median=c1median,stage=stage,predict=predict(ef1))
library(reshape2)
meltPlot1=melt(plot1)
library(ggplot2)
#ggplot(plot1,aes(x=stage,y=median,group=1))+geom_line(colour="#87CEFA",size=1)+geom_point(shape=20,colour="#87CEFA",size=5)
ggplot(meltPlot1,aes(x=stage,y=value,colour=variable,group=variable))+geom_line()+geom_point(aes(x=stage,y=value,colour=variable),shape=20,colour="#87CEFA",size=5)+scale_color_manual(values=c("#87CEFA","#e55c5c"))+theme_bw()

c2median=apply(c2, 2, median)
stage=colnames(c2)
ef2=lm(c2median~x)
#b=-1.3417  k=0.4612
predict(ef2)
plot2=data.frame(median=c2median,stage=stage,predict=predict(ef2))
library(reshape2)
meltPlot2=melt(plot2)
library(ggplot2)
#ggplot(plot1,aes(x=stage,y=median,group=1))+geom_line(colour="#87CEFA",size=1)+geom_point(shape=20,colour="#87CEFA",size=5)
ggplot(meltPlot2,aes(x=stage,y=value,colour=variable,group=variable))+geom_line()+geom_point(aes(x=stage,y=value,colour=variable),shape=20,colour="#87CEFA",size=5)+scale_color_manual(values=c("#87CEFA","#e55c5c"))+theme_bw()

c3median=apply(c3, 2, median)
stage=colnames(c3)
ef3=lm(c3median~x)
#b=0.17501  k=-0.09361
predict(ef3)
plot3=data.frame(median=c3median,stage=stage,predict=predict(ef3))
library(reshape2)
meltPlot3=melt(plot3)
library(ggplot2)
#ggplot(plot1,aes(x=stage,y=median,group=1))+geom_line(colour="#87CEFA",size=1)+geom_point(shape=20,colour="#87CEFA",size=5)
ggplot(meltPlot3,aes(x=stage,y=value,colour=variable,group=variable))+geom_line()+geom_point(aes(x=stage,y=value,colour=variable),shape=20,colour="#87CEFA",size=5)+scale_color_manual(values=c("#87CEFA","#e55c5c"))+theme_bw()

c4median=apply(c4, 2, median)
stage=colnames(c4)
ef4=lm(c4median~x)
#b=0.9572  k=-0.3275
predict(ef4)
plot4=data.frame(median=c4median,stage=stage,predict=predict(ef4))
library(reshape2)
meltPlot4=melt(plot4)
library(ggplot2)
#ggplot(plot1,aes(x=stage,y=median,group=1))+geom_line(colour="#87CEFA",size=1)+geom_point(shape=20,colour="#87CEFA",size=5)
ggplot(meltPlot4,aes(x=stage,y=value,colour=variable,group=variable))+geom_line()+geom_point(aes(x=stage,y=value,colour=variable),shape=20,colour="#87CEFA",size=5)+scale_color_manual(values=c("#87CEFA","#e55c5c"))+theme_bw()

c5median=apply(c5, 2, median)
stage=colnames(c5)
ef5=lm(c5median~x)
#b=1.616  k=-0.535
predict(ef5)
plot5=data.frame(median=c5median,stage=stage,predict=predict(ef5))
library(reshape2)
meltPlot5=melt(plot5)
library(ggplot2)
#ggplot(plot1,aes(x=stage,y=median,group=1))+geom_line(colour="#87CEFA",size=1)+geom_point(shape=20,colour="#87CEFA",size=5)
ggplot(meltPlot5,aes(x=stage,y=value,colour=variable,group=variable))+geom_line()+geom_point(aes(x=stage,y=value,colour=variable),shape=20,colour="#87CEFA",size=5)+scale_color_manual(values=c("#87CEFA","#e55c5c"))+theme_bw()

c6median=apply(c6, 2, median)
stage=colnames(c6)
ef6=lm(c6median~x)
#b=0.7102  k=-0.2341
predict(ef6)
plot6=data.frame(median=c6median,stage=stage,predict=predict(ef6))
library(reshape2)
meltPlot6=melt(plot6)
library(ggplot2)
#ggplot(plot1,aes(x=stage,y=median,group=1))+geom_line(colour="#87CEFA",size=1)+geom_point(shape=20,colour="#87CEFA",size=5)
ggplot(meltPlot6,aes(x=stage,y=value,colour=variable,group=variable))+geom_line()+geom_point(aes(x=stage,y=value,colour=variable),shape=20,colour="#87CEFA",size=5)+scale_color_manual(values=c("#87CEFA","#e55c5c"))+theme_bw()






x=1:6
lm(c6median~x)
#b=-1.495 k=0.423
#b=1.1206 k=-0.3293
#b=-0.3747 k=0.1024
#b=1.4782 k=-0.4246
#b=-1.2927 k=0.3633
#b=0.09361 k=-0.03205






1####################################################################
load("pltData.Rdata")
lineData=norm_plotData[which(kmeans$cluster==4),]
dim(lineData)
range(lineData)
x=1:6
plot(x,lineData[1,],type = "l",col=grey(0.6),lwd=1,ylim = c(-1.9,1.8))

for(i in 2:nrow(lineData)){
  lines(x,lineData[i,],type = "l",col=grey(0.6),lwd=1)
}
medianCol=apply(lineData,2,median)
lines(x,medianCol,type = "l",col=rgb(255,48,48,maxColorValue = 255),lwd=2)
ef=lm(medianCol~x)
lines(x,predict(ef),type = "l",col="blue",lwd=2)
lm(medianCol~x)
#####################################################################