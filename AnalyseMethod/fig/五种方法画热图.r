### 五种方法画热图
# 读取数据并整理
Deg_exprData <- read.table("Deg_exprData.txt",sep = "\t", stringsAsFactors = F,header = T)
rownames(Deg_exprData) <- Deg_exprData[,1]
Deg_exprData <- Deg_exprData[,-1]
temp <- as.matrix(Deg_exprData) #矩阵化
#1、heatmap函数
heatmap(temp,col=colorRampPalette(c("blue","white","red"))(30),ColSideColors=colorRampPalette(c("blue","white","red"))(30),Colv=NA,cexRow=0.8,cexCol=1.2)
# (30)指的是x轴（即样本）有30个

# Col指定热图所用颜色：
# ColSideColors/RowSideColors代表列、行边是否显示颜色bar
# ColV/RowV表示是否按照列/行聚类，默认均为真值
# cexCol/cexRow分别表示列/行标签字体大小
# scale=c("row","column","none")#设置是否归一化
# margins=c(5,5)#设置热图下方及右方宽度
# 若热图中对于样本有多种分组，需要在行/列帝显示两行或多行颜色bar，可选择heatmap.plus或pheatmap包

#2、pheatmap
library(pheatmap)
# cellwidth、cellheigh t#小方格宽度、高度
# scale="none" ##是否归一化
# cluster_rows、cluster_cols #是否按行、列聚类
# treeheight_row、treeheight_col #横向、纵向树高度
# legend、annotation（设置分组）等高级选项见后面
# display_numbers=TRUE #在小方格中显示数字
###注意：以下一行代码生成一张图
pheatmap(temp) #无参数
pheatmap(temp,treeheight_row=120,treeheight_col=20) #设置col、row方向的树高
pheatmap(temp,treeheight_row=120,treeheight_col=20,cluster_cols=FALSE) # 单一方向聚类
pheatmap(temp,treeheight_row=120,treeheight_col=20,cluster_cols=FALSE,color=colorRampPalette(c("green","black","red"))(1000))#更改颜色
pheatmap(temp,treeheight_row=120,treeheight_col=20,cluster_cols=FALSE,color=colorRampPalette(c("green","black","red"))(5),border_color=NA,fontsize=10,fontsize_row=8,fontsize_col=16,filename="fig1.pdf",width=6,height=14) #去掉方格边框、调节字体大小、保存
# treeheight_col、treeheight_row分别为纵向、横向树形高度
# border_color方格边框颜色
# fontsize为所有字体大小，fontsize_row、fontsize_col分别为row、col方向标签字体大小
pheatmap(temp,cluster_cols=FALSE,cluster_rows=TRUE,legend=TRUE,color=colorRampPalette(c("green","black","red"))(1000),border_color=FALSE,fontsize=10,fontsize_row=12,fontsize_col=12,annotation=anno,annotation_legend=TRUE,annotation_colors=anno_colors) #设置多个分组
# annotation是一个dataframe，它的每行代表一个列（样本）的信息，包括每一列（样本）所属的group名，它的colnames代表分组名，rownames为列（样本）名；annotation_colors是list，每个列表元素代表一个级别，包括分组中各组名称以及对应的颜色值。
pheatmap(temp,cluster_col=FALSE,cluster_rows=FALSE,legend=FALSE,color=colorRampPalette(c("green","black","red"))(1000),border_color=TRUE,fontsize=10,fontsize_row=12,fontsize_col=12,width=6,height=14) #不聚类
pheatmap(temp,cluster_col=FALSE,cluster_rows=FALSE,legend=FALSE,color=colorRampPalette(c("blue","white","red"))(1000),border_color=TRUE,fontsize=10,fontsize_row=12,fontsize_col=12,display_numbers=TRUE,number_format="%.2f",width=6,height=14) #不聚类、填充数字

#3、ggplot2包
library("ggplot2")
library("reshape2")
##首先分别对数据temp在行/列方向上聚类，并保存聚类后的行列顺序
hc=hclust(dist(temp))
row_order=hc$order
temp1=temp[row_order,] ###是否对行列进行聚类
temp1=melt(temp1) ###数据变换，从matrix到ggplot可以识别的类型
p<-ggplot(temp1,aes(x=Var2,y=Var1,fill=value))+
xlab("")+ylab("")+labs(title="")+geom_tile(colour="white",size=0)+scale_fill_grandient(low="green",high="red")+
geom_text(aes(label=round(value,2)),angel=45,size=3)###加数字
print(p)

#4、gplots包（heatmap.2）
#scale按行或列均一化（"col"或"row"）
#Rowv、Colv是否对行、列聚类
#dendrogram是否绘制树形（none、both、row、col）
#margins设置下方、右方label宽度(eg:margins=c(5,7))
#ColSideColors，RowSideColors设置列、行分组颜色（颜色分别对应各列、各行）
library("gplots")
group=colorRampPalette(c("green","red"))(12)
#12种表示分组关系的颜色
heatmap.2(temp,col=redgreen,Colv=FALSE,ColSideColors=group,key=TRUE,symkey=FALSE,density.info="none",Rowv = F,trace="none")

#5、lattice包
# library("lattice")
# library("latticeExtra")
# colorkey=list(space="left",width=1.5) #设置颜色条的宽度和位置（top、left、right、bottom，但不可以与树形放置在同一侧）
# legend=list(...)#在图上的顶部或右部添加聚类树（仅可以在右方或上方加），并定义树形状（type参数指定：三角形triangle或矩形<默认>）
# levelplot(t(data),aspect="fill",...) ###data为matrix

library("lattice")
library(latticeExtra)
hc=hclust(dist(temp)) ###按行聚类
dd.row=as.dendrogram(hc)###保存行聚类树形
row.ord=order.dendrogram(dd.row) ###保存行聚类顺序
hc=hclust(dist(t(temp))) ###按列聚类
dd.col=as.dendrogram(hc) ###保存列聚类树形
col.rod=order.dendrogram(dd.col) ###保存列聚类顺序
temp1=temp[row.ord,] ###只对行聚类（是否对行、列聚类）
levelplot(t(temp1),aspect="fill",colorkey=list(space="left",width=1.5),xlab="",ylab="",legend=list(left=list(fun=dendrogramGrob,args=list(x=dd.row,rod=row.ord,side='left',size=5)),
scales=list(x=list(rot=60))))###x轴标签旋转60度