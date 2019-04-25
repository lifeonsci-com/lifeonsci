### 读入核心基因和差异文本
HubGenes <- read.table("C:\\Users\\waiting\\Desktop\\mo\\HubGenes.txt",sep = "\t",stringsAsFactors = F,header = T,fill = TRUE)
DEG_sort_0.05 <- read.table("F:\\Projects2018\\December\\Liver Cancer\\1.Difference analysis\\DEG_sort_0.05.txt",sep = "\t",stringsAsFactors = F,header = T,fill = TRUE)

volcanoData = merge(HubGenes, DEG_sort_0.05, by.x = "HubGenes",by.y = "DEG",all=FALSE)
# data$threshold <- as.factor(ifelse(data$P.Value < 0.05 & abs(data$logFC) >=0.5,ifelse(data$logFC > 0.5 ,'Up','Down'),'Not'))
direction=as.factor(ifelse(abs(volcanoData$log2FoldChange)>0.5,ifelse(volcanoData$log2FoldChange>0.5,"UP","DW"),"NO"))
volcanoData$direct=direction;






library(ggplot2)
library(ggplot2)
p <- ggplot(volcanoData, aes(x=Module, y=log2FoldChange))
p <- p + geom_point()
# 前面是给p不断添加图层的过程
# 单输入一个p是真正作图
# 前面有人说，上面都输完了，怎么没出图
# 就因为差了一个p
p



# xlim调整下X轴的区间使图对称
p <- ggplot(volcanoData,  aes(x=Module,  y=log2FoldChange)) +
     geom_point() +
     xlim(-2.5, 2.5)
p





## https://ggplot2.tidyverse.org/reference/geom_text.html
library(ggbeeswarm)
p <- ggplot(volcanoData, aes(x=log2FoldChange, y=Module)) + geom_quasirandom()
# 使用geom_text增加点的标记
# label表示标记哪一列的数值
# position_quasirandom 获取点偏移后的位置
# xjust调整对齐方式; hjust是水平的对齐方式，0为左，1为右，0.5居中，0-1之间可以取任意值。vjust是垂直对齐方式，0底对齐，1为顶对齐，0.5居中，0-1之间可以取任意值。
# check_overlap检查名字在图上是否重叠
p <- p + geom_text(aes(label=HubGenes),position=position_quasirandom(),hjust=0, check_overlap=T)
p



###  点的标签位置
# nudge_x, nudge_y	Horizontal and vertical adjustment to nudge labels by. Useful for offsetting text from points, particularly on discrete scales.
p <- ggplot(volcanoData, aes(x=log2FoldChange, y=Module)) + geom_quasirandom()
p <- p + geom_text(aes(label=HubGenes),hjust=0,nudge_x = 0.05, check_overlap=T)
p




## 完整版
p <- ggplot(volcanoData, aes(x=log2FoldChange, y=Module,color=direct)) + geom_quasirandom()+xlim(-2.5, 2.5)+
	geom_point(size=5,shape=20) + geom_vline(xintercept = -0.5,linetype="dashed",color="#9B9898")+
	geom_vline(xintercept = 0.5,linetype="dashed",color="#9B9898")+geom_vline(xintercept = 0,color="#9B9898")+
	scale_color_manual(values=c("#4388F7","#514F4F","#F72F8F"))

p <- p + geom_text(aes(label=HubGenes),hjust=0.5,xjust=0.5,nudge_x = 0.05, check_overlap=T)
p





