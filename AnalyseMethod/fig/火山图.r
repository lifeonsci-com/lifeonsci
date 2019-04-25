# install.packages("Cairo")
library(ggplot2)
library(Cairo)

setwd("G:\\") # 设置工作空间
# 读取文件，data为limma筛选差异后的表格
# data=final #继续limma工作data就是final
data <- read.table("DegData_limma.txt",sep = "\t", stringsAsFactors = F,header = T)

# 设置颜色域 
data$threshold <- as.factor(ifelse(data$P.Value < 0.05 & abs(data$logFC) >=1.5,ifelse(data$logFC > 1.5 ,'Up','Down'),'Not'))
# 说明：注意logFC的参数，正常用1.5，若实验倍数差异不明显，可使用0.5（后面参数同样修改）
##Construct the plot object
# with legend
# 生成文件
# Cairo(file="Volcan1.pdf", type="png",units="in",bg="white",width=5.5, height=5, pointsize=12, dpi=300)
ggplot(data=data, aes(x=logFC, y =-log10(P.Value), colour=threshold,fill=threshold)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_point(alpha=0.4, size=1.2) +
  xlim(c(-4, 4)) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=12),
        axis.text.y = element_text(face="bold",  color="black", size=12),
        axis.title.x = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(face="bold",color="black", size=12))+
  labs(x="log2 (fold change)",y="-log10 (p-value)",title="Volcano picture of DEG")
# dev.off()
 
 
 
# without legend（没有图例）
# Cairo(file="Volcan.png", type="png",units="in",bg="white",width=6, height=5, pointsize=12, dpi=300)
ggplot(data=data, aes(x=logFC, y =-log10(P.Value), colour=threshold,fill=threshold)) +
  scale_color_manual(values=c("blue", "grey","red"))+
  geom_point(alpha=0.4, size=1.2) +
  xlim(c(-4, 4)) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="grey",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.6)+
  theme(legend.position="none",
        panel.grid=element_blank(),
        # legend.title = element_blank(),
        # legend.text= element_text(face="bold", color="black",family = "Times", size=8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face="bold", color="black", size=12),
        axis.text.y = element_text(face="bold",  color="black", size=12),
        axis.title.x = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(face="bold",color="black", size=12))+
  labs(x="log2 (fold change)",y="-log10 (p-value)",title="Volcano picture of DEG")
# dev.off()

# 上调下调颜色一样
data$threshold = as.factor(data$P.Value < 0.05 & abs(data$logFC) >=1.5)
 
Cairo(file="volcan.png", 
      type="png",
      units="in",
      bg="white",
      width=8, 
      height=6, 
      pointsize=12, 
      dpi=300)
ggplot(data=data, aes(x=logFC, y =-log10(P.Value), colour=threshold)) +
  geom_point(alpha=0.4, size=1.75) +
  xlim(c(-4, 4)) +
  xlab("log2 fold change") + ylab("-log10 p-value") +
  ggtitle("Volcano picture of DEG")+
  theme_gray() +
  geom_vline(xintercept=c(-1.5,1.5),lty=4,col="grey",lwd=0.5)+ # 在x轴-1.5与1.5的位置画两根竖线
  geom_hline(yintercept = -log10(0.05),lty=4,col="grey",lwd=0.5)+ #在p value 0.05的位置画一根横线
  theme(legend.position="none")
 
dev.off()

# 将下调，上调的分别改变颜色体现
data$threshold <- as.factor(ifelse(data$P.Value < 0.05 & abs(data$logFC) >=1.5,ifelse(data$logFC > 1.5 ,'Up','Down'),'Not'))
scale_color_manual(values=c("blue", "grey","red"))
theme(legend.position="right")

