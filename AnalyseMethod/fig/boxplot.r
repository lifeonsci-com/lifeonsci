### ICC

setwd("F:\\Projects\\ICC\\GSE45001")
GPL14550 <- read.table("GPL14550.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE)
GSE45001 <- read.table("GSE45001.txt",sep = "\t",comment.char = "!", stringsAsFactors = F,header = T, fill=TRUE)
GSM <- read.table("GSM.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE)

Zhushi = GPL14550[,c(1,7)]
Exp = merge(Zhushi,GSE45001,by.x = "ID",by.y = "ID_REF",all=FALSE)

Exp_1 = Exp[,-1]
Exp_2 = Exp_1[Exp_1$GENE_SYMBOL != "",]

exprData_result = Exp_2
exprData_result[,1] <- toupper(exprData_result[,1])
colnames(exprData_result)[1] <- "Symbol"
meanfun <- function(x) {
    x1 <- data.frame(unique(x[,1]))
	colnames(x1) <- c("Symbol")
    for (i in 2:ncol(x)){
        x2 <- data.frame(tapply(x[,i],x[,1],mean))
		x2[,2] <- rownames(x2)
		colnames(x2) <- c(colnames(x)[i], "Symbol")
        x1 <- merge(x1,x2,by.x = "Symbol",by.y = "Symbol",all=FALSE)
    }
    return(x1)
}
exprData_result <- meanfun(exprData_result)
write.table(exprData_result, 'exprData_result.txt',sep = '\t', quote=F, row.names=F)


PPM1D = exprData_result[exprData_result$Symbol == "PPM1D",]
GSM_Con = GSM[GSM$Type == "ICC_Non_Tumoral_Stroma",]
GSM_Case = GSM[GSM$Type == "ICC_Tumoral_Stroma",]

Exp_PPM1D_Con = data.frame(Con = " ")
write.table(Exp_PPM1D_Con, 'Exp_PPM1D_Con.txt',sep = '\t', quote=F, row.names=F)
for (i in 1:nrow(GSM_Con)){
  for (j in 1:ncol(PPM1D)){
    if (GSM_Con[i,1] == colnames(PPM1D)[j]){
	  Con = data.frame(PPM1D[1,j])
	  colnames(Con) = "Con"
	  write.table(Con, 'Exp_PPM1D_Con.txt', append=T, sep = '\t', quote=F, row.names=F, col.names=F)
	}
  }
}

Exp_PPM1D_Case = data.frame(Case = " ")
write.table(Exp_PPM1D_Case, 'Exp_PPM1D_Case.txt',sep = '\t', quote=F, row.names=F)
for (i in 1:nrow(GSM_Case)){
  for (j in 1:ncol(PPM1D)){
    if (GSM_Case[i,1] == colnames(PPM1D)[j]){
	  Case = data.frame(PPM1D[1,j])
	  colnames(Case) = "Case"
	  write.table(Case, 'Exp_PPM1D_Case.txt', append=T, sep = '\t', quote=F, row.names=F, col.names=F)
	}
  }
}

Exp_PPM1D_Con<- read.table("Exp_PPM1D_Con.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE)
Exp_PPM1D_Case <- read.table("Exp_PPM1D_Case.txt",sep = "\t",comment.char = "#", stringsAsFactors = F,header = T, fill=TRUE)
Exp_PPM1D = cbind(Exp_PPM1D_Con, Exp_PPM1D_Case)
Exp_PPM1D = Exp_PPM1D[-1,]

# 去除离群值(均值加减2倍标准差）
Ave_Con = mean(Exp_PPM1D$Con)
Sd_Con = sd(Exp_PPM1D$Con)
Max_Con = Ave_Con + 2 * Sd_Con
Min_Con = Ave_Con - 2 * Sd_Con
Con = data.frame(Exp_PPM1D[Exp_PPM1D$Con < Max_Con & Exp_PPM1D$Con > Min_Con,])
Con_1 = data.frame(Con[,1])
colnames(Con_1) = "Con"
Con_2 = data.frame(Group = "ICC_Non_Tumoral_Stroma", PPM1D = Con_1[,1])

Ave_Case = mean(Exp_PPM1D$Case)
Sd_Case = sd(Exp_PPM1D$Case)
Max_Case = Ave_Case+ 2 * Sd_Case
Min_Case = Ave_Case - 2 * Sd_Case
Case = data.frame(Exp_PPM1D[Exp_PPM1D$Case < Max_Case & Exp_PPM1D$Case > Min_Case,])
Case_1 = data.frame(Case[,1])
colnames(Case_1) = "Case"
Case_2 = data.frame(Group = "ICC_Tumoral_Stroma", PPM1D = Case_1[,1])

Exp_PPM1D_1 = rbind(Con_2, Case_2)

# 不去离群值
# colnames(Exp_PPM1D) = c("ICC_Non_Tumoral_Stroma", "ICC_Tumoral_Stroma")
# Exp_PPM1D_1 <- melt(Exp_PPM1D)
# colnames(Exp_PPM1D_1) = c("Group", "PPM1D")


# 差异检验
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("limma", version = "3.8")

library(limma)

# 把ExprData_result.txt编辑成rep("normal",5),rep("asthma",5))形式
# Expdata = data.frame(t(Exp_PPM1D_1)) # 去除离群值

Expdata = PPM1D


# 开始差异分析
rownames(Expdata)<-Expdata[,1]
Expdata <-Expdata[,-1]
Expdata <-as.matrix(Expdata)

Expdata_1 = data.frame(t(Expdata))
Expdata_2 = data.frame(GSM = rownames(Expdata_1),Expdata_1)
Expdata_3 = merge(Expdata_2,GSM,by.x = "GSM",by.y = "GSM",all=FALSE)
samps = factor(Expdata_3[,3])
design <- model.matrix(~0+samps)
design = design[,c(2,1)]
colnames(design) <- c("Case","Con")
fit <- lmFit(Expdata,design)
cont.matrix<-makeContrasts(Con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
final<-topTable(fit2, coef=1, number=dim(Expdata)[1], adjust.method="BH", sort.by="B", resort.by="M")
DEG = data.frame(DEG=rownames(final),final)
DEG_1 = data.frame(Gene = "PPM1D", DEG)
# DEG_sort <- DEG[DEG$adj.P.Val < 0.01 ,]
write.table(DEG_1,"DegData_limma.txt",quote=FALSE,sep="\t")

titile_plot = paste("Violin Plot p = ",DEG[1,5],", logFC = ",DEG[1,2])

# 箱式图
# Group和PPM1D为矩阵melt后的两列的名字，内部变量, Group代表了点线的属性，PPM1D代表对应的值。
p <- ggplot(Exp_PPM1D_1, aes(x=Group, y=PPM1D),color=Group) + 
geom_boxplot() + 
theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) +
theme(legend.position="none")
p


# 颜色箱式图
p <- ggplot(Exp_PPM1D_1, aes(x=Group, y=PPM1D),color=Group) + 
geom_boxplot(aes(fill=factor(Group))) + 
theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) +
theme(legend.position="none")
p
# 散点箱式图
p = ggplot(Exp_PPM1D_1,aes(x=Group, y=PPM1D,fill=Group))+geom_boxplot()+geom_jitter()+ 
theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) + labs(x="Group",y="PPM1D exp",title = "BoxPlot p = 0.4532922, logFC = -0.217")
p

# 小提琴图
p <- ggplot(Exp_PPM1D_1, aes(x=Group, y=PPM1D),color=Group) + 
geom_violin(aes(fill=factor(Group))) + 
theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) +
theme(legend.position="none")
p
# 散点小提琴图
p = ggplot(Exp_PPM1D_1,aes(x=Group, y=PPM1D,fill=Group))+geom_violin()+geom_jitter()+ 
theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) + labs(x="Group",y="PPM1D exp",title = "Violin Plot p = 0.4532922, logFC = -0.217")
p
