## DEG (P.Value < 0.05)
library(limma)
# 输入文本是case-control的表达谱（前面是case，后面是control）
exp = DEG_EXP
for (i in 2:ncol(exp)){
	exp[,i] = as.numeric(exp[,i])
}
# 开始差异分析
rownames(exp)<-exp[,1]
exp<-exp[,-1]
# exp = exp[,1:12]
exp<-as.matrix(exp)

samps<-factor(c(rep("Case",5),rep("Control_sample",5))
design <- model.matrix(~0+samps);
colnames(design) <- c("Case","Control_sample")
rownames(design)=colnames(exp)
design #检查design正误
fit <- lmFit(exp,design)
cont.matrix<-makeContrasts(Case-Control_sample,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
final<-topTable(fit2, coef=1, number=dim(exp)[1], adjust.method="BH", sort.by="B", resort.by="M")

DEG = data.frame(DEG=rownames(final),final)
DEG_sort_p <- DEG[DEG$P.Value < 0.05 ,]
