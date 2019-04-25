setwd('E:/project')
Expdata <- read.table("Expdata.txt",sep = "\t",stringsAsFactors = F,header = T,fill = TRUE)
GPL3802 <- read.table("GPL3802.txt",sep = "\t",stringsAsFactors = F,header = T,fill = TRUE)
Expdata_Result <- merge(GPL3802,Expdata,by.x = "ID",by.y = "ID",all=FALSE)
write.table(Expdata_Result,"Expdata_Result.txt",row.names = F,quote = F,sep="\t")#生成txt文件

Expdata_Result1 <- Expdata_Result[,-1]
write.table(Expdata_Result1,"Expdata_Result1.txt",row.names = F,quote = F,sep="\t")#生成txt文件

Expdata_Result2 <- read.table("Expdata_Result2.txt",sep = "\t",stringsAsFactors = F,header = T,fill = TRUE)

Deg_exprData <- Expdata_Result2
rownames(Deg_exprData) <- Deg_exprData[,1]
Deg_exprData <- Deg_exprData[,-1]
temp <- as.matrix(Deg_exprData) 


hv <- heatmap.2(x, col=cm.colors(255), scale="column",margin=c(5, 10),
xlab="Sample", ylab= "Gene",labRow=F,
main="heatmap(Sample_Gene)",
tracecol="green", density="density", colCol=cc,
srtCol=45, adjCol=c(0.5,1))
