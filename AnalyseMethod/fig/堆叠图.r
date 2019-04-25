setwd("D:\\projects\\liver2\\堆叠图")
Enrichment_Group <- read.table("Enrichment_Group.txt",sep = "\t",stringsAsFactors = F,header = T,fill=T,quote = "")
dr = Enrichment_Group
row.names(dr)=dr[,1]
dr = dr[,-1]
sum(dr)

##  绘图  ##


par(mfrow = c(1,1))
x <- barplot(t(dr),legend=rownames(t(dr)),col=c("black","purple","cyan","green","blue","yellow","orange","red"),xlab="Module")




