setwd("D:\\projects\\xiaohuadao2\\Enrichment2\\gastric_cancer")
bar <- read.table("bar.txt",sep = "\t",stringsAsFactors = F,header = T,fill=T,quote = "")
# ggplot(bar,aes(x=interaction(ID),y=Count))+geom_bar(stat="identity")+geom_text(aes(label=Count), vjust=-0.2)#显示在上面
ggplot(bar,aes(x=ID,y=Count))+geom_bar(stat="identity")+geom_text(aes(label=Count), vjust=-0.2)#显示在上面
