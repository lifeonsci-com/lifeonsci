## 富集
# source("http://bioconductor.org/biocLite.R")
# biocLite("clusterProfiler")
# biocLite("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
Enrichment <- read.table("Enrichment.txt",sep = "\t",stringsAsFactors = F,header = T)
Module <- data.frame(unique(Enrichment[,1]))#对模块列去重

for (i in 1:nrow(Module)){
  Enrich_Source <- Enrichment[Enrichment$Module == Module[i,1],]
  gene <- as.character(Enrich_Source[,2])

  ego_CC <- enrichGO(gene = gene,
                     OrgDb=org.Hs.eg.db,
                     ont = "CC",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
  CC = as.data.frame(ego_CC @ result)
  write.table(CC, paste('GO_CC_', Module[i,1], '.txt', sep=""),row.names = F,quote = F,sep="\t")

  ego_BP <- enrichGO(gene = gene,
                     OrgDb=org.Hs.eg.db,
                     ont = "BP",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
  BP = as.data.frame(ego_BP @ result)
  write.table(BP, paste('GO_BP_', Module[i,1], '.txt', sep=""),row.names = F,quote = F,sep="\t")
  
  ego_MF <- enrichGO(gene = gene,
                     OrgDb=org.Hs.eg.db,
                     ont = "MF",
                     pAdjustMethod = "BH",
                     minGSSize = 1,
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05,
                     readable = TRUE)
  MF = as.data.frame(ego_MF @ result)
  write.table(MF, paste('GO_MF_',Module[i,1], '.txt', sep=""),row.names = F,quote = F,sep="\t")

# barplot(ego_CC, showCategory=15,title="EnrichmentGO_CC")#条状图，按p从小到大排的

# dotplot(ego_BP,title="EnrichmentGO_BP_dot",showCategory = 20)#点图，按富集的数从大到小的

  KEGG <- enrichKEGG(gene = gene, organism = "hsa", keyType = "kegg",
  pvalueCutoff = 0.05, pAdjustMethod = "BH", minGSSize = 1,
  qvalueCutoff = 0.05, use_internal_data = FALSE)
  KEGG_Pathway <- as.data.frame(KEGG @ result)
  write.table(KEGG_Pathway , paste('KEGG_Pathway_', Module[i,1], '.txt', sep=""),row.names = F,quote = F,sep="\t")
}

## 获取当前路径的所有文件
dir = data.frame(list.files(getwd()))
colnames(dir) = "Module"
dir = data.frame(dir[grep("CC|MF|BP|Pathway",dir$"Module"),]) 
colnames(dir) = "Module"

library(splitstackshape)
dir3 = concat.split.multiple(dir, "Module", ".")
dir3 = data.frame(dir3)

data_all <- data.frame("Module" = "0","ID"="0","Description"="0","GeneRatio"="0","BgRatio"="0","pvalue"="0","p.adjust"="0","qvalue"="0","geneID"="0","Count"="0",stringsAsFactors = F)

for (i in 1:nrow(dir3)){
    data1 <- read.table(paste(dir3[i,1],".txt",sep=""), sep = "\t",stringsAsFactors = F, header = T, fill=TRUE, quote = "")
    if (nrow(data1) > 0 ){
        Module = data.frame(Module=dir3[i,1])
        data2 = cbind(Module=dir3[i,1],data1)
        data_all <- rbind(data_all,data2,stringsAsFactors = F)
    }
}
Enrichment_all = data_all[-1,]
