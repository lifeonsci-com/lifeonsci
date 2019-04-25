setwd('E:/HMU/project1/metaData/refData/Rtest')
library(affydb,character.only=TRUE)
raw_geneid = as.matrix(getEG(rownames(rma_exp),affydb))
entrezID = unique(unlist(raw_geneid))
entrezID=na.omit(entrezID)
dif = read.table('degData/pancreatic_deg.txt', header = FALSE)
dif = as.matrix(as.character(dif))

##GO富集分析
library(GOstats)
params = new('GOHyperGParams', geneIds = dif, universeGeneIds = entrezID, 
annotation = affydb, ontology = 'BP', pvalueCutoff = 0.001, conditional = FALSE,
testDirection = 'over')

hgOver = hyperGTest(params)
bp = summary(hgOver)
htmlReport(hgOver, file = 'pancreatic_go.html')
head(bp)

##KEGG富集分析
library(GeneAnswers)

##KEGG富集和GO富集统一的包
library(clusterProfiler)
entrezID = read.table('entrezID.txt')
entrezID = as.character(as.matrix(entrezID))

gene = c('6879','117143','112869','6878','93624','55689')
name = c('m63')

#GO富集
enrichGO_MF = enrichGO(gene = gene, 'org.Hs.eg.db', universe = entrezID, 
ont = 'MF', pvalueCutoff = 0.01, readable = TRUE)
enrichGO_CC = enrichGO(gene = gene, 'org.Hs.eg.db', universe = entrezID, 
ont = 'CC', pvalueCutoff = 0.01, readable = TRUE)
enrichGO_BP = enrichGO(gene = gene, 'org.Hs.eg.db', universe = entrezID, 
ont = 'BP', pvalueCutoff = 0.01, readable = TRUE)

write.table(name, 'ppi/enrichment2/module_GO_MF.txt', sep = '\t', row.names = FALSE,
col.names = FALSE, quote = F, append = TRUE)
write.table(enrichGO_MF, 'ppi/enrichment2/module_GO_MF.txt', sep = '\t', row.names = FALSE,
col.names = FALSE, quote = F, append = TRUE)
write.table(name, 'ppi/enrichment2/module_GO_CC.txt', sep = '\t', row.names = FALSE,
col.names = FALSE, quote = F, append = TRUE)
write.table(enrichGO_CC, 'ppi/enrichment2/module_GO_CC.txt', sep = '\t', row.names = FALSE,
col.names = FALSE, quote = F, append = TRUE)
write.table(name, 'ppi/enrichment2/module_GO_BP.txt', sep = '\t', row.names = FALSE,
col.names = FALSE, quote = F, append = TRUE)
write.table(enrichGO_BP, 'ppi/enrichment2/module_GO_BP.txt', sep = '\t', row.names = FALSE,
col.names = FALSE, quote = F, append = TRUE)


#KEGG富集
enrichKEGG = enrichKEGG(gene = gene, pvalueCutoff = 0.05)

write.table(name, 'ppi/enrichment2/module_kegg.txt', sep = '\t', row.names = FALSE,
col.names = FALSE, quote = F, append = TRUE)
write.table(enrichKEGG, 'ppi/enrichment2/module_kegg.txt', sep = '\t', row.names = FALSE,
col.names = FALSE, quote = F, append = TRUE)

##结果可视化
#GO/KEGG富集结果的可视化
barplot(enrichKEGG, font.size = 12, showCategory = 4)
#按照基因集所富集的功能分类
compare = compareCluster(geneCluster = dif, fun = 'enrichKEGG')



