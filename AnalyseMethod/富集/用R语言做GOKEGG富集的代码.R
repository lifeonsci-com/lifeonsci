#参考文件夹下的  R做通路富集分析.pdf
#示例是给出了三个symbol做BP富集分析，KEGG富集分析，并给出了富集到通路的图
#代码实现的平台是R-3.2.2 镜像选择beijing3(https) 其他全选就可
#要到的包("org.Hs.eg.db")("GSEABase")("GOstats")("KEGG.db")("pathview")

#三个symbol 一定要大写toupper() 转换成entrzID 这是第一种转换方法 我看网页（R做GO KEGG富集分析enrichment analysis - 生物信息 - 生物秀.htm）-在文件夹下
#是用另一种方法转换的，详细请看   R 包描述及操作/各种geneID之间转换
gene=c("APC", "KRAS", "TP53")
egids <- select(org.Hs.eg.db, gene, "ENTREZID", "SYMBOL")[,2]


###################################  网页转换ID   ####################################
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)
head(entrezIDs)
entrezIDs <- as.character(entrezIDs)
###########################################################################

#准备好要富集的id集合

#构建背景集合 gsc
library("org.Hs.eg.db")
frame = toTable(org.Hs.egGO)
goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)
head(goframeData)

goFrame=GOFrame(goframeData,organism="Homo sapiens")
goAllFrame=GOAllFrame(goFrame)
library("GSEABase")
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())


#universe所有的ID
#!GSEAGOHyperGParams构建要的一下参数对象，geneSetCollection背景集，geneIds要富集的id集，universeGeneIds所有ids,ontology要富集的分支BP、MF、CC，pvalueCutoff要卡的阈值可以选着0.1然后自己卡

library("GOstats")
universe = Lkeys(org.Hs.egGO)
params <- GSEAGOHyperGParams(name="My Custom GSEA based annot Params",geneSetCollection=gsc,geneIds = egids,universeGeneIds = universe,ontology = "BP",pvalueCutoff = 0.05,conditional = FALSE,testDirection = "over")
Over <- hyperGTest(params)
head(summary(Over))
#!

#!对富集的结果添加一些项
library(Category)
glist <- geneIdsByCategory(Over)
glist <- sapply(glist, function(.ids) {.sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA);.sym[is.na(.sym)] <- .ids[is.na(.sym)];paste(.sym, collapse=";");	})

head(glist)
bp <- summary(Over)
#添加p值FDR校正的值  添加富集到的Symbols
bp$p.adjust <- p.adjust(bp$Pvalue,method="fdr",n=length(bp$Pvalue));
bp$Symbols <- glist[as.character(bp$GOBPID)];
head(bp)
#!
write.table(bp,file=paste("D:/others/gx/crcBPEnrichmentRes/",crcfs[i],"BP_EnrichmentRes.txt"),quote=F,col.names=T,row.names=F,sep="\t");
#以上bp写出即可


#################################################################################
#清除上面一下重复变量
rm("frame","gsc","universe","glist","Over","params")
gc()
#要导入("KEGG.db")包 

library("KEGG.db")

#！同样要构建背景集合gsc
frame = toTable(org.Hs.egPATH)
keggframeData = data.frame(frame$path_id, frame$gene_id)
head(keggframeData)
keggFrame=KEGGFrame(keggframeData,organism="Homo sapiens")
gsc <- GeneSetCollection(keggFrame, setType = KEGGCollection())
universe = Lkeys(org.Hs.egGO)
#！
#!设置GSEAKEGGHyperGParams参数对象 geneIds = egids，选着合适的pvaluecuttoff 0.1然后自己卡

params <- GSEAKEGGHyperGParams(name="My Custom GSEA based annot Params", geneSetCollection=gsc, geneIds = egids, universeGeneIds = universe, pvalueCutoff = 0.05, testDirection = "over")
kOver <- hyperGTest(params)
head(summary(kOver))
kegg <- summary(kOver)
#!
#！同样添加一些项 p值FDR校正的值  添加富集到的Symbols
kegg$p.adjust <- p.adjust(kegg$Pvalue,method="fdr",n=length(kegg$Pvalue));
library(Category)
glist <- geneIdsByCategory(kOver)
glist <- sapply(glist, function(.ids) { .sym <- mget(.ids, envir=org.Hs.egSYMBOL, ifnotfound=NA);	.sym[is.na(.sym)] <- .ids[is.na(.sym)];	paste(.sym, collapse=";"); 	})
kegg$Symbols <- glist[as.character(kegg$KEGGID)]
head(kegg)
write.table(kegg,file=paste("D:/others/gx/crcKeggEnrichmentRes/",crcfs[i],"Kegg_EnrichmentRes.txt"),quote=F,col.names=T,row.names=F,sep="\t");
#！
#输出富集到通路的图片 这里选着了前三个，当然可以选着kegg对象中所有
#out.suffix="pathview"可添加输出文件后缀
####	for(j in 1:nrow(kegg)){pv.out <- pathview(gene.data, pathway.id=as.character(kegg$KEGGID)[j], species="hsa", out.suffix=paste("pathview_",crcfs[i],sep=""), kegg.native=T)};

library("pathview")
gene.data <- rep(1, length(egids))
names(gene.data) <- egids
for(i in 1:3){pv.out <- pathview(gene.data, pathway.id=as.character(kegg$KEGGID)[i], species="hsa", out.suffix="pathview", kegg.native=T)}
