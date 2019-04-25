### https://www.sohu.com/a/245475759_777125
library(biomaRt)

#Biomart目前提供了四种数据库，可以使用listMarts()函数查看：
listMarts()
               # biomart               version
# 1 ENSEMBL_MART_ENSEMBL      Ensembl Genes 95
# 2   ENSEMBL_MART_MOUSE      Mouse strains 95
# 3     ENSEMBL_MART_SNP  Ensembl Variation 95
# 4 ENSEMBL_MART_FUNCGEN Ensembl Regulation 95

#选定数据库
my_mart=useMart("ENSEMBL_MART_ENSEMBL")

# 在ensembl数据库中包含了77个数据集，大家可以使用以下代码查看：
datasets=listDatasets(my_mart)

# 为了进一步操作，我们必须从这77个数据集中选择一个数据集，此处当然是选择人类基因的ensembl数据集咯~代码如下：
my_datasets=useDataset("hsapiens_gene_ensembl",mart=my_mart)

# 至此，数据集的工作已经完成了，下面我们就可以在此数据集中，根据需要进行ID的转换。
# 此工作由getBM()函数完成。该函数有四个参数，通常情况下都需要指定！
# 1）attributes参数：该参数可以接受一个字符串向量，用来指定输出的数据类型，就是你要什么，比如entrezgene，hgnc_id。
# 2）filters参数：该参数也可以接受一个字符串向量，用来指定数据的输入类型，比如你的原始信息是基因的ensembl ID，并且有这些基因的染色体位置信息，那么此处的filter就是ensembl ID和chromosome_name等。
# 3）values参数：与filters参数一一对应，即设定filters的具体的值，也就是你的原始数据。
# 4）mart参数：此前定义的数据库，此处就是my_dataset。

# attributes和filters参数具体的取值可以通过listAttributes()和listFilters()函数查看。
# 只需要在函数中定义mart参数为我们当前选择的my_dataset即可。
ensemblID <- protein
gene_names <- getBM(
    filters= "ensembl_gene_id", 
    attributes= c("ensembl_gene_id","hgnc_symbol","entrezgene","description"),
    values= ensemblID,
    mart= my_datasets)










### 实例：

# 用String构建PPIs
# 在String下载人类蛋白质互作数据，String.txt，构建PPIs网络
# String_Symbol转换

library(biomaRt)
### test
# mart = useMart(host = 'grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='hsapiens_gene_ensembl')
# mart=useDataset("hsapiens_gene_ensembl","" mart = mart)

# ensembl_genes <- "ENSP00000000233"

# gene_names <- getBM(
    # filters= "ensembl_peptide_id", 
    # attributes= c("ensembl_peptide_id","hgnc_symbol","description"),
    # values= ensembl_genes,
    # mart= mart)
### String_ID转换
setwd("G:\\Projects\\Shanghai_jinshi\\1. Data\\5. String")
String_mmu <- read.table("String_mmu.txt",sep = " ",stringsAsFactors = F,header = T,fill = TRUE,quote = "")
protein1=data.frame(unlist(strsplit(String_mmu$'protein1', "[.]")))
colnames(protein1) <- c("protein")
protein1 <- unique(protein1)
protein2=data.frame(unlist(strsplit(String_mmu$'protein2', "[.]")))
colnames(protein2) <- c("protein")
protein2 <- unique(protein2)
protein <- rbind(protein1,protein2)
protein <- unique(protein)
protein <- data.frame(protein[-1,])
colnames(protein) <- "protein"

mart = useMart(host = 'grch37.ensembl.org', biomart='ENSEMBL_MART_ENSEMBL', dataset='mmusculus_gene_ensembl')

ensembl_genes <- protein

gene_names <- getBM(
    filters= "ensembl_peptide_id", 
    attributes= c("ensembl_peptide_id","mgi_symbol","entrezgene","description"),
    values= ensembl_genes,
    mart= mart)

write.table(gene_names,"String_Symbol_change.txt",row.names = F,quote = F,sep="\t")

String_mmu1 <- String_mmu
String_mmu1[,3] <- paste("10090.",String_mmu1[,3],sep="")
String_mmu2 <- data.frame(unlist(strsplit(String_mmu1$'protein1', "[.]")),unlist(strsplit(String_mmu1$'protein2', "[.]")),unlist(strsplit(String_mmu1$'combined_score', "[.]")))
colnames(String_mmu2) <- c("protein1","protein2","combined_score")
String_mmu_PPIs_head <- String_mmu2[String_mmu2$protein1 != "10090",]
write.table(String_mmu_PPIs_head,"String_mmu_PPIs_head.txt",row.names = F,quote = F,sep="\t")

String_Symbol_change <- gene_names
String_Symbol_change1 <- data.frame(unique(String_Symbol_change[,1:2]))
String_mmu_PPIs1 <- merge(String_mmu_PPIs_head,String_Symbol_change1,by.x = "protein1",by.y = "ensembl_peptide_id",all=FALSE)
String_mmu_PPIs1 <- String_mmu_PPIs1[,2:4]
colnames(String_mmu_PPIs1)[3] <- c("protein1")
String_mmu_PPIs2 <- merge(String_mmu_PPIs1,String_Symbol_change1,by.x = "protein2",by.y = "ensembl_peptide_id",all=FALSE)
String_mmu_PPIs2 <- String_mmu_PPIs2[,c(3,4,2)]
colnames(String_mmu_PPIs2)[2] <- c("protein2")

String_mmu_PPIs3 <- String_mmu_PPIs2[String_mmu_PPIs2$protein1 != "" & String_mmu_PPIs2$protein2 != "" ,] #去除空值
String_mmu_PPIs3 <- String_mmu_PPIs3[,c(2,1,3)]
colnames(String_mmu_PPIs3)[1:2] <- c("protein1","protein2")
write.table(String_mmu_PPIs3,"String_mmu_PPIs_Symbol.txt",row.names = F,quote = F,sep="\t")
String_mmu_PPIs_Symbol = String_mmu_PPIs3