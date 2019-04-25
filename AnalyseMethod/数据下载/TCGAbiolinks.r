### TCGAbiolinks
## 肺癌mRNA数据的下载。可推广到miR、CNV、甲基化等数据。

setwd("C:\\Users\\waiting\\Desktop")

library(TCGAbiolinks)

projectid <- "TCGA-PAAD"

query.count <- GDCquery(project= projectid,

                       data.category = "Transcriptome Profiling",

                       data.type = "Gene Expression Quantification",

                       workflow.type = "HTSeq - Counts")   # 需注意“-”前后的空格

					   
					   
# 下载数据

GDCdownload(query.count)

# 获得表达矩阵

dataAssay = GDCprepare(query.count, summarizedExperiment = F)

rownames(dataAssay) = as.character(dataAssay[,1])

# dataAssay就是矩阵了，它此时在R的环境变量里、也就是在计算机内存中。你可以在R中对它进行进一步的分析。

# 也可以用write.table或write.csv命令把它从R里保存出来到硬盘，并保存为csv的格式，就可以用excel打开了。

write.csv(dataAssay, "TCGA-matrix.csv")  # 此时，保存的文件名为“TCGA-matrix.csv”





















