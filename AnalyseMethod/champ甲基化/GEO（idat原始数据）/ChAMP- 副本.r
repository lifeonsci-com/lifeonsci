### GEO数据
### arraytype="450K" 数据为850K改 arraytype="850K"或者arraytype = 'EPIC'
## ChAMP
# https://blog.csdn.net/herokoking/article/details/82345263

# GEO数据下载地址https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE89648
# 甲基化芯片的分组也跟表达芯片一样，两组各4个样本，在GEO网站上找到其对应的GSM号，
# 然后参照ChAMP包给的测试例子中的lung_test_set.csv作为样式，手动做个sample.csv，以便能让ChAMP包识别





setwd("C:\\Users\\waiting\\Desktop\\ChAMP")
library(ChAMP)

##数据过滤
# 准备好上述文件后，接下来则是读入数据，
# 使用ChAMP包的champ.load函数，
# 按照其说法，champ.load函数其实做了两个步骤，
# 第一是读入数据（champ.import函数），
# 第二是过滤数据（champ.filter函数），
# 所以我们需要了解champ.load以什么标准过滤：
# 过滤掉detection p-value大于0.01的探针
# 过滤掉在至少5%样本中bead count小于3的探针
# 过滤掉非GpC位点的探针
# 过滤掉所有SNP相关的探针
# 过滤掉multi-hit探针，即映射到多个位置的
# 过滤掉X和Y染色体上的探针
testDir = "C:/Users/waiting/Desktop/ChAMP/450K"
myLoad <- champ.load(testDir,arraytype="450K")



# ChAMP包还提供了CpG.GUI()函数用于展示distribution of CpGs on chromosome, CpG island, TSS reagions. e.g，
# 其实是几张简单展示myLoad$beta分布的图片，但是界面很友好，也难怪依赖包里要shiny包了
CpG.GUI(CpG = rownames(myLoad$beta))



## 数据质控QC
# champ.QC()和QC.GUI函数用于检查看看过滤后的数据是否符合要求，能否进行下一步的分析。
# champ.QC()函数结果三张图，
# 如有 mdsPlot (Multidimensional Scaling Plot)，主要看看不同分组样本是否分开；
# densityPlot ，每个样本的beta值的分布图，主要看看有无异常的样本；
# dendrogram ，样本的聚类图
champ.QC(beta = myLoad$beta, pheno = myLoad$pd$Sample_Name)

# MDS 1000 most variable positions.png
QC.GUI(beta=myLoad$beta)



## 数据标准化
# 简单的说就是消除测序芯片上的两种测序技术之间的差异。
# ChAMP包标准化方法有BMIQ、SWAN、PBC以及FunctionalNormliazation
# 默认是采用BMIQ
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5)



## SVD plot
# SVD(singular value decomposition) 只为检测变异组分的显著性
# 主要用于查看变异的主要成分是生物学处理等影响的，还是技术因素所造成的
# 如果是后者则需要进行后续的批次校正。从下面代码也可看出其主要根据myNorm（校正后的探针的beta）和pd（样本信息）来计算变异程度
champ.SVD(beta=myNorm,pd=myLoad$pd)



## 批次效应校正
# 如果上面SVD Plot显示红色块并不在Sample_Group行，而是较多的在Array行，
# 那么则需要对这次数据进行批次校正，
# 使用champ.runCombat函数，batchname则是根据上述SVD plot出现的红块位置而定，决定哪个batch factor需要被校正
# 这个过程比较耗时，做完后最好再用champ.SVD检查下
myCombat <- champ.runCombat(beta=myNorm,pd=myLoad$pd,batchname=c("Slide"))
champ.SVD(beta=myNorm,pd=myLoad$pd)



# 甲基化探针差异分析
# 比较常见的差异甲基化分析是DMP（Differential Methylation Probe）
# 用于找出差异的甲基化位点，然后ChAMP包已经将分析过程（主要还是基于Limma包，这个包是专门用于分析芯片差异表达的包）
# 然后根据你的分组信息，获得差异甲基化位点；最后用DMP.GUI查看下结果
myDMP <- champ.DMP(beta = myNorm,pheno=myLoad$pd$Sample_Group)
DMP.GUI(DMP=myDMP[[1]],beta=myNorm,pheno=myLoad$pd$Sample_Group)



# 除了可以分析差异甲基化位点外
# 还可以分析DMR（Differential Methylation Regions）
# 用于找出差异甲基化片段
# 用的函数也很好记champ.DMR()
# 该函数现在只支持两组样本进行差异分析
myDMR <- champ.DMR(beta=myNorm,pheno=myLoad$pd$Sample_Group,method="Bumphunter")

# 这里需要注意的是method的选择
# 主要有Bumphunter, ProbeLasso and DMRcate三者可以选择
# Bumphunter ，groups all probes into small clusters (or regions), then applies random permutation method to estimate candidate DMRs，不用依赖上述差异位点分析的结果
# ProbeLasso ，the final data frame is distilled from a limma output of probes and their association statistics，必须要有champ.DMP()的结果作为输入
# DMRcate ，a new method recently incoporated into ChAMP,pass the square of the moderated t statistic calculated on each 450K probe to DMR-finding function（这个还是看文档吧）



# 除了DMR外还有一个DMB（Differential Methylation Blocks），也是用来识别大片段的差异甲基化
myBlock <- champ.Block(beta=myNorm,pheno=myLoad$pd$Sample_Group,arraytype="450K")











