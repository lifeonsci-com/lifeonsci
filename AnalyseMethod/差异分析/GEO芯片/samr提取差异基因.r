#R-3.2.4
rm(list=ls())
gc()

setwd("D:\\projects_a\\LIHC\\DATA\\GSE115018\\lncRNA")
exp=read.table("Exp.txt",header=T,sep="\t",row.names = 1)
geneName=rownames(exp)
length(unique(geneName))
exp=apply(exp,2,as.numeric)
dim(exp)
# install.packages("samr")
suppressPackageStartupMessages(library(samr))

help(package="samr")
colnames(exp)

group=c(rep(1,times=10),rep(2,times=10))
data=list(x=exp,y=group, 
          geneid=as.character(1:nrow(exp)),
          genenames=geneName, 
          logged2=TRUE
)

#resp.type=c("Quantitative","Two class unpaired","Survival","Multiclass","One class", "Two class paired","Two class unpaired timecourse","One class timecourse","Two class paired timecourse", "Pattern discovery")
#可以用的数据类型


samr.obj<-samr(data, resp.type="Two class unpaired", nperms=1000)
#the "min.foldchange=0.1" means that the fold change for the two groups should be >0.1
delta.table <- samr.compute.delta.table(samr.obj, min.foldchange=0.1,nvals=200)
siggenes.table <- samr.compute.siggenes.table(samr.obj, del=0, data, delta.table,all.genes=TRUE)
str(siggenes.table)

a <- siggenes.table$genes.up; # all up regulated genes
b <- siggenes.table$genes.lo; # all down regulated genes
c <- rbind(a,b)
head(c)
head(b,100)
#Let's write all FDR < 5% DEGs into a file
#the 8th column of siggenes.table is FDR
#the 7th column of siggenes.table is FC=B/A
#the 9th column of siggenes.table is FC=B/A
lo <- c[as.numeric(c[,8])<5,]
lo=as.data.frame(lo)
lo$logFC=log(as.numeric(as.matrix(lo[,7])),2)
head(lo)
write.csv(lo,"DEGs_samr.csv")
