### 韦恩图

## 1. 安装并加载包VennDiagram
install.packages("VennDiagram")
library(VennDiagram)


## 2. 画维恩图的函数--venn.diagram()
# VennDiagram包中画维恩图的函数是venn.diagram(). 下面我们来看一下venn.diagram()函数的使用及参数说明。
 venn.diagram(x, filename, height = 3000, width = 3000, resolution =
    500, imagetype = "tiff", units = "px", compression =
    "lzw", na = "stop", main = NULL, sub = NULL, main.pos
    = c(0.5, 1.05), main.fontface = "plain",
    main.fontfamily = "serif", main.col = "black",
    main.cex = 1, main.just = c(0.5, 1), sub.pos = c(0.5,
    1.05), sub.fontface = "plain", sub.fontfamily =
    "serif", sub.col = "black", sub.cex = 1, sub.just =
    c(0.5, 1), category.names = names(x), force.unique =
    TRUE, print.mode = "raw", sigdigs = 3, direct.area =
    FALSE, area.vector = 0, hyper.test = FALSE, total.population = NULL, 
    lower.tail = TRUE, ...)

# 部分参数说明：
# x: a list of vectors, e.g: list(A=1:10, B=3:8, C=5:13)
# filename: 设置图形输出文件名
# resolution: 输出图形的清晰度，DPI数值
# imagetype: 输出图形的格式，tiff, png, svg 等
# alpha: 设置每个区块的透明度
# main: 图形标题
# main.fontface: 字体样式，比如斜体，粗体等
# main.fontfamily: 字体，比如Time New Roman等
# 关于调解文字的，不仅可以针对标题调节，还有参数分别针对子标题，维恩图中每个部分（类别）的名字进行字体，大小，和字体样式的设置




# 下面举个我自己刚刚画的一个例子：
# 我有三个向量，分别是wdspWD40, smartWD40, 和pfamWD40. 它们之间可能会有交集，我想用维恩图来可视化，代码如下：
venn.diagram(list(WDSP=wdspWD40,Pfam=pfamWD40,SMART=smartWD40),
    resolution = 300, imagetype = "tiff", alpha=c(0.5,0.5,0.5),
    fill=c("red","yellow","blue"), cat.fontface=4,fontfamily=3,
    main="Intersection of WD40 genes identified by different methods",
    main.cex = 2, main.fontface = 2, main.fontfamily = 3,
    filename = "VennDiagram.tif")
	
	
	
	
	
	
	
	
	
	
	


