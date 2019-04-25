# 相关性分析

pairs(iris[1:4], main = "Anderson's Iris Data -- 3 species",

pch = 21, bg = c("red", "green3", "blue")[unclass(iris$Species)])

# 等同于 pairs(~Sepal.Length+Sepal.Width+Petal.Length+Petal.Width, data=iris,main = "Anderson's Iris Data -- 3 species",pch = 21,
bg = c("red", "green3", "blue")[unclass(iris$Species)]

# # 1. 自定义函数pannel.cor：显示两两变量间的相关系数，相关系数越大字号越大。 
# panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor){ 
# usr = par("usr");on.exit(par(usr)) par(usr = c(0, 1, 0, 1)) r <- abs(cor(x, y)) txt <- format(c(r, 0.123456789), digits = digits)[1] txt <- paste0(prefix, txt) if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt) text(0.5, 0.5, txt, cex = cex.cor * r)
# }

# # 2. 自定义函数pannel.hist：展示各个变量的直方图 
# panel.hist <- function(x, ...) { usr <- par("usr"); on.exit(par(usr)) par(usr = c(usr[1:2], 0, 1.5) ) h <- hist(x, plot = FALSE) breaks <- h$breaks; nB <- length(breaks) y <- h$counts; y <- y/max(y) rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...) } 

# # 3. 自定义函数panel.ls：绘制散点图，并为其添加线性拟合直线 
# panel.lm<-function(x,y,col=par("col"),bg=NA,pch=par("pch"), cex=1,col.smooth="black",...){ points(x,y,pch=pch,col=col,bg=bg,cex=cex) abline(stats::lm(y~x),col=col.smooth,...) } 

# # 4. 用相关系数（pannel.cor）替代默认图形上三角的散点图，用直方图（pannel.hist）替代默认图形对角线的变量名称， 用添加线性拟合线的散点图（panel.ls）代替默认图形下三角的散点图。 
# pairs(iris[1:4], main = "Anderson's Iris Data -- 3 species", pch = 21, bg = c("red", "green3", "blue")[unclass(iris$Species)], diag.panel=panel.hist, upper.panel=panel.cor, lower.panel=panel.lm)

# gpairs {gpairs}
install.packages("gpairs")
library(gpairs)
gpairs(iris, upper.pars = list(scatter = 'stats'), scatter.pars = list(pch = substr(as.character(iris$Species), 1, 1), col = as.numeric(iris$Species)), stat.pars = list(verbose = TRUE))

# corrgram {corrgram}
install.packages("corrgram")
library(corrgram)
vars2 <- c("Assists","Atbat","Errors","Hits","Homer","logSal", "Putouts","RBI","Runs","Walks","Years")
corrgram(baseball[vars2], order=TRUE, main="Baseball data PC2/PC1 order", lower.panel=panel.shade, upper.panel=panel.pie)
corrgram(auto, order=TRUE, main="Auto data (PC order)",lower.panel=corrgram::panel.ellipse,upper.panel=panel.bar, diag.panel=panel.minmax,col.regions=colorRampPalette(c("darkgoldenrod4", "burlywood1", "darkkhaki", "darkgreen")))

# corrplot {corrplot}
install.packages("corrplot")
library(corrplot)

data(mtcars)

M <- cor(mtcars)

corrplot(M, order = "AOE", type = "upper", tl.pos = "d")

corrplot(M, add = TRUE, type = "lower", method = "number", order = "AOE",

diag = FALSE, tl.pos = "n", cl.pos = "n")

###
res1 <- cor.mtest(mtcars, conf.level = 0.95)

corrplot(M, method="ellipse",p.mat = res1$p, sig.level = 0.2,order = "AOE", type = "upper", tl.pos = "d")

corrplot(M, add = TRUE, p.mat = res1$p, sig.level = 0.2,type = "lower", method = "number", order = "AOE",

diag = FALSE, tl.pos = "n", cl.pos = "n")

###
dat = matrix(c(1:10,10:1), nrow = 10)

wb <- c("white", "black")

corrplot(t(dat), method="pie", is.corr = F, cl.pos = "n", tl.pos = "n",

cl.lim = c(1,10),col = wb, bg = "gold2")


# ggpairs {GGally}
install.packages("GGally")
library(GGally)

ggpairs(flea, columns = 2:4, ggplot2::aes(colour=species))

###
library(ggplot2)

diamonds.samp <- diamonds[sample(1:dim(diamonds)[1], 1000), ]

ggpairs(diamonds.samp[, c(1:2,5,7)],mapping = aes(color = cut),lower = list(continuous = wrap("density", alpha = 0.5), combo = "dot_no_facet"),title = "Diamonds")


# # coplot {graphics}
# par(mar = rep(0, 4), mgp = c(2, 0.5, 0))
# install.packages("maps")
# library(maps)

# coplot(lat ~ long | depth, data = quakes, number = 4,ylim = c(-45, -10.72), panel = function(x, y){map("world2", regions = c("New Zealand","Fiji"),add = TRUE, lwd = 0.1, fill = TRUE,col = "lightgray")text(180, -13, "Fiji", adj = 1)text(170, -35, "NZ")points(x, y, col = rgb(0.2, 0.2, 0.2, 0.5))})


# splom {lattice}
library(lattice)

super.sym <- trellis.par.get("superpose.symbol")

splom(~iris[1:4], groups = Species, data = iris,panel = panel.superpose,key = list(title = "Three Varieties of Iris",columns = 3,points = list(pch = super.sym$pch[1:3],col = super.sym$col[1:3]),text = list(c("Setosa", "Versicolor", "Virginica"))))

# scatterplotMatrix {car}
install.packages("car")
library(car)

scatterplotMatrix(~ income + education + prestige | type, data=Duncan)

# ggscatmat {GGally}
library(GGally)

data(flea)

ggscatmat(flea, columns = 2:4, color = "species")


# cpairs {gclus}
library(gclus)

data(USJudgeRatings)

judge.cor <- cor(USJudgeRatings)

judge.color <- dmat.color(judge.cor)

cpairs(USJudgeRatings,panel.colors=judge.color,pch=".",gap=.5)






















