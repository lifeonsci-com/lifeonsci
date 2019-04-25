
library(psych)
### fun_to_corr: 
#### input: heat_in_1, heat_in_2: f_ra/g_ra/biochem/kegg_in...you can select the feature first
#### output: list(), t_cor is correlation index(-1~1), t_p is raw p-value, t_p_sig is formatted with significant levels
#### formatted significant levels: 0~0.001: ***; 0.001~0.01: **; 0.05~0.1: *; >0.1 +
fun_to_corr <- function(heat_in_1, heat_in_2) {
  t <- corr.test(heat_in_1, heat_in_2, use = "pairwise", method = "spearman", adjust = "fdr")
  t_cor <- data.frame(t$r, check.names = FALSE)
  t_p <- data.frame(t$p, check.names = FALSE)
  cut_sig <- function(p) {
    out <- cut(p, breaks = c(0, 0.001,0.01,0.05,1), include.lowest = T, labels = c("***", "**", "*", "+"))
    return(out)
  }
  t_p_sig <- apply(t_p, 2, cut_sig)
  rownames(t_p_sig) <- rownames(t_p)
  return(list(t_cor = t_cor, t_p_sig = t_p_sig, t_p = t_p))
}
corr_test=fun_to_corr(heat_in_1,heat_in_2)






### 如何将相关系数以及显著性水平p-value整合进一个矩阵内，可以自定义一个函数flattenCorrMatrix。
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values


flattenCorrMatrix <- function(cormat, pmat) {

	ut <- upper.tri(cormat) 
	data.frame( row = rownames(cormat)[row(cormat)[ut]],column = rownames(cormat)[col(cormat)[ut]], cor =(cormat)[ut], p = pmat[ut])

}










flattenCorrMatrix <- function(cormat, pmat, sigmat) {

	ut <- upper.tri(cormat) 
	data.frame( row = rownames(cormat)[row(cormat)[ut]],column = rownames(cormat)[col(cormat)[ut]], cor =(cormat)[ut], p = pmat[ut],sig=sigmat[ut] )

}
corr_test=fun_to_corr(heat_in_1,heat_in_2)
flattenCorrMatrix(corr_test$t_cor, corr_test$t_p, corr_test$t_p_sig)

