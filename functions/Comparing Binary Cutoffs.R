
load('data/esets_bin_STAT117.RData')
probeset = "TGFBR2"
probeset = "CXCL12"
esets_probeset = esets_bin[sapply(esets_bin, function(x) probeset %in% featureNames(x))]
XX = sapply(1:length(esets_probeset), function(i) { exprs(esets_probeset[[i]])[probeset,]})
YY = sapply(1:length(esets_probeset), function(i) { pData(esets_probeset[[i]])[,"debulking"] == "optimal" })
XX_new <- XX[17]
YY_new <- YY[17]

XX_new <- XX_new[[1]]
YY_new <- YY_new[[1]]


XX_orig <- XX_orig[[1]]
YY_orig <- YY_orig[[1]]

XXd = sapply(1:length(esets_probeset), function(i) { 1 * ( XX[[i]] > 0 ) })
TAB = as.data.frame(matrix(nrow=length(esets_probeset),ncol=4))
colnames(TAB) = c("x1","n1","x0","n0")
for (i in 1:length(esets_probeset)) {
  if (sum (!is.na(YY[[i]])) > 0) {
  TAB[i,"x1"] = sum( XXd[[i]] == 1 & YY[[i]] == 1 & !is.na(YY[[i]]))
  TAB[i,"x0"] = sum( XXd[[i]] == 1 & YY[[i]] == 0 & !is.na(YY[[i]]))
  TAB[i,"n1"] = sum( YY[[i]] == 1 & !is.na(YY[[i]]))
  TAB[i,"n0"] = sum( YY[[i]] == 0 & !is.na(YY[[i]]))
  }
}
TAB = TAB[rowSums(is.na(TAB)) != ncol(TAB), ]
TAB



#############################################################

load('data/esets.RData')
probeset = "TGFBR2"
esets_probeset = esets[sapply(esets, function(x) probeset %in% featureNames(x))]
XX = sapply(1:length(esets_probeset), function(i) { exprs(esets_probeset[[i]])[probeset,] })
YY = sapply(1:length(esets_probeset), function(i) { pData(esets_probeset[[i]])[,"debulking"] == "optimal" })
XX_orig <- XX[15]
YY_orig <- YY[15]