# which(pData(esets_bin$TCGA.RNASeqV2_eset)[,"unique_patient_ID"] %in% pData(esets_bin$TCGA_eset)[,"unique_patient_ID"]) 
# TCGA.RNASeqV2_eset subset of TCGA_eset
setwd('/Users/albert/Dropbox/GP_Teaching/STA117/Textbook117/figures/Gene_Graphs_Abundance_Histograms_for_Meta_Analysis')
load('/Users/albert/Dropbox/GP_Teaching/STA117/Textbook117/data/esets_meta_STAT117.RData')
library(curatedOvarianData)

ind.rm <- which(names(esets_meta) %in% c("GSE17260_eset","GSE32062.GPL6480_eset","GSE32063_eset","GSE49997_eset","GSE51088_eset","GSE8842_eset"))
if(length(ind.rm)>0){
	esets_meta <- esets_meta[-ind.rm]
}
ind.rm.debulk <- which(names(esets_meta) %in% c("GSE13876_eset","GSE14764_eset","GSE19829.GPL570_eset","GSE19829.GPL8300_eset","GSE26193_eset"))
if(length(ind.rm.debulk)>0){
	esets_meta <- esets_meta[-ind.rm.debulk]
}

CompSummaries = function(XX,YY){
NGenes = nrow(XX)
SSS = data.frame(matrix(NA,NGenes,4))
colnames(SSS) = 
  c("abundance","difference","variance","SNratio")
for (gg in 1:NGenes){
  SSS[gg,"abundance"] = mean(XX[gg,])
  SSS[gg,"difference"] = 
    mean( XX[gg,YY==1] ) - mean ( XX[gg,YY==0] )
  SSS[gg,"variance"] = var( XX[gg,] )
  SSS[gg,"SNratio"] = 
    SSS[gg,"difference"] / sqrt( SSS[gg,"variance"] )
}
return(SSS)
}

names(esets_meta)
par(mfrow = c(2,3))

# GSE18520_eset
XX = as.matrix(cbind(exprs(esets_meta$GSE18520_eset)))
YY = 1 * as.vector(pData(esets_meta$GSE18520_eset)[,"debulking"]=="optimal")
XX = XX[,!is.na(YY)];
YY = YY[!is.na(YY)]
GeneSummaries = CompSummaries(XX,YY)
hist(GeneSummaries[,"abundance"],nclass=100, main = "Abundance Histogram for GSE18520_eset")

# GSE26712_eset
XX = as.matrix(cbind(exprs(esets_meta$GSE26712_eset)))
YY = 1 * as.vector(pData(esets_meta$GSE26712_eset)[,"debulking"]=="optimal")
XX = XX[,!is.na(YY)];
YY = YY[!is.na(YY)]
GeneSummaries = CompSummaries(XX,YY)
hist(GeneSummaries[,"abundance"],nclass=100, main = "Abundance Histogram for GSE26712_eset")

# GSE30161_eset
XX = as.matrix(cbind(exprs(esets_meta$GSE30161_eset)))
YY = 1 * as.vector(pData(esets_meta$GSE30161_eset)[,"debulking"]=="optimal")
XX = XX[,!is.na(YY)];
YY = YY[!is.na(YY)]
GeneSummaries = CompSummaries(XX,YY)
hist(GeneSummaries[,"abundance"],nclass=100, main = "Abundance Histogram for GSE30161_eset")

# GSE9891_eset
XX = as.matrix(cbind(exprs(esets_meta$GSE9891_eset)))
YY = 1 * as.vector(pData(esets_meta$GSE9891_eset)[,"debulking"]=="optimal")
XX = XX[,!is.na(YY)];
YY = YY[!is.na(YY)]
GeneSummaries = CompSummaries(XX,YY)
hist(GeneSummaries[,"abundance"],nclass=100, main = "Abundance Histogram for GSE9891_eset")

# PMID17290060_eset
XX = as.matrix(cbind(exprs(esets_meta$PMID17290060_eset)))
YY = 1 * as.vector(pData(esets_meta$PMID17290060_eset)[,"debulking"]=="optimal")
XX = XX[,!is.na(YY)];
YY = YY[!is.na(YY)]
GeneSummaries = CompSummaries(XX,YY)
hist(GeneSummaries[,"abundance"],nclass=100, main = "Abundance Histogram for PMID17290060_eset")

# TCGA_eset
XX = as.matrix(cbind(exprs(esets_meta$TCGA_eset)))
YY = 1 * as.vector(pData(esets_meta$TCGA_eset)[,"debulking"]=="optimal")
XX = XX[,!is.na(YY)];
YY = YY[!is.na(YY)]
GeneSummaries = CompSummaries(XX,YY)
hist(GeneSummaries[,"abundance"],nclass=100, main = "Abundance Histogram for TCGA_eset")

















