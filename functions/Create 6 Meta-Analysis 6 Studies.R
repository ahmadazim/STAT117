# which(pData(esets_bin$TCGA.RNASeqV2_eset)[,"unique_patient_ID"] %in% pData(esets_bin$TCGA_eset)[,"unique_patient_ID"]) 
# TCGA.RNASeqV2_eset subset of TCGA_eset
load('Textbook117/data/esets_meta_STAT117.RData')
library(curatedOvarianData)

ind.rm <- which(names(esets_meta) %in% c("GSE17260_eset","GSE32062.GPL6480_eset","GSE32063_eset","GSE49997_eset","GSE51088_eset","GSE8842_eset"))
if(length(ind.rm)>0){
	esets_meta <- esets_meta[-ind.rm]
}
ind.rm.debulk <- which(names(esets_meta) %in% c("GSE13876_eset","GSE14764_eset","GSE19829.GPL570_eset","GSE19829.GPL8300_eset","GSE26193_eset"))
if(length(ind.rm.debulk)>0){
	esets_meta <- esets_meta[-ind.rm.debulk]
}

esets <- list()

for(i in 1:length(names(esets_meta))){
  XX = as.matrix(cbind(exprs(esets_meta[[i]])))
  YY = 1 * as.vector(pData(esets_meta[[i]])[,"debulking"]=="optimal")
  XX = XX[,!is.na(YY)]
  YY = YY[!is.na(YY)]
  esets_temp <- list(XX=XX,YY=YY)
  esets <- append(esets, list(esets_temp))
  names(esets)[i] <- names(esets_meta)[i]
}

# To call out XX and YY for first study
esets[[1]]$XX
esets[[1]]$YY

# To call out XX and YY for TCGA_eset
esets$TCGA_eset$XX
esets$TCGA_eset$YY

save(esets, file = "List_of_6_Studies.RData")








