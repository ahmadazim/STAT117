
load('/Users/albert/Dropbox/GP_Teaching/STA117/Textbook117/data/esets_meta_STAT117.RData')
#probesets = c("POSTN", "CXCL12", "CXCL14", "BRCA1")
library(curatedOvarianData)

# POSTN
probeset = "POSTN"
esets_probeset = esets_meta[sapply(esets_meta, function(x) probeset %in% featureNames(x))]
ind.rm <- which(names(esets_probeset) %in% c("GSE17260_eset","GSE32062.GPL6480_eset","GSE32063_eset","GSE49997_eset","GSE51088_eset","GSE8842_eset"))
if(length(ind.rm)>0){
	esets_probeset <- esets_probeset[-ind.rm]
}
ind.rm.debulk <- which(names(esets_probeset) %in% c("GSE13876_eset","GSE14764_eset","GSE19829.GPL570_eset","GSE19829.GPL8300_eset","GSE26193_eset"))
if(length(ind.rm.debulk)>0){
	esets_probeset <- esets_probeset[-ind.rm.debulk]
}
#esets_probeset <- esets_probeset[-which(names(esets_probeset) == "GSE30161_eset")]
if(probeset == "BRCA1"){
	esets_probeset <- esets_probeset[-which(names(esets_probeset) == "GSE30009_eset")]
}

XX = sapply(1:length(esets_probeset), function(i) { exprs(esets_probeset[[i]])[probeset,]})
YY = sapply(1:length(esets_probeset), function(i) { pData(esets_probeset[[i]])[,"debulking"] == "optimal" })

for (i in 1:length(esets_probeset)) {
  if (sum (is.na(YY[[i]])) > 0) {
    XX[[i]] <- XX[[i]][!is.na(YY[[i]])]
    YY[[i]] <- YY[[i]][!is.na(YY[[i]])]
  }
}

par(mfrow = c(2,3))
for(i in 1:length(XX)){
	if(i == 8){
		ffx = density(XX[[i]][YY[[i]]==FALSE],from=2,to=14)$x
		ff0 = density(XX[[i]][YY[[i]]==FALSE],from=2,to=14)$y
		ff1 = density(XX[[i]][YY[[i]]==TRUE],from=2,to=14)$y
	} else {
		ffx = density(XX[[i]][YY[[i]]==FALSE],from=0,to=15)$x
		ff0 = density(XX[[i]][YY[[i]]==FALSE],from=0,to=15)$y
		ff1 = density(XX[[i]][YY[[i]]==TRUE],from=0,to=15)$y
	}
	plot(ffx,ff0,xlab=c("Log Expression",probeset),ylab="Density", main = paste(names(esets_probeset)[i], " N=", length(YY[[i]])),ylim=c(0,1), col = "blue",type="l",lwd=2)
	lines(ffx,ff1, col = "green",lwd=2)
	rug(XX[[i]][YY[[i]]==FALSE],ticksize = .05,lwd=2, col = "blue")
	rug(XX[[i]][YY[[i]]==TRUE],ticksize = .05,lwd=2, col = "green")
}



# CXCL12
probeset = "CXCL12"
esets_probeset = esets_meta[sapply(esets_meta, function(x) probeset %in% featureNames(x))]
ind.rm <- which(names(esets_probeset) %in% c("GSE17260_eset","GSE32062.GPL6480_eset","GSE32063_eset","GSE49997_eset","GSE51088_eset","GSE8842_eset"))
if(length(ind.rm)>0){
	esets_probeset <- esets_probeset[-ind.rm]
}
ind.rm.debulk <- which(names(esets_probeset) %in% c("GSE13876_eset","GSE14764_eset","GSE19829.GPL570_eset","GSE19829.GPL8300_eset","GSE26193_eset"))
if(length(ind.rm.debulk)>0){
	esets_probeset <- esets_probeset[-ind.rm.debulk]
}
#esets_probeset <- esets_probeset[-which(names(esets_probeset) == "GSE30161_eset")]
XX = sapply(1:length(esets_probeset), function(i) { exprs(esets_probeset[[i]])[probeset,]})
YY = sapply(1:length(esets_probeset), function(i) { pData(esets_probeset[[i]])[,"debulking"] == "optimal" })
for (i in 1:length(esets_probeset)) {
  if (sum (is.na(YY[[i]])) > 0) {
    XX[[i]] <- XX[[i]][!is.na(YY[[i]])]
    YY[[i]] <- YY[[i]][!is.na(YY[[i]])]
  }
}
par(mfrow = c(2,3))
for(i in 1:length(XX)){
	if(i == 8){
		ffx = density(XX[[i]][YY[[i]]==FALSE],from=0,to=15)$x
		ff0 = density(XX[[i]][YY[[i]]==FALSE],from=0,to=15)$y
		ff1 = density(XX[[i]][YY[[i]]==TRUE],from=0,to=15)$y
	} else {
		ffx = density(XX[[i]][YY[[i]]==FALSE],from=2,to=15)$x
		ff0 = density(XX[[i]][YY[[i]]==FALSE],from=2,to=15)$y
		ff1 = density(XX[[i]][YY[[i]]==TRUE],from=2,to=15)$y
	}
	plot(ffx,ff0,xlab=c("Log Expression",probeset),ylab="Density", main = paste(names(esets_probeset)[i], " N=", length(YY[[i]])),ylim=c(0,1), col = "blue",type="l",lwd=2)
	lines(ffx,ff1, col = "green",lwd=2)
	rug(XX[[i]][YY[[i]]==FALSE],ticksize = .05,lwd=2, col = "blue")
	rug(XX[[i]][YY[[i]]==TRUE],ticksize = .05,lwd=2, col = "green")
}



# CXCL14
probeset = "CXCL14"
esets_probeset = esets_meta[sapply(esets_meta, function(x) probeset %in% featureNames(x))]
ind.rm <- which(names(esets_probeset) %in% c("GSE17260_eset","GSE32062.GPL6480_eset","GSE32063_eset","GSE49997_eset","GSE51088_eset","GSE8842_eset"))
if(length(ind.rm)>0){
	esets_probeset <- esets_probeset[-ind.rm]
}
ind.rm.debulk <- which(names(esets_probeset) %in% c("GSE13876_eset","GSE14764_eset","GSE19829.GPL570_eset","GSE19829.GPL8300_eset","GSE26193_eset"))
if(length(ind.rm.debulk)>0){
	esets_probeset <- esets_probeset[-ind.rm.debulk]
}
#esets_probeset <- esets_probeset[-which(names(esets_probeset) == "GSE30161_eset")]

XX = sapply(1:length(esets_probeset), function(i) { exprs(esets_probeset[[i]])[probeset,]})
YY = sapply(1:length(esets_probeset), function(i) { pData(esets_probeset[[i]])[,"debulking"] == "optimal" })
for (i in 1:length(esets_probeset)) {
  if (sum (is.na(YY[[i]])) > 0) {
    XX[[i]] <- XX[[i]][!is.na(YY[[i]])]
    YY[[i]] <- YY[[i]][!is.na(YY[[i]])]
  }
}
par(mfrow = c(2,3))
for(i in 1:length(XX)){
	if(i == 8){
		ffx = density(XX[[i]][YY[[i]]==FALSE],from=2,to=15)$x
		ff0 = density(XX[[i]][YY[[i]]==FALSE],from=2,to=15)$y
		ff1 = density(XX[[i]][YY[[i]]==TRUE],from=2,to=15)$y
	} else {
		ffx = density(XX[[i]][YY[[i]]==FALSE],from=0,to=15)$x
		ff0 = density(XX[[i]][YY[[i]]==FALSE],from=0,to=15)$y
		ff1 = density(XX[[i]][YY[[i]]==TRUE],from=0,to=15)$y
	}
	plot(ffx,ff0,xlab=c("Log Expression",probeset),ylab="Density", main = paste(names(esets_probeset)[i], " N=", length(YY[[i]])),ylim=c(0,1), col = "blue",type="l",lwd=2)
	lines(ffx,ff1, col = "green",lwd=2)
	rug(XX[[i]][YY[[i]]==FALSE],ticksize = .05,lwd=2, col = "blue")
	rug(XX[[i]][YY[[i]]==TRUE],ticksize = .05,lwd=2, col = "green")
}


# BRCA1
probeset = "BRCA1"
esets_probeset = esets_meta[sapply(esets_meta, function(x) probeset %in% featureNames(x))]
ind.rm <- which(names(esets_probeset) %in% c("GSE17260_eset","GSE32062.GPL6480_eset","GSE32063_eset","GSE49997_eset","GSE51088_eset","GSE8842_eset"))
if(length(ind.rm)>0){
	esets_probeset <- esets_probeset[-ind.rm]
}
ind.rm.debulk <- which(names(esets_probeset) %in% c("GSE13876_eset","GSE14764_eset","GSE19829.GPL570_eset","GSE19829.GPL8300_eset","GSE26193_eset"))
if(length(ind.rm.debulk)>0){
	esets_probeset <- esets_probeset[-ind.rm.debulk]
}
#esets_probeset <- esets_probeset[-which(names(esets_probeset) == "GSE30161_eset")]

XX = sapply(1:length(esets_probeset), function(i) { exprs(esets_probeset[[i]])[probeset,]})
YY = sapply(1:length(esets_probeset), function(i) { pData(esets_probeset[[i]])[,"debulking"] == "optimal" })
for (i in 1:length(esets_probeset)) {
  if (sum (is.na(YY[[i]])) > 0) {
    XX[[i]] <- XX[[i]][!is.na(YY[[i]])]
    YY[[i]] <- YY[[i]][!is.na(YY[[i]])]
  }
}
par(mfrow = c(2,3))
for(i in 1:length(XX)){
	if(i %in% c(6)){
		ffx = density(XX[[i]][YY[[i]]==FALSE],from=2,to=8)$x
		ff0 = density(XX[[i]][YY[[i]]==FALSE],from=2,to=8)$y
		ff1 = density(XX[[i]][YY[[i]]==TRUE],from=2,to=8)$y
	} else {
		ffx = density(XX[[i]][YY[[i]]==FALSE],from=4,to=9)$x
		ff0 = density(XX[[i]][YY[[i]]==FALSE],from=4,to=9)$y
		ff1 = density(XX[[i]][YY[[i]]==TRUE],from=4,to=9)$y
	}
	if(i %in% c(1,2,3,5)){
		plot(ffx,ff0,xlab=c("Log Expression",probeset),ylab="Density", main = paste(names(esets_probeset)[i], " N=", length(YY[[i]])),ylim=c(0,4), col = "blue",type="l",lwd=2)
	} else{
		plot(ffx,ff0,xlab=c("Log Expression",probeset),ylab="Density", main = paste(names(esets_probeset)[i], " N=", length(YY[[i]])),ylim=c(0,1), col = "blue",type="l",lwd=2)
	}
	lines(ffx,ff1, col = "green",lwd=2)
	rug(XX[[i]][YY[[i]]==FALSE],ticksize = .05,lwd=2, col = "blue")
	rug(XX[[i]][YY[[i]]==TRUE],ticksize = .05,lwd=2, col = "green")
}





