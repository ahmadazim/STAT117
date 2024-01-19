makeBIO = function(probeset,small=FALSE) {
  esets_probeset_m = esets_meta[sapply(esets_meta, function(x) probeset %in% featureNames(x))]
  XX = sapply(1:length(esets_probeset_m), function(i) { exprs(esets_probeset_m[[i]])[probeset,] })
  YY = sapply(1:length(esets_probeset_m), function(i) { 1 * ( pData(esets_probeset_m[[i]])[,"debulking"] == "optimal" ) })
  allNA = sapply (1:length(esets_probeset_m), function(i) { sum(is.na(YY[[i]])) == length(YY[[i]]) } )
  XX = XX[ (1:length(esets_probeset_m))[!allNA] ]
  YY = YY[ (1:length(esets_probeset_m))[!allNA] ]
  SS = length(YY)
  XX = sapply (1:SS, function(i) { XX[[i]][ !is.na(YY[[i]]) ] } )
  YY = sapply (1:SS, function(i) { YY[[i]][ !is.na(YY[[i]]) ] } )
  NS = sapply (1:SS, function(i) { length(YY[[i]]) } )
  XXX = YYY = ZZZ = NULL
  if (small) {
    for (s in 1:SS) {
      Nss = floor(NS[s]/4)
      XXX = c(XXX,XX[[s]][1:Nss])
      YYY = c(YYY,YY[[s]][1:Nss])
      ZZZ = c(ZZZ,rep(s,Nss))
    }}
  else{
    for (s in 1:SS) {
      XXX = c(XXX,XX[[s]])
      YYY = c(YYY,YY[[s]])
      ZZZ = c(ZZZ,rep(s,NS[s]))
    }}
  BIO = list (XX=XXX,YY=as.vector(YYY),ZZ=ZZZ,SS=SS,N=length(XXX))
  # BIO = list (XX=XXX,YY=as.vector(YYY),ZZ=ZZZ,SS=SS,N=length(XXX),genename=probeset)
return(BIO)}
