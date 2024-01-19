masomenos.train = function(X,   # pxn data frame of predictors
                           Y,   # nx1 vectors of labels or responses
                           P,   # number of predictors
                           training.criterion="AUC",
                           filtering.fraction=.5)
{
  # eliminate unclassified samples
  YY = Y[!is.na(Y)]
  XX = X[,!is.na(Y)]
  
  Nvar = nrow(XX) # Number of genes
  crite = rep(NA,Nvar) # These are to store the criterion for each gene.
  if (training.criterion=="AUC"){
    for (pp in 1:Nvar){
      crite[pp] <- as.numeric( wilcox.test(XX[pp,]~YY)$statistic / (sum(YY==0)*sum(YY==1)) )
      #XX[pp,] gives data values, YY gives levels of the groups to be compared.
    }
  }
  cutoff = sort(abs(crite),decreasing = TRUE)[P]
  # sort the absolute value of the criterion from largest to smallest, then look at the corresponding cutoff value of the Pth largest value.
  cutoff
  variables = (1:Nvar)[abs(crite) >= cutoff]
  # retrieve the variable name/number of the top P.
  variables
  variables.signs = ( 2 * ( crite > 0.5 ) - 1 ) [variables]
  scores = apply ( XX[ variables, ] * variables.signs, 2, mean )
  # this is known as the risk score, or the mean expression value for the top P genes for each patient.

  if (training.criterion=="AUC"){
    crite.mom = as.numeric( wilcox.test(scores~YY)$statistic / (sum(YY==0)*sum(YY==1)) )
  }
  
  MoM = list(XX=XX,YY=YY,cutoff=cutoff,
             training.criterion=training.criterion,
             variables = variables,
             variables.signs=variables.signs,
             variables.criterion=crite[variables],
             scores=scores,
             criterion.mom=crite.mom)
  return(MoM)
}

masomenos.test = function(X,   # pxn data frame of predictors
                          Y,   # nx1 vectors of labels or responses
                          MoM.out # output form masomenos.train
){
  # eliminate unclassified samples
  YY = Y[!is.na(Y)]
  XX = X[,!is.na(Y)]

  Nvar = nrow(XX)
  scores = apply ( XX[ MoM.out$variables, ] * MoM.out$variables.signs, 2, mean )
  
  if (MoM.out$training.criterion=="AUC"){ 
    crite.mom = as.numeric( wilcox.test(scores~YY)$statistic / (sum(YY==0)*sum(YY==1)) ) 
  }
  
  MoM = list(XX=XX,YY=YY,
             scores=scores,
             criterion.mom=crite.mom)
  return(MoM)
}


