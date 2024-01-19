
# Helper function to quickly calculate the AUC without worrying about the ROC curve
fast_auc = function(x, y) {
  rank_value = rank(x, ties.method="average")
  n_pos = sum(y == 1)
  rank_sum = sum(rank_value[y == 1])
  u_value = rank_sum - (n_pos * (n_pos + 1)) / 2
  auc = u_value / (n_pos * (length(x) - n_pos))
  if (auc >= 0.50) auc else 1.0 - auc 
}

# Helper function for fast t-test
fast_t_test = function(x, y) {
  nx = length(x); mx = mean(x); vx = mean((x - mx)^2)/(nx-1);
  ny = length(y); my = mean(y); vy = mean((y - my)^2)/(ny-1);
  stderrx <- sqrt(vx); stderry <- sqrt(vy);
  stderr <- sqrt(stderrx^2 + stderry^2)
  df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny -  1))
  tstat = (mx - my)/stderr
  common.sd = sqrt( ( sum((x - mx)^2) + sum((y - my)^2) ) / (nx+ny-2) )
  list(p_value = 2 * pt(-abs(tstat), df),
       pooled_se = stderr, common_sd = common.sd, tstat = tstat )
}

# Calculate summary statistics for each gene in the expression matrix
# AUC = AUC
# nlpvalueT = log p-value for t test
# FoldChange = difference in log means
# ResiErr = pooled residual error estimate from t test
# 
# Use it as `scores = CompScores(XX, YY)`
CompScores = function(XX, YY) {
  NGenes = nrow(XX)
  # create empty vectors to store results in
  AUC = numeric(NGenes)
  nlpvalueT = numeric(NGenes)
  resi = numeric(NGenes)
  
  for (gg in 1:NGenes) { #gg is index
    # AUC
    AUC[gg] = fast_auc(XX[gg, ], YY)
    # Tests difference in means for XX under Y=1 and XX under Y=0
    # T-test negative log p-value
    t_test = fast_t_test(XX[gg, YY==1], XX[gg, YY==0])
    nlpvalueT[gg] = -log(t_test$p_value)
    resi[gg] = t_test$common_sd
  }
  
  # Difference in means (aka fold change as data are on log scale)
  FoldChange = rowMeans(XX[, YY==1]) - rowMeans(XX[, YY==0])
  
  scores = data.frame(AUC=AUC, nlpvalueT=nlpvalueT, 
                      FoldChange=FoldChange, ResiErr=resi)
  
  return(scores)
}

