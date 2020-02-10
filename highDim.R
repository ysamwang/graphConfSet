########################################################
source("debiasedLasso.r")

oneResamplingHD <- function(X, yHat, errs, type = "moment", phil = "perm", K = 3){
  n <- dim(X)[1]
  p <- dim(X)[2]
  if(phil == "jointApprox"){
    
    ## Joint bootstrap approximation
    errs <- errs - mean(errs)
    ind <- sample(n, replace = T)
    X <- X[ind, , drop = F]
    yCheck <- yHat[ind] + sample(errs, replace = T)
    
  } else if(phil == "perm"){
    
    ## Residual randomization
    yCheck <- yHat + errs[nonIntersectPerm(n)]
    
  } else if (phil == "residBS") {
    
    ## Residual Bootstrap approximation of just error distribution
    errs <- errs - mean(errs)
    yCheck <- yHat + errs[sample(n, replace = T)]
  }
  
  refit <- SSLasso(X, yCheck)  
  ## Return the testStatistic from using yCheck
  return(getStat(X, refit$res, type, K = K))
}


nonIntersectPerm <- function(p){
  out <- sample(p)
  while(any(out == 1:p)){
    out <- sample(p)
  }
  
  return(out)
}



getPvalHD <- function(X, Y, type = "moment", phil = "perm", K = 3, bs = 1000, pValOnly = T){
  
  ## Regression with Observed data
  ### UPDATE ###
  regOutput <- SSLasso(X, Y)
  testStat <- getStat(X, regOutput$res, type, K = K)
  
  ### Estimate Null Distribution
  nullDist <- replicate(bs, oneResamplingHD(X, regOutput$fitted, regOutput$res, type = type, phil = phil, K = K))
  
  if(pValOnly){
    return(list(pval = mean(testStat < nullDist))) 
  } else {
    return(list(pval = mean(testStat < nullDist), testStat = testStat, nullDist = nullDist))
  }
}






