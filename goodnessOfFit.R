
## Calculate test Statistic
getStat <- function(X, errs, type = "moment", K = 3){
  n <- dim(X)[1]
  if(type == "moment"){
    
    return(max(abs(errs %*% X^K / n)))

  } else if (type == "hsic"){
    
    return(n * dHSIC::dhsic(X, errs)$dHSIC)
  } 
}


oneResampling <- function(X, yHat, errs, type = "moment", phil = "perm", K = 3){
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
  
  ## Return the testStatistic from using yCheck
  return(getStat(X, RcppArmadillo::fastLm(X = X, y = yCheck)$res, type, K = K))
}


nonIntersectPerm <- function(p){
  out <- sample(p)
  while(any(out == 1:p)){
    out <- sample(p)
  }
  
  return(out)
}






getPval <- function(X, Y, type = "moment", phil = "perm", K = 3, bs = 1000, pValOnly = T){
  
  ## Regression with Observed data
  regOutput <- RcppArmadillo::fastLm(X = X, y = Y)
  testStat <- getStat(X, regOutput$res, type, K = K)
  
  # preComp <- t(X^K) %*% (diag(n) - X %*% solve(t(X) %*% X , t(X)))
  
  ### Estimate Null Distribution
  nullDist <- replicate(bs, oneResampling(X, regOutput$fitted, regOutput$res, type = type, phil = phil, K = K))
  
  if(pValOnly){
    return(list(pval = mean(testStat < nullDist))) 
  } else {
    return(list(pval = mean(testStat < nullDist), testStat = testStat, nullDist = nullDist))
  }
}










