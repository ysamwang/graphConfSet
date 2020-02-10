getPvalSpeed <- function(X, Y, type = "moment", phil = "perm", K = 3, bs = 1000, pValOnly = T){
  n <- dim(X)[1]
  
  ## Regression with Observed data
  regOutput <- RcppArmadillo::fastLm(X = X, y = Y)
  errs <- regOutput$res
  preComp <- t(X^K) %*% (diag(n) - X %*% solve(t(X) %*% X , t(X)))
  
  testStat <- max(abs(preComp %*% errs) / n)
  
  n <- dim(X)[1]
  oneResampSpeed <- function(phil){
    if(phil == "perm"){
      ind <- nonIntersectPerm(n)
    } else {
      ind <- sample(n, replace = T)
    }
    max(abs(preComp %*% errs[ind]) / n)
  }
  
  nullDist <- replicate(bs, oneResampSpeed(phil))
  
  if(pValOnly){
    return(list(pval = mean(testStat < nullDist))) 
  } else {
    return(list(pval = mean(testStat < nullDist), testStat = testStat, nullDist = nullDist))
  }
}


nonIntersectPerm <- function(p){
  out <- sample(p)
  while(any(out == 1:p)){
    out <- sample(p)
  }
  
  return(out)
}

testOrdering <- function(Y, topOrder, type = "moment", phil = "perm", K = 3, bs = 1000, aggMethod = "fisher") {
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  
  ## reordering Y
  Y <- Y[, topOrder, drop = F]
  
  if(aggMethod == "fisher"){
    testStat <- 0
  } else {
    testStat <- 1
  }
  
  
  
  
  j <- 1
  while(j < p){
    
    j <- j + 1
    
    ### Test ancestral relationships
    u <- getPvalSpeed(Y[, 1:(j-1), drop = F], Y[, j, drop = F], type = type, phil = phil, K = K, bs = bs, pValOnly = T)$p
    
    if(aggMethod == "fisher"){
      testStat <- testStat - 2 * log(max(u, 10e-5))
      pValUB <- pchisq(testStat, df = 2*(p-1), lower.tail = F)
    } else {
      testStat <- min(testStat, u)
      pValUB <- mean(testStat > minOfUnif(p-1))
    }
    
  }
  
  if(aggMethod == "fisher"){
    pVal <- pchisq(testStat, df = 2*(p-1), lower.tail = F)
  } else {
    pVal <- mean(testStat > minOfUnif(p-1))
  }


  
  return(list(stat = testStat, pVal = pVal))
  
}


minOfUnif <-function(numRv, draws = 1000){
  return(apply(matrix(runif(draws * numRv), nrow = draws), MAR = 1, min))
} 



exhaustiveBB <- function(Y, phi = "perm", K = 3, bs = 1000, alpha = .1, verbose = F, aggMethod = "fisher"){
  n <- dim(Y)[1]
  p <- dim(Y)[2]
  
  possibleOrders <- matrix(1:p, nrow = p)
  if(aggMethod == "fisher"){
    testStat <- rep(0, p)
  } else {
    testStat <- rep(1, p)
    refDist <- minOfUnif(p-1, draws = 10000)
  }
  
   
  j <- 1
  while(j < p){
    if(verbose == T){
      print(paste("j =", j))
      print("Possible Orders:")
      print(possibleOrders)
    }
    
    j <- j + 1
    
    newPossibleOrders <- matrix(0, nrow = 0, ncol = j)
    newTestStat <- c()
    
    
    ## Go through all possible orders and see if you can add one more
    for(z in 1:dim(possibleOrders)[1]){
      
      X <- Y[, possibleOrders[z, ], drop = F]
      if(verbose){
        cat("\n")
        cat("=================")
        cat("Testing: ")
        cat(possibleOrders[z, ])
        cat("\n")
      }
      
      ## nodes that could be next in ordering
      possibleRoots <- setdiff(1:p, possibleOrders[z, ])
      
      for(m in 1:length(possibleRoots)){
        u <- getPvalSpeed(X, Y[, possibleRoots[m]], phil = phil, K = K, bs = bs)$p
        
        
        
        
        if(aggMethod == "fisher"){
          singleStat <-  testStat[z] - 2 * log(max(u, 10e-10))
          singleP <- pchisq(singleStat, df = 2 * (p-1), lower = F)
        } else {
          singleStat <-  min(testStat[z], u)
          singleP <- mean(singleStat >= refDist)
        }
        
        if(verbose){
          cat(paste(possibleRoots[m], " : ", singleStat, sep = ""))
          cat("\n")
        }
        
        
        
        
        ## If passes cutoff, add to next round
        if (singleP > alpha){
          newPossibleOrders <- rbind(newPossibleOrders, c(possibleOrders[z, ], possibleRoots[m]))
          newTestStat[length(newTestStat) + 1] <- singleStat
        }
        
      } # end search through possible roots
    
    } # end search through currently possible orderings
    
    possibleOrders <- newPossibleOrders
    testStat <- newTestStat
    
  } # end while loop
  
  
  if(aggMethod == "fisher"){
    pVals <- sapply(testStat, FUN = function(x){pchisq(x, df = 2 * (p-1), lower = F)})
  } else {
    pVals <- sapply(testStat, FUN = function(x){mean(x >= refDist)})
  }
    return(list(orders = possibleOrders, pVals = pVals))
}

