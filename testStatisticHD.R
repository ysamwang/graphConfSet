runInd <- 1
args <- commandArgs(TRUE)
for(i in 1:length(args)){
  eval(parse(text = args[[i]]))
}
print(runInd)


source("goodnessOfFit.R")

run.onceHD <- function(seed, n, p, distro, error.sd = 1, corrParam = .8, s = 3){
  K <- 3
  rec <- matrix(0, nrow = 1, ncol = 6)
  
  beta <-  c(rnorm(s), rep(0, p-s))
  if(distro == "gamma"){
    
    ## Independent Gamma's
    X <- matrix(rgamma(n * p, 1, 1) - 1, nrow = n, ncol = p)
    Y <- X %*%beta + (rgamma(n, 1, 1) - 1) * error.sd
    
  } else if(distro == "gammaCor"){
    
    ## Correlated Gamma's
    M <- toeplitz(corrParam^(0:(p-1)))
    X <- lcmix::rmvgamma(n, corr = M) - 1
    Y <- X %*%beta + (rgamma(n, 1, 1) - 1) * error.sd
    
  } else if(distro == "unif"){
    
    X <- matrix(runif(n * p, -sqrt(3), sqrt(3)), nrow = n, ncol = p)
    Y <- X %*%beta + runif(n, -sqrt(3), sqrt(3)) * error.sd
    
  } else if(distro == "unifCor"){
    
    M <- toeplitz(corrParam^(0:(p-1)))
    X <- (MultiRNG::draw.d.variate.uniform(1000, p, M) - .5) * 2  * sqrt(3)
    Y <- X %*%beta + runif(n, -sqrt(3), sqrt(3)) * error.sd
    
  } else if (distro == "laplace"){
    
    X <- matrix(rmutil::rlaplace(n * p, 0, 1/sqrt(2)), nrow = n, ncol = p)
    Y <- X %*%beta + rmutil::rlaplace(n, 0, 1/sqrt(2)) * error.sd
    
  } else if (distro == "laplaceCor"){

    M <- toeplitz(corrParam^(0:(p-1)))
    X <- LaplacesDemon::rmvl(n, mu = rep(0, p), Sigma = M)
    Y <- X %*%beta + rmutil::rlaplace(n, 0, 1/sqrt(2)) * error.sd
    
  } else if (distro == "t10") {
    
    M <- toeplitz(corrParam^(0:(p-1)))
    X <- mvtnorm::rmvt(n, sigma = M * 8 / 10, df = 10)
    Y <- X %*%beta + (rt(n, df = 10) * sqrt(8/10)) * error.sd
    
  } else if (distro == "weibull"){
    
    X <- matrix((rweibull(n * p, .5, scale = 1) - 2) / sqrt(20), nrow = n, ncol = p)
    Y <- X %*%beta + (rweibull(n, .5, scale = 1) - 2) / sqrt(20) * error.sd
    
  } else if (distro == "weibullCor"){
    
    M <- toeplitz(corrParam^(0:(p-1)))
    X <- (lcmix::rmvweisd(n, shape = .5, decay =1, corr = M) - 2) / sqrt(20)
    Y <- X %*%beta + (rweibull(n, .5, scale = 1) - 2) / sqrt(20) * error.sd
    
  }
  
  
  
  if(p > 150){
    
    ## if p > 150, then don't run HSIC
    bs <- 100
    rec[, 1] <- getPvalHD(X, Y, type = "moment", phil = "perm", K = K, bs = bs)$p
    rec[, 3] <- getPvalHD(X, Y, type = "moment", phil = "jointApprox", K = K, bs = bs)$p
    rec[, 5] <- getPvalHD(X, Y, type = "moment", phil = "residBS", K = K, bs = bs)$p
  } else {
    bs <- 100
    rec[, 1] <- getPvalHD(X, Y, type = "moment", phil = "perm", K = K, bs = bs)$p
    rec[, 2] <- getPvalHD(X, Y, type = "hsic", phil = "perm", bs = bs)$p
    rec[, 3] <- getPvalHD(X, Y, type = "moment", phil = "jointApprox", K = K, bs = bs)$p
    rec[, 4] <- getPvalHD(X, Y, type = "hsic", phil = "jointApprox", bs = bs)$p
    rec[, 5] <- getPvalHD(X, Y, type = "moment", phil = "residBS", K = K, bs = bs)$p
    rec[, 6] <- getPvalHD(X, Y, type = "hsic", phil = "residBS", bs = bs)$p
  }


  return(rec)
}



##################
library(parallel)
runInd <- 1
ss <- 10
p.list <- c(10)
esd.list <- c(1)
d.list <- c("gamma")
param.grid <- expand.grid(p.list, esd.list, d.list)


# cl <- makeCluster(31)



p <- param.grid[runInd, 1]
n <- p^2
error.sd <- param.grid[runInd, 2]
distro <- param.grid[runInd, 3]
if(p > 150){
    sim.size <- ss/2
} else {
    sim.size <- ss
}
  
# clusterExport(cl, ls())

run.onceHD(0, 100, 150, "gamma", 1, s = 3)



out <- parSapply(cl, 1:sim.size, run.onceHD, n, p, as.character(distro), error.sd)
write.csv(t(out), paste("hvmRes/res_HVM_n",n, "_p", p,"_esd", error.sd * 100, "_",distro, ".csv", sep = ""))

stopCluster(cl)


