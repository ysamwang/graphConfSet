p <- 8
n <- 500
dat <- rDAG(p, n, maxInDegree = p, dist = "gamma", lowScale = .5, highScale = 1, lowEdge = .8, highEdge = 1)
Y <- dat$Y
cov(Y)


system.time(testOrdering(Y, 1:p, type = "moment", phil = "perm", K = 3, bs = 1000, aggMethod = "fisher"))
system.time(testOrdering(Y, 1:p, type = "moment", phil = "perm", K = 3, bs = 1000, aggMethod = "other"))



system.time(outFisher <- exhaustiveBB(Y, verbose = F, alpha = .1, aggMethod = "fisher"))
system.time(outOther <- exhaustiveBB(Y, verbose = F, alpha = .1, aggMethod = "other"))




dim(outFisher$orders)[1] / factorial(p)
dim(outOther$orders)[1] / factorial(p)

