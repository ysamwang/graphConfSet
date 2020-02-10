Y <- as.matrix(read.csv("data/Bollentestdata.txt", sep = "\t", header = T))
Y <- scale(Y)

par(mfrow = c(2,3))
for(i in 1:6){
hist(Y[,i])
}


outFisher <- exhaustiveBB(Y, verbose = T, aggMethod = "fisher", K = 3, alpha = .1)
outOther <- exhaustiveBB(Y, verbose = T, aggMethod = "other", K = 3, alpha = .1)



getPvalSpeed(X = Y[, c(1,3,6)], Y = Y[, 5])
getPvalSpeed(X = Y[, c(1,3,6, 5)], Y = Y[, 4])
getPvalSpeed(X = Y[, c(1,3,6, 5, 4)], Y = Y[, 2])


pchisq(-2 * (log(.023) + log(.751) + log(.098)), df = 4, lower = F)



quantile(minOfUnif(3), .1)
